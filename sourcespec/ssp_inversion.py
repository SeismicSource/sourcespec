# -*- coding: utf8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Spectral inversion routines for sourcespec.

:copyright:
    2012 Claudio Satriano <satriano@ipgp.fr>

    2013-2014 Claudio Satriano <satriano@ipgp.fr>,
              Emanuela Matrullo <matrullo@geologie.ens.fr>,
              Agnes Chounet <chounet@ipgp.fr>

    2015-2022 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
import logging
import numpy as np
from collections import OrderedDict
from scipy.optimize import curve_fit, minimize, basinhopping
from scipy.signal import argrelmax
from obspy import Stream
from obspy.geodetics import gps2dist_azimuth
from sourcespec.ssp_spectral_model import (
    spectral_model, objective_func, callback)
from sourcespec.ssp_util import (
    mag_to_moment, source_radius, bsd, quality_factor, select_trace, smooth)
from sourcespec.ssp_data_types import (
    InitialValues, Bounds, StationSourceParameters, SourceParameters)
from sourcespec.ssp_grid_sampling import GridSampling
logger = logging.getLogger(__name__.split('.')[-1])


def _curve_fit(config, spec, weight, yerr, initial_values, bounds):
    """
    Curve fitting.

    Available algorithms:
      - Levenberg-Marquardt (LM, via `curve_fit()`). Automatically switches to
        Trust Region Reflective algorithm if bounds are provided.
      - Truncated Newton algorithm (TNC) with bounds.
      - Basin-hopping (BH)
      - Grid search (GS)
    """
    freq_log = spec.freq_log
    ydata = spec.data_log_mag
    minimize_func = objective_func(freq_log, ydata, weight)
    if config.inv_algorithm == 'TNC':
        res = minimize(
            minimize_func,
            x0=initial_values.get_params0(), method='TNC',
            callback=callback, bounds=bounds.bounds
        )
        params_opt = res.x
        # trick: use curve_fit() bounded to params_opt
        # to get the covariance
        _, params_cov = curve_fit(
            spectral_model, freq_log, ydata,
            p0=params_opt, sigma=yerr,
            bounds=(params_opt-(1e-10), params_opt+(1e-10))
        )
        err = np.sqrt(params_cov.diagonal())
        # symmetric error
        params_err = ((e, e) for e in err)
    elif config.inv_algorithm == 'LM':
        bnds = bounds.get_bounds_curve_fit()
        if bnds is not None:
            logger.info(
                'Trying to use using Levenberg-Marquardt '
                'algorithm with bounds. Switching to the '
                'Trust Region Reflective algorithm.'
            )
        params_opt, params_cov = curve_fit(
            spectral_model, freq_log, ydata,
            p0=initial_values.get_params0(), sigma=yerr,
            bounds=bnds
        )
        err = np.sqrt(params_cov.diagonal())
        # symmetric error
        params_err = ((e, e) for e in err)
    elif config.inv_algorithm == 'BH':
        res = basinhopping(
            minimize_func, x0=initial_values.get_params0(), niter=100,
            accept_test=bounds
        )
        params_opt = res.x
        # trick: use curve_fit() bounded to params_opt
        # to get the covariance
        _, params_cov = curve_fit(
            spectral_model, freq_log, ydata,
            p0=params_opt, sigma=yerr,
            bounds=(params_opt-(1e-10), params_opt+(1e-10))
        )
        err = np.sqrt(params_cov.diagonal())
        # symmetric error
        params_err = ((e, e) for e in err)
    elif config.inv_algorithm in ['GS', 'IS']:
        nsteps = (20, 150, 150)  # we do fewer steps in magnitude
        sampling_mode = ('lin', 'log', 'lin')
        params_name = ('Mw', 'fc', 't_star')
        params_unit = ('', 'Hz', 's')
        grid_sampling = GridSampling(
            minimize_func, bounds.bounds, nsteps,
            sampling_mode, params_name, params_unit)
        if config.inv_algorithm == 'GS':
            grid_sampling.grid_search()
        elif config.inv_algorithm == 'IS':
            grid_sampling.kdtree_search()
        params_opt = grid_sampling.params_opt
        params_err = grid_sampling.params_err
        spec_label = '{} {}'.format(spec.id, spec.stats.instrtype)
        grid_sampling.plot_conditional_misfit(config, spec_label)
        # fc-t_star
        plot_par_idx = (1, 2)
        grid_sampling.plot_misfit_2d(config, plot_par_idx, spec_label)
        # fc-Mw
        plot_par_idx = (1, 0)
        grid_sampling.plot_misfit_2d(config, plot_par_idx, spec_label)
        # tstar-Mw
        plot_par_idx = (2, 0)
        grid_sampling.plot_misfit_2d(config, plot_par_idx, spec_label)
    misfit = minimize_func(params_opt)
    return params_opt, params_err, misfit


def _freq_ranges_for_Mw0_and_tstar0(config, weight, freq_log, statId):
    """Find the frequency range to compute Mw_0 and, possibly, t_star_0."""
    if config.weighting == 'noise':
        # we start where signal-to-noise becomes strong
        idx0 = np.where(weight > 0.5)[0][0]
        # we stop at the first max of signal-to-noise (proxy for fc)
        idx_max = argrelmax(weight)[0]
        # just keep the indexes for maxima > 0.5
        idx_max = [idx for idx in idx_max if weight[idx] > 0.5]
        if not idx_max:
            # if idx_max is empty, then the source and/or noise spectrum
            # is most certainly "strange". In this case, we simply give up.
            msg = '{}: unable to find a frequency range to compute Mw_0. '
            msg += 'This is possibly due to an uncommon spectrum '
            msg += '(e.g., a resonance).'
            msg = msg.format(statId)
            raise RuntimeError(msg)
        idx1 = idx_max[0]
        if idx1 == idx0:
            try:
                idx1 = idx_max[1]
            except IndexError:
                # if there are no other maxima, just take 5 points
                idx1 = idx0+5
    elif config.weighting == 'frequency':
        idx0 = 0
        # the closest index to f_weight:
        idx1 = np.where(freq_log <= config.f_weight)[0][-1]
    else:
        idx0 = 0
        idx1 = int(len(weight)/2)
    return idx0, idx1


def _spec_inversion(config, spec, spec_weight):
    """Invert one spectrum."""
    # azimuth computation
    coords = spec.stats.coords
    hypo = spec.stats.hypo
    stla = coords.latitude
    stlo = coords.longitude
    evla = hypo.latitude
    evlo = hypo.longitude
    geod = gps2dist_azimuth(evla, evlo, stla, stlo)
    az = geod[1]

    freq_log = spec.freq_log
    ydata = spec.data_log_mag
    statId = '{} {}'.format(spec.id, spec.stats.instrtype)
    weight = spec_weight.data_log

    # 'curve_fit' interprets 'yerr' as standard deviation array
    # and calculates weights as 1/yerr^2 .
    # Therefore we build yerr as:
    yerr = 1./np.sqrt(weight)

    # Find frequency range (indexes) to compute Mw_0 and t_star_0
    # When using noise weighting, idx1 is the first maximum in
    # signal-to-noise ratio
    idx0, idx1 = _freq_ranges_for_Mw0_and_tstar0(
        config, weight, freq_log, statId)

    # first maximum is a proxy for fc, we use it for fc_0:
    fc_0 = freq_log[idx1]

    t_star_min = t_star_max = None
    if config.invert_t_star_0:
        # fit t_star_0 and Mw on the initial part of the spectrum,
        # corrected for the effect of fc
        ydata_corr = ydata - spectral_model(freq_log, Mw=0, fc=fc_0, t_star=0)
        ydata_corr = smooth(ydata_corr, window_len=18)
        slope, Mw_0 = np.polyfit(
            freq_log[idx0: idx1], ydata_corr[idx0: idx1], deg=1)
        t_star_0 = -3./2 * slope / (np.pi * np.log10(np.e))
        t_star_min = t_star_0 * (1 - config.t_star_0_variability)
        t_star_max = t_star_0 * (1 + config.t_star_0_variability)
    if not config.invert_t_star_0 or t_star_0 < 0:
        # we calculate the initial value for Mw as an average
        Mw_0 = np.nanmean(ydata[idx0: idx1])
        t_star_0 = config.t_star_0

    initial_values = InitialValues(Mw_0, fc_0, t_star_0)
    logger.info('{}: initial values: {}'.format(statId, str(initial_values)))
    bounds = Bounds(config, spec, initial_values)
    bounds.Mw_min = Mw_0 - config.Mw_0_variability
    bounds.Mw_max = Mw_0 + config.Mw_0_variability
    if t_star_min is not None:
        bounds.t_star_min = t_star_min
    if t_star_max is not None:
        bounds.t_star_max = t_star_max
    logger.info('{}: bounds: {}'.format(statId, str(bounds)))
    try:
        params_opt, params_err, misfit = _curve_fit(
            config, spec, weight, yerr, initial_values, bounds)
    except (RuntimeError, ValueError) as m:
        msg = str(m) + '\n'
        msg += '{}: unable to fit spectral model'.format(statId)
        raise RuntimeError(msg)

    params_name = ('Mw', 'fc', 't_star')
    par = OrderedDict(zip(params_name, params_opt))
    par_str = '; '.join(['{}: {:.4f}'.format(key, par[key]) for key in par])
    logger.info('{}: optimal values: {}'.format(statId, par_str))
    logger.info('{}: misfit: {:.3f}'.format(statId, misfit))

    if np.isclose(par['fc'], bounds.fc_min, rtol=0.1):
        msg = '{}: optimal fc within 10% of fc_min: {:.3f} ~= {:.3f}: '
        msg += 'ignoring inversion results'
        msg = msg.format(statId, par['fc'], bounds.fc_min)
        raise ValueError(msg)

    if np.isclose(par['fc'], bounds.fc_max, rtol=1e-4):
        msg = '{}: optimal fc within 10% of fc_max: {:.3f} ~= {:.3f}: '
        msg += 'ignoring inversion results'
        msg = msg.format(statId, par['fc'], bounds.fc_max)
        raise ValueError(msg)

    misfit_max = config.pi_misfit_max or np.inf
    if misfit > misfit_max:
        msg = '{}: misfit larger than pi_misfit_max: {:.3f} > {:.3f}: '
        msg += 'ignoring inversion results'
        msg = msg.format(statId, misfit, misfit_max)
        raise ValueError(msg)

    # Check post-inversion bounds for t_star and fc
    t_star = par['t_star']
    pi_t_star_min, pi_t_star_max =\
        config.pi_t_star_min_max or (-np.inf, np.inf)
    if not (pi_t_star_min <= t_star <= pi_t_star_max):
        msg = '{}: t_star: {:.3f} not in allowed range [{:.3f}, {:.3f}]: '
        msg += 'ignoring inversion results'
        msg = msg.format(statId, t_star, pi_t_star_min, pi_t_star_max)
        raise ValueError(msg)
    fc = par['fc']
    pi_fc_min, pi_fc_max = config.pi_fc_min_max or (-np.inf, np.inf)
    if not (pi_fc_min <= fc <= pi_fc_max):
        msg = '{}: fc: {:.3f} not in allowed range [{:.3f}, {:.3f}]: '
        msg += 'ignoring inversion results'
        msg = msg.format(statId, fc, pi_fc_min, pi_fc_max)
        raise ValueError(msg)

    par['hyp_dist'] = spec.stats.hypo_dist
    par['epi_dist'] = spec.stats.epi_dist
    par['az'] = az
    par['lon'] = spec.stats.coords.longitude
    par['lat'] = spec.stats.coords.latitude

    # additional parameters, computed from fc, Mw and t_star
    vs = config.vs_source
    travel_time = spec.stats.travel_times[config.wave_type[0]]
    # seismic moment
    par['Mo'] = mag_to_moment(par['Mw'])
    # source radius in meters
    par['ra'] = source_radius(par['fc'], vs*1e3)
    # Brune stress drop in MPa
    par['bsd'] = bsd(par['Mo'], par['ra'])
    # quality factor
    par['Qo'] = quality_factor(travel_time, par['t_star'])

    # Check post-inversion bounds for bsd
    pi_bsd_min, pi_bsd_max = config.pi_bsd_min_max or (-np.inf, np.inf)
    if not (pi_bsd_min <= par['bsd'] <= pi_bsd_max):
        msg = '{}: bsd: {:.3e} not in allowed range [{:.3e}, {:.3e}]: '
        msg += 'ignoring inversion results'
        msg = msg.format(statId, par['bsd'], pi_bsd_min, pi_bsd_max)
        raise ValueError(msg)

    # additional parameter errors, computed from fc, Mw and t_star
    par_err = OrderedDict(zip(params_name, params_err))
    # seismic moment
    Mw_min = par['Mw'] - par_err['Mw'][0]
    Mw_max = par['Mw'] + par_err['Mw'][1]
    Mo_min = mag_to_moment(Mw_min)
    Mo_max = mag_to_moment(Mw_max)
    par_err['Mo'] = (par['Mo'] - Mo_min, Mo_max - par['Mo'])
    # source radius in meters
    fc_min = par['fc'] - par_err['fc'][0]
    if fc_min <= 0:
        fc_min = freq_log[0]
    fc_max = par['fc'] + par_err['fc'][1]
    ra_min = source_radius(fc_max, vs*1e3)
    ra_max = source_radius(fc_min, vs*1e3)
    par_err['ra'] = (par['ra']-ra_min, ra_max-par['ra'])
    # Brune stress drop in MPa
    bsd_min = bsd(Mo_min, ra_max)
    bsd_max = bsd(Mo_max, ra_min)
    par_err['bsd'] = (par['bsd']-bsd_min, bsd_max-par['bsd'])
    # quality factor
    t_star_min = par['t_star'] - par_err['t_star'][0]
    if t_star_min <= 0:
        t_star_min = 0.001
    t_star_max = par['t_star'] + par_err['t_star'][1]
    Qo_min = quality_factor(travel_time, t_star_max)
    Qo_max = quality_factor(travel_time, t_star_min)
    par_err['Qo'] = (par['Qo']-Qo_min, Qo_max-par['Qo'])

    return par, par_err


def _synth_spec(config, spec, par, par_err):
    """Return a stream with one or more synthetic spectra."""
    spec_st = Stream()
    params_opt = [par[key] for key in ('Mw', 'fc', 't_star')]

    freq = spec.get_freq()
    freq_log = spec.freq_log
    spec_synth = spec.copy()
    spec_synth.stats.channel = spec.stats.channel[:-1] + 'S'
    spec_synth.stats.par = par
    spec_synth.stats.par_err = par_err
    spec_synth.data_mag = spectral_model(freq, *params_opt)
    spec_synth.data = mag_to_moment(spec_synth.data_mag)
    spec_synth.data_log_mag = spectral_model(freq_log, *params_opt)
    spec_synth.data_log = mag_to_moment(spec_synth.data_log_mag)
    spec_st.append(spec_synth)

    # Add an extra spectrum with no attenuation
    if config.plot_spectra_no_attenuation:
        spec_synth = spec.copy()
        spec_synth.stats.channel = spec.stats.channel[:-1] + 's'
        _params = list(params_opt)
        _params[-1] = 0
        spec_synth.data_mag = spectral_model(freq, *_params)
        spec_synth.data = mag_to_moment(spec_synth.data_mag)
        spec_synth.data_log_mag = spectral_model(freq_log, *_params)
        spec_synth.data_log = mag_to_moment(spec_synth.data_log_mag)
        spec_st.append(spec_synth)

    if config.plot_spectra_no_fc:
        spec_synth = spec.copy()
        spec_synth.stats.channel = spec.stats.channel[:-1] + 't'
        _params = list(params_opt)
        _params[1] = 1e999
        spec_synth.data_mag = spectral_model(freq, *_params)
        spec_synth.data = mag_to_moment(spec_synth.data_mag)
        spec_synth.data_log_mag = spectral_model(freq_log, *_params)
        spec_synth.data_log = mag_to_moment(spec_synth.data_log_mag)
        spec_st.append(spec_synth)
    return spec_st


def spectral_inversion(config, spec_st, weight_st):
    """Inversion of displacement spectra."""
    logger.info('Inverting spectra...')
    weighting_messages = {
        'noise': 'Using noise weighting for inversion.',
        'frequency': 'Using frequency weighting for inversion.',
        'no_weight': 'Using no weighting for inversion.'
    }
    logger.info(weighting_messages[config.weighting])
    algorithm_messages = {
        'TNC': 'Using truncated Newton algorithm for inversion.',
        'LM': 'Using Levenberg-Marquardt algorithm for inversion.',
        'BH': 'Using basin-hopping algorithm for inversion.',
        'GS': 'Using grid search for inversion.',
        'IS': 'Using k-d tree importance sampling for inversion.'
    }
    logger.info(algorithm_messages[config.inv_algorithm])

    sourcepar = SourceParameters()
    stations = set(x.stats.station for x in spec_st)
    spectra = [sp for sta in stations for sp in spec_st.select(station=sta)]
    for spec in sorted(spectra, key=lambda sp: sp.id):
        if spec.stats.channel[-1] != 'H':
            continue
        if spec.stats.ignore:
            continue
        spec_weight = select_trace(weight_st, spec.id, spec.stats.instrtype)
        try:
            par, par_err = _spec_inversion(config, spec, spec_weight)
        except (RuntimeError, ValueError) as msg:
            logger.warning(msg)
            continue
        spec_st += _synth_spec(config, spec, par, par_err)
        statId = '{} {}'.format(spec.id, spec.stats.instrtype)
        par = StationSourceParameters(statId, par, par_err)
        sourcepar.station_parameters[statId] = par

    logger.info('Inverting spectra: done')
    return sourcepar
