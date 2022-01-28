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
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import logging
import math
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
from sourcespec.ssp_radiated_energy import radiated_energy
from sourcespec.ssp_inversion_types import InitialValues, Bounds
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
    if config.inv_algorithm == 'TNC':
        minimize_func = objective_func(freq_log, ydata, weight)
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
        params_err = np.sqrt(params_cov.diagonal())
        # symmetric error
        params_err = ((e, e) for e in err)
    elif config.inv_algorithm == 'BH':
        minimize_func = objective_func(freq_log, ydata, weight)
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
        params_err = np.sqrt(params_cov.diagonal())
        # symmetric error
        params_err = ((e, e) for e in err)
    elif config.inv_algorithm in ['GS', 'IS']:
        minimize_func = objective_func(freq_log, ydata, weight)
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
    return params_opt, params_err


def _spec_inversion(config, spec, noise_weight):
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

    noise_weight = noise_weight.data_log
    if config.weighting == 'noise':
        weight = noise_weight
        # 'curve_fit' interprets 'yerr' as standard deviation vector
        # and calculates weights as 1/yerr^2 .
        # Therefore we build yerr as:
        yerr = 1./np.sqrt(weight)
    elif config.weighting == 'frequency':
        # frequency weighting:
        #   config.weight for f<=f_weight
        #   1      for f> f_weight
        yerr = np.ones(len(ydata))
        yerr[freq_log <= config.f_weight] = 1./math.sqrt(config.weight)
        weight = 1./np.power(yerr, 2)
    elif config.weighting is None:
        weight = yerr = np.ones(len(ydata))

    # Find the frequency range to compute Mw_0 and, possibly, t_star_0:
    # we start where signal-to-noise becomes strong
    idx0 = np.where(noise_weight > 0.5)[0][0]
    # we stop at the first max of signal-to-noise (proxy for fc)
    idx_max = argrelmax(noise_weight)[0]
    # just keep the indexes for maxima > 0.5
    idx_max = [idx for idx in idx_max if noise_weight[idx] > 0.5]
    if not idx_max:
        # if idx_max is empty, then the source and/or noise spectrum
        # is most certainly "strange". In this case, we simply give up.
        logger.warning(
            '{}: unable to find a frequency range to compute Mw_0'.format(
                statId))
        logger.warning(
            '   This is possibly due to an uncommon '
            'spectrum for the trace (e.g., a resonance).')
        raise RuntimeError
    idx1 = idx_max[0]
    if idx1 == idx0 and len(idx_max) > 1:
        idx1 = idx_max[1]
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
        params_opt, params_err = _curve_fit(
            config, spec, weight, yerr, initial_values, bounds)
    except (RuntimeError, ValueError) as m:
        logger.warning(m)
        logger.warning('{}: unable to fit spectral model'.format(statId))
        raise

    params_name = ('Mw', 'fc', 't_star')
    par = OrderedDict(zip(params_name, params_opt))
    par_str = '; '.join(['{}: {:.4f}'.format(key, par[key]) for key in par])
    logger.info('{}: optimal values: {}'.format(statId, par_str))

    # Ignore spectra with negative fc or t_star
    t_star = par['t_star']
    if t_star < 0:
        logger.warning(
            '{}: t_star: {:.3f} < 0: ignoring inversion results'.format(
                statId, t_star))
        raise ValueError
    fc = par['fc']
    if fc < 0:
        logger.warning(
            '{}: fc: {:.3f} < 0: ignoring inversion results'.format(
                statId, fc))
        raise ValueError

    par['Mo'] = mag_to_moment(par['Mw'])
    par['hyp_dist'] = spec.stats.hypo_dist
    par['epi_dist'] = spec.stats.epi_dist
    par['az'] = az
    par['lon'] = spec.stats.coords.longitude
    par['lat'] = spec.stats.coords.latitude

    # additional parameters, computed from fc, Mw and t_star
    vs = config.hypo.vs
    # source radius in meters
    par['ra'] = source_radius(par['fc'], vs*1e3)
    # Brune stress drop in MPa
    par['bsd'] = bsd(par['Mo'], par['ra'])
    # quality factor
    par['Qo'] = quality_factor(par['hyp_dist'], vs, par['t_star'])

    # Check if Brune stress drop is acceptable
    if config.bsd_min_max is not None:
        bsd_min, bsd_max = config.bsd_min_max
        if not (bsd_min <= par['bsd'] <= bsd_max):
            msg = '{}: bsd: {:.3e} not in allowed range [{:.3e}, {:.3e}]: '
            msg += 'ignoring inversion results'
            msg = msg.format(statId, par['bsd'], bsd_min, bsd_max)
            logger.warning(msg)
            raise ValueError

    par_err = OrderedDict(zip(params_name, params_err))
    return par, par_err


def _synth_spec(config, spec, par, par_err):
    """Return a stream with one or more synthetic spectra."""
    spec_st = Stream()
    params_opt = [par[key] for key in ('Mw', 'fc', 't_star')]

    freq = spec.get_freq()
    spec_synth = spec.copy()
    spec_synth.stats.channel = spec.stats.channel[:-1] + 'S'
    spec_synth.stats.par = par
    spec_synth.stats.par_err = par_err
    spec_synth.data_mag = spectral_model(freq, *params_opt)
    spec_synth.data = mag_to_moment(spec_synth.data_mag)
    spec_st.append(spec_synth)

    # Add an extra spectrum with no attenuation
    if config.plot_spectra_no_attenuation:
        spec_synth = spec.copy()
        spec_synth.stats.channel = spec.stats.channel[:-1] + 's'
        _params = list(params_opt)
        _params[-1] = 0
        spec_synth.data_mag = spectral_model(freq, *_params)
        spec_synth.data = mag_to_moment(spec_synth.data_mag)
        spec_st.append(spec_synth)

    if config.plot_spectra_no_fc:
        spec_synth = spec.copy()
        spec_synth.stats.channel = spec.stats.channel[:-1] + 't'
        _params = list(params_opt)
        _params[1] = 1e999
        spec_synth.data_mag = spectral_model(freq, *_params)
        spec_synth.data = mag_to_moment(spec_synth.data_mag)
        spec_st.append(spec_synth)
    return spec_st


def spectral_inversion(config, spec_st, weight_st):
    """Inversion of displacement spectra."""
    logger.info('Inverting spectra...')
    weighting_messages = {
        'noise': 'Using noise weighting for inversion.',
        'frequency': 'Using frequency weighting for inversion.',
        None: 'Using no weighting for inversion.'
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

    sourcepar = OrderedDict()
    sourcepar_err = OrderedDict()
    stations = set(x.stats.station for x in spec_st)
    spectra = [sp for sta in stations for sp in spec_st.select(station=sta)]
    for spec in sorted(spectra, key=lambda sp: sp.id):
        if spec.stats.channel[-1] != 'H':
            continue
        if spec.stats.ignore:
            continue
        noise_weight = select_trace(weight_st, spec.id, spec.stats.instrtype)
        try:
            par, par_err = _spec_inversion(config, spec, noise_weight)
        except (RuntimeError, ValueError):
            continue
        spec_st += _synth_spec(config, spec, par, par_err)
        statId = '{} {}'.format(spec.id, spec.stats.instrtype)
        sourcepar[statId] = par
        sourcepar_err[statId] = par_err

    radiated_energy(config, spec_st, sourcepar)

    logger.info('Inverting spectra: done')
    return sourcepar, sourcepar_err
