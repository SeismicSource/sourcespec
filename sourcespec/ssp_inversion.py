# -*- coding: utf8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Spectral inversion routines for sourcespec.

:copyright:
    2012 Claudio Satriano <satriano@ipgp.fr>

    2013-2014 Claudio Satriano <satriano@ipgp.fr>,
              Emanuela Matrullo <matrullo@geologie.ens.fr>,
              Agnes Chounet <chounet@ipgp.fr>

    2015-2023 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
import logging
import numpy as np
from scipy.optimize import curve_fit, minimize, basinhopping
from scipy.signal import argrelmax
from obspy import Stream
from obspy.geodetics import gps2dist_azimuth
from sourcespec.ssp_spectral_model import (
    spectral_model, objective_func, callback)
from sourcespec.ssp_util import (
    mag_to_moment, source_radius, bsd, quality_factor, select_trace, smooth)
from sourcespec.ssp_data_types import (
    InitialValues, Bounds, SpectralParameter, StationParameters,
    SourceSpecOutput)
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
        spec_label = f'{spec.id} {spec.stats.instrtype}'
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
            raise RuntimeError(
                f'{statId}: unable to find a frequency range to compute Mw_0. '
                'This is possibly due to an uncommon spectrum '
                '(e.g., a resonance).'
            )
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
        idx1 = len(weight) // 2
    return idx0, idx1


def _spec_inversion(config, spec, spec_weight):
    """Invert one spectrum, return a Station() object."""
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
    statId = f'{spec.id} {spec.stats.instrtype}'
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
    logger.info(f'{statId}: initial values: {initial_values}')
    bounds = Bounds(config, spec, initial_values)
    bounds.Mw_min = np.nanmin(ydata[idx0: idx1]) * 0.9
    bounds.Mw_max = np.nanmax(ydata[idx0: idx1]) * 1.1
    if t_star_min is not None:
        bounds.t_star_min = t_star_min
    if t_star_max is not None:
        bounds.t_star_max = t_star_max
    logger.info(f'{statId}: bounds: {bounds}')
    try:
        params_opt, params_err, misfit = _curve_fit(
            config, spec, weight, yerr, initial_values, bounds)
    except (RuntimeError, ValueError) as m:
        raise RuntimeError(f'{m}\n{statId}: unable to fit spectral model')

    Mw, fc, t_star = params_opt
    Mw_err, fc_err, t_star_err = params_err
    inverted_par_str = f'Mw: {Mw:.4f}; fc: {fc:.4f}; t_star: {t_star:.4f}'
    logger.info(f'{statId}: optimal values: {inverted_par_str}')
    logger.info(f'{statId}: misfit: {misfit:.3f}')

    if np.isclose(fc, bounds.fc_min, rtol=0.1):
        raise ValueError(
            f'{statId}: optimal fc within 10% of fc_min: '
            f'{fc:.3f} ~= {bounds.fc_min:.3f}: ignoring inversion results'
        )

    if np.isclose(fc, bounds.fc_max, rtol=1e-4):
        raise ValueError(
            f'{statId}: optimal fc within 0.1% of fc_max: '
            f'{fc:.3f} ~= {bounds.fc_max:.3f}: ignoring inversion results'
        )

    misfit_max = config.pi_misfit_max or np.inf
    if misfit > misfit_max:
        raise ValueError(
            f'{statId}: misfit larger than pi_misfit_max: '
            f'{misfit:.3f} > {misfit_max:.3f}: ignoring inversion results'
        )

    # Check post-inversion bounds for t_star and fc
    pi_t_star_min, pi_t_star_max =\
        config.pi_t_star_min_max or (-np.inf, np.inf)
    if not (pi_t_star_min <= t_star <= pi_t_star_max):
        raise ValueError(
            f'{statId}: t_star: {t_star:.3f} not in allowed range '
            f'[{pi_t_star_min:.3f}, {pi_t_star_max:.3f}]: '
            'ignoring inversion results'
        )
    pi_fc_min, pi_fc_max = config.pi_fc_min_max or (-np.inf, np.inf)
    if not (pi_fc_min <= fc <= pi_fc_max):
        raise ValueError(
            f'{statId}: fc: {fc:.3f} not in allowed range '
            f'[{pi_fc_min:.3f}, {pi_fc_max:.3f}]: ignoring inversion results'
        )

    station_pars = StationParameters(
        id=spec.id, instrument_type=spec.stats.instrtype,
        latitude=stla, longitude=stlo,
        hypo_dist_in_km=spec.stats.hypo_dist,
        epi_dist_in_km=spec.stats.epi_dist,
        azimuth=az)
    station_pars.Mw = SpectralParameter(
        id='Mw', value=Mw,
        lower_uncertainty=Mw_err[0], upper_uncertainty=Mw_err[1],
        confidence_level=68.2, format='{:.2f}')
    station_pars.fc = SpectralParameter(
        id='fc', value=fc,
        lower_uncertainty=fc_err[0], upper_uncertainty=fc_err[1],
        confidence_level=68.2, format='{:.3f}')
    station_pars.t_star = SpectralParameter(
        id='t_star', value=t_star,
        lower_uncertainty=t_star_err[0], upper_uncertainty=t_star_err[1],
        confidence_level=68.2, format='{:.3f}')

    # additional parameters, computed from fc, Mw and t_star
    vs = config.hypo.vs
    travel_time = spec.stats.travel_times[config.wave_type[0]]
    # seismic moment
    station_pars.Mo = SpectralParameter(
        id='Mw', value=mag_to_moment(Mw), format='{:.3e}')
    # source radius in meters
    station_pars.radius = SpectralParameter(
        id='radius', value=source_radius(fc, vs*1e3), format='{:.3f}')
    # Brune stress drop in MPa
    station_pars.bsd = SpectralParameter(
        id='bsd', value=bsd(station_pars.Mo.value, station_pars.radius.value),
        format='{:.3e}')
    # quality factor
    station_pars.Qo = SpectralParameter(
        id='Qo', value=quality_factor(travel_time, t_star), format='{:.1f}')

    # Check post-inversion bounds for bsd
    pi_bsd_min, pi_bsd_max = config.pi_bsd_min_max or (-np.inf, np.inf)
    if not (pi_bsd_min <= station_pars.bsd.value <= pi_bsd_max):
        raise ValueError(
            f'{statId}: bsd: {station_pars.bsd.value:.3e} '
            f'not in allowed range [{pi_bsd_min:.3e}, {pi_bsd_max:.3e}]: '
            'ignoring inversion results'
        )

    # additional parameter errors, computed from fc, Mw and t_star
    # seismic moment
    Mw_min = Mw - Mw_err[0]
    Mw_max = Mw + Mw_err[1]
    Mo_min = mag_to_moment(Mw_min)
    Mo_max = mag_to_moment(Mw_max)
    station_pars.Mo.lower_uncertainty = station_pars.Mo.value - Mo_min
    station_pars.Mo.upper_uncertainty = Mo_max - station_pars.Mo.value
    station_pars.Mo.confidence_level = 68.2
    # source radius in meters
    fc_min = fc - fc_err[0]
    if fc_min <= 0:
        fc_min = freq_log[0]
    fc_max = fc + fc_err[1]
    radius_min = source_radius(fc_max, vs*1e3)
    radius_max = source_radius(fc_min, vs*1e3)
    station_pars.radius.lower_uncertainty =\
        station_pars.radius.value - radius_min
    station_pars.radius.upper_uncertainty =\
        radius_max - station_pars.radius.value
    station_pars.radius.confidence_level = 68.2
    # Brune stress drop in MPa
    bsd_min = bsd(Mo_min, radius_max)
    bsd_max = bsd(Mo_max, radius_min)
    station_pars.bsd.lower_uncertainty = station_pars.bsd.value - bsd_min
    station_pars.bsd.upper_uncertainty = bsd_max - station_pars.bsd.value
    station_pars.bsd.confidence_level = 68.2
    # quality factor
    t_star_min = t_star - t_star_err[0]
    if t_star_min <= 0:
        t_star_min = 0.001
    t_star_max = t_star + t_star_err[1]
    Qo_min = quality_factor(travel_time, t_star_max)
    Qo_max = quality_factor(travel_time, t_star_min)
    station_pars.Qo.lower_uncertainty = station_pars.Qo.value - Qo_min
    station_pars.Qo.upper_uncertainty = Qo_max - station_pars.Qo.value
    station_pars.Qo.confidence_level = 68.2

    return station_pars


def _synth_spec(config, spec, station_pars):
    """Return a stream with one or more synthetic spectra."""
    par = station_pars._params
    par_err = station_pars._params_err
    spec_st = Stream()
    params_opt = [par[key] for key in ('Mw', 'fc', 't_star')]

    chan_no_orientation = spec.stats.channel[:-1]

    freq = spec.get_freq()
    freq_log = spec.freq_log
    spec_synth = spec.copy()
    spec_synth.stats.channel = f'{chan_no_orientation}S'
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
        spec_synth.stats.channel = f'{chan_no_orientation}s'
        _params = list(params_opt)
        _params[-1] = 0
        spec_synth.data_mag = spectral_model(freq, *_params)
        spec_synth.data = mag_to_moment(spec_synth.data_mag)
        spec_synth.data_log_mag = spectral_model(freq_log, *_params)
        spec_synth.data_log = mag_to_moment(spec_synth.data_log_mag)
        spec_st.append(spec_synth)

    # Add an extra spectrum with no corner frequency
    if config.plot_spectra_no_fc:
        spec_synth = spec.copy()
        spec_synth.stats.channel = f'{chan_no_orientation}t'
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

    stations = {x.stats.station for x in spec_st}
    spectra = [sp for sta in stations for sp in spec_st.select(station=sta)]

    sspec_output = SourceSpecOutput()
    hypo = config.hypo
    sspec_output.event_info.event_id = hypo.evid
    sspec_output.event_info.longitude = hypo.longitude
    sspec_output.event_info.latitude = hypo.latitude
    sspec_output.event_info.depth_in_km = hypo.depth
    sspec_output.event_info.origin_time = hypo.origin_time
    for spec in sorted(spectra, key=lambda sp: sp.id):
        if spec.stats.channel[-1] != 'H':
            continue
        if spec.stats.ignore:
            continue
        spec_weight = select_trace(weight_st, spec.id, spec.stats.instrtype)
        try:
            station_pars = _spec_inversion(config, spec, spec_weight)
        except (RuntimeError, ValueError) as msg:
            logger.warning(msg)
            continue
        spec_st += _synth_spec(config, spec, station_pars)
        sspec_output.station_parameters[station_pars._id] = station_pars

    logger.info('Inverting spectra: done')
    return sspec_output
