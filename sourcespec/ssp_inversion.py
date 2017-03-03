# -*- coding: utf8 -*-
"""
Spectral inversion routines for sourcespec.

:copyright:
    2012 Claudio Satriano <satriano@ipgp.fr>

    2013-2014 Claudio Satriano <satriano@ipgp.fr>,
              Emanuela Matrullo <matrullo@geologie.ens.fr>,
              Agnes Chounet <chounet@ipgp.fr>

    2015-2017 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement, Version 2.1
    (http://www.cecill.info/index.en.html)
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import logging
import math
import numpy as np
from scipy.optimize import curve_fit, minimize
from scipy.signal import argrelmax
from obspy.geodetics import gps2dist_azimuth
from sourcespec.ssp_spectral_model import (spectral_model, objective_func,
                                           callback)
from sourcespec.ssp_util import mag_to_moment, select_trace
from sourcespec.ssp_radiated_energy import radiated_energy
from sourcespec.ssp_inversion_types import InitialValues, Bounds


def _curve_fit(config, spec, weight, yerr, initial_values, bounds):
    """
    Curve fitting.

    Uses "curve_fit()" (Levenberg-Marquardt algorithm if no bounds
    or Trust Region Reflective algorithm if bounds)
    or the truncated Newton algorithm (TNC) with bounds.
    """
    freq_log = spec.freq_log
    ydata = spec.data_log_mag
    try:
        if config.inv_algorithm == 'TNC':
            minimize_func = objective_func(freq_log, ydata, weight)
            res =\
                minimize(minimize_func,
                         x0=initial_values.get_params0(), method='TNC',
                         callback=callback, bounds=bounds.bounds)
            params_opt = res.x
            # trick: use curve_fit() bounded to params_opt
            # to get the covariance
            _, params_cov = curve_fit(spectral_model, freq_log, ydata,
                                      p0=params_opt, sigma=yerr,
                                      bounds=(params_opt-(1e-10),
                                              params_opt+(1e-10)))
        elif config.inv_algorithm == 'LM':
            logging.info('%s %s: initial values: %s' %
                         (spec.id, spec.stats.instrtype,
                          str(initial_values)))
            bnds = bounds.get_bounds_curve_fit()
            if bnds is not None:
                logging.info('Trying to use using Levenberg-Marquardt '
                             'algorithm with bounds. Switching to the '
                             'Trust Region Reflective algorithm.')
            params_opt, params_cov =\
                curve_fit(spectral_model,
                          freq_log, ydata,
                          p0=initial_values.get_params0(),
                          sigma=yerr,
                          bounds=bounds.get_bounds_curve_fit())
        elif config.inv_algorithm == 'BH':
            from scipy.optimize import basinhopping
            minimize_func = objective_func(freq_log, ydata, weight)
            res = basinhopping(minimize_func,
                               x0=initial_values.get_params0(),
                               niter=100, accept_test=bounds)
            params_opt = res.x
    except RuntimeError:
        logging.warning('%s %s: unable to fit spectral model' %
                        (spec.id, spec.stats.instrtype))
        return None, None
    return params_opt, params_cov


def spectral_inversion(config, spec_st, weight_st, Ml):
    """Inversion of displacement spectra."""
    if config.weighting == 'noise':
        logging.info('Using noise weighting for inversion.')
    elif config.weighting == 'frequency':
        logging.info('Using frequency weighting for inversion.')
    elif config.weighting is None:
        logging.info('Using no weighting for inversion.')
    if config.inv_algorithm == 'TNC':
        logging.info('Using truncated Newton algorithm for inversion.')
    elif config.inv_algorithm == 'LM':
        logging.info('Using Levenberg-Marquardt algorithm for inversion.')
    elif config.inv_algorithm == 'BH':
        logging.info('Using basin-hopping algorithm for inversion.')
    else:
        raise ValueError('Invalid choice of inversion algorithm.')

    sourcepar = dict()
    sourcepar_err = dict()
    stations = set(x.stats.station for x in spec_st)
    spectra = [sp for sta in stations for sp in spec_st.select(station=sta)]
    for spec in spectra:
        if spec.stats.channel[2] != 'H':
            continue

        # azimuth computation
        coords = spec.stats.coords
        hypo = spec.stats.hypo
        stla = coords.latitude
        stlo = coords.longitude
        evla = hypo.latitude
        evlo = hypo.longitude
        geod = gps2dist_azimuth(evla, evlo, stla, stlo)
        az = geod[1]

        freq = spec.get_freq()
        freq_log = spec.freq_log
        ydata = spec.data_log_mag

        noise_weight = select_trace(
            weight_st, spec.id, spec.stats.instrtype)
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

        # Find the frequency range to compute Mw_0:
        # we start where signal-to-noise becomes strong
        idx0 = np.where(noise_weight > 0.5)[0][0]
        # we stop at the first max of signal-to-noise (proxy for fc)
        idx_max = argrelmax(noise_weight)[0]
        # just keep the indexes for maxima > 0.5
        idx_max = [idx for idx in idx_max if noise_weight[idx] > 0.5]
        idx1 = idx_max[0]

        # We calculate the initial value for Mw as an average
        Mw_0 = np.nanmean(ydata[idx0: idx1])
        # We try to retrieve fc_0 from the configuration...
        fc_0 = config.fc_0
        # ...if it is not available, we calculate it
        if fc_0 is None:
            # fc_0 = freq_log[idx1]
            log_m0 = math.log10(mag_to_moment(Mw_0))
            log_beta = math.log10(config.hypo.vs*1000.)
            log_bsd = math.log10(config.bsd)
            log_fc = log_bsd - log_m0 + 3*log_beta - 0.935
            log_fc /= 3.
            fc_0 = math.pow(10, log_fc)
            logging.info('%s fc_0 autoset to: %.2f' % (spec.id, fc_0))
        # initial value for t_star
        t_star_0 = config.t_star_0
        initial_values = InitialValues(Mw_0, fc_0, t_star_0)
        logging.info('%s %s: initial values: %s' %
                     (spec.id, spec.stats.instrtype,
                      str(initial_values)))

        bounds = Bounds(config, spec, initial_values)
        bounds.Mw_min = Mw_0 - 0.1
        bounds.Mw_max = Mw_0 + 0.1
        logging.info('%s %s: bounds: %s' %
                     (spec.id, spec.stats.instrtype, str(bounds)))

        params_opt, params_cov = _curve_fit(
            config, spec, weight, yerr, initial_values, bounds)
        if params_opt is None:
            continue

        params_name = ('Mw', 'fc', 't_star')
        par = dict(zip(params_name, params_opt))
        par['Mo'] = mag_to_moment(par['Mw'])
        par['hyp_dist'] = spec.stats.hypo_dist
        par['az'] = az
        par['Ml'] = Ml  # FIXME: this is the network magnitude!

        error = np.sqrt(params_cov.diagonal())
        par_err = dict(zip(params_name, error))

        statId = '%s %s' % (spec.id, spec.stats.instrtype)
        sourcepar[statId] = par
        sourcepar_err[statId] = par_err

        spec_synth = spec.copy()
        spec_synth.stats.channel = spec.stats.channel[0:2] + 'S'
        spec_synth.stats.par = par
        spec_synth.stats.par_err = par_err
        spec_synth.data_mag = spectral_model(freq, *params_opt)
        spec_synth.data = mag_to_moment(spec_synth.data_mag)
        spec_st.append(spec_synth)

        # Add an extra spectrum with no attenuation
        if config.plot_spectra_no_attenuation:
            spec_synth = spec.copy()
            spec_synth.stats.channel = spec.stats.channel[0:2] + 's'
            params_opt[-1] = 0
            spec_synth.data_mag = spectral_model(freq, *params_opt)
            spec_synth.data = mag_to_moment(spec_synth.data_mag)
            spec_st.append(spec_synth)

    # Filter stations with negative t_star
    # or with anomalous corner frequencies
    f1 = config.min_corner_freq
    f2 = config.max_corner_freq
    # Make a copy of sourcepar.keys() since the dictionary
    # may change during iteration
    for statId in list(sourcepar.keys()):
        par = sourcepar[statId]
        t_star = par['t_star']
        if t_star < 0:
            logging.warning('Ignoring station: %s t_star: %f' %
                            (statId, t_star))
            sourcepar.pop(statId, None)
        fc = par['fc']
        if fc < f1 or fc > f2:
            logging.warning('Ignoring station: %s fc: %f' % (statId, fc))
            sourcepar.pop(statId, None)

    radiated_energy(config, spec_st, sourcepar)

    return sourcepar, sourcepar_err
