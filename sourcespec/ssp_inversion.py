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


class InitialValues():
    """Initial values for spectral inversion."""

    def __init__(self, Mw_0=None, fc_0=None, t_star_0=None):
        self.Mw_0 = Mw_0
        self.fc_0 = fc_0
        self.t_star_0 = t_star_0

    def __str__(self):
        """String representation."""
        s = 'Mw_0: %s; ' % round(self.Mw_0, 4)
        s += 'fc_0: %s; ' % round(self.fc_0, 4)
        s += 't_star_0: %s' % round(self.t_star_0, 4)
        return s

    def get_params0(self):
        return (self.Mw_0, self.fc_0, self.t_star_0)


class Bounds():
    """Bounds for bounded spectral inversion."""

    def __init__(self, config, spec, initial_values):
        self.config = config
        self.spec = spec
        self.hd = spec.stats.hypo_dist
        self.ini_values = initial_values
        self.Mw_min, self.Mw_max = self._check_minmax(config.Mw_min_max)
        self.fc_min, self.fc_max = self._check_minmax(config.fc_min_max)
        if config.Qo_min_max is None:
            self.t_star_min, self.t_star_max =\
                self._check_minmax(config.t_star_min_max)
        else:
            self.t_star_min, self.t_star_max = self._Qo_to_t_star()
        self.set_bounds()
        self._fix_initial_values()

    def __str__(self):
        """String representation."""
        s = 'Mw: %s, %s; ' %\
            tuple(round(x, 4) if x is not None else x for x in self.bounds[0])
        s += 'fc: %s, %s; ' %\
            tuple(round(x, 4) if x is not None else x for x in self.bounds[1])
        s += 't_star: %s, %s' %\
            tuple(round(x, 4) if x is not None else x for x in self.bounds[2])
        return s

    def _check_minmax(self, minmax):
        if minmax is None:
            return (None, None)
        else:
            return minmax

    def _Qo_to_t_star(self):
        t_star_max, t_star_min =\
            self.hd/(self.config.hypo.vs*np.array(self.config.Qo_min_max))
        return self._nan_to_none((t_star_min, t_star_max))

    def _fix_initial_values(self):
        if (self.Mw_min is not None and
            self.Mw_max is not None and
            (self.ini_values.Mw_0 is None or
             self.ini_values.Mw_0 <= self.Mw_min or
             self.ini_values.Mw_0 >= self.Mw_max)):
            Mw_0 = (self.Mw_max + self.Mw_min) / 2.
            logging.warning('%s %s: initial Mw value: %s outside '
                            'bounds. Using bound average: %s' %
                            (self.spec.id, self.spec.stats.instrtype,
                             self.ini_values.Mw_0, round(Mw_0, 4)))
            self.ini_values.Mw_0 = Mw_0
        if (self.fc_min is not None and
            self.fc_max is not None and
            (self.ini_values.fc_0 is None or
             self.ini_values.fc_0 <= self.fc_min or
             self.ini_values.fc_0 >= self.fc_max)):
            fc_0 = (self.fc_max + self.fc_min) / 2.
            logging.warning('%s %s: initial fc value: %s outside '
                            'bounds. Using bound average: %s' %
                            (self.spec.id, self.spec.stats.instrtype,
                             self.ini_values.fc_0, round(fc_0, 4)))
            self.ini_values.fc_0 = fc_0
        if (self.t_star_min is not None and
            self.t_star_max is not None and
            (self.ini_values.t_star_0 is None or
             self.ini_values.t_star_0 <= self.t_star_min or
             self.ini_values.t_star_0 >= self.t_star_max)):
            t_star_0 = (self.t_star_max + self.t_star_min) / 2.
            logging.warning('%s %s: initial t_star value: %s outside '
                            'bounds. Using bound average: %s' %
                            (self.spec.id, self.spec.stats.instrtype,
                             self.ini_values.t_star_0, round(t_star_0, 4)))
            self.ini_values.t_star_0 = t_star_0

    def __call__(self, **kwargs):
        """Interface for basin-hopping."""
        params = kwargs['x_new']
        tmin = bool(np.all(params >= self.params_min))
        tmax = bool(np.all(params <= self.params_max))
        return tmin and tmax

    def set_bounds(self):
        self.bounds = ((self.Mw_min, self.Mw_max),
                       (self.fc_min, self.fc_max),
                       (self.t_star_min, self.t_star_max))
        self.params_min = (self.Mw_min, self.fc_min, self.t_star_min)
        self.params_max = (self.Mw_max or 1e300, self.fc_max or 1e300,
                           self.t_star_max or 1e300)

    def get_bounds(self):
        """Get bounds for minimize()."""
        return self.bounds

    def get_bounds_curve_fit(self):
        """Get bounds for curve-fit()."""
        bnds = np.array(self.bounds, dtype=float).T
        if np.all(np.isnan(bnds)):
            return None
        bnds[0, np.isnan(bnds[0])] = -1e100
        bnds[1, np.isnan(bnds[1])] = 1e100
        return bnds


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
    for station in set(x.stats.station for x in spec_st.traces):
        spec_st_sel = spec_st.select(station=station)
        for spec in spec_st_sel.traces:
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
            bounds.set_bounds()
            logging.info('%s %s: bounds: %s' %
                         (spec.id, spec.stats.instrtype, str(bounds)))

            # Curve fitting using "curve_fit()" (Levenberg-Marquardt algorithm
            # if no bounds or Trust Region Reflective algorithm if bounds)
            # or the truncated Newton algorithm (TNC) with bounds.
            try:
                if config.inv_algorithm == 'TNC':
                    minimize_func = objective_func(freq_log, ydata, weight)
                    res =\
                        minimize(minimize_func,
                                 x0=initial_values.get_params0(), method='TNC',
                                 callback=callback, bounds=bounds.get_bounds())
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
