# -*- coding: utf8 -*-
"""
Spectral inversion routines for sourcespec.

:copyright:
    2012 Claudio Satriano <satriano@ipgp.fr>

    2013-2014 Claudio Satriano <satriano@ipgp.fr>,
              Emanuela Matrullo <matrullo@geologie.ens.fr>,
              Agnes Chounet <chounet@ipgp.fr>

    2015-2016 Claudio Satriano <satriano@ipgp.fr>
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
from obspy.geodetics import gps2dist_azimuth
from sourcespec.ssp_setup import dprint
from sourcespec.ssp_spectral_model import (spectral_model, objective_func,
                                           callback)
from sourcespec.ssp_util import mag_to_moment, select_trace


class InitialValues():
    """Initial values for spectral inversion."""

    def __init__(self, Mw_0=None, fc_0=None, t_star_0=None):
        self.Mw_0 = Mw_0
        self.fc_0 = fc_0
        self.t_star_0 = t_star_0

    def __str__(self):
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
        self.bounds = ((self.Mw_min, self.Mw_max),
                       (self.fc_min, self.fc_max),
                       (self.t_star_min, self.t_star_max))
        self._fix_initial_values()

    def __str__(self):
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
            self.hd/(self.config.vs*np.array(self.config.Qo_min_max))
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

    def get_bounds(self):
        return self.bounds


def spectral_inversion(config, spec_st, weight_st, Ml):
    """Inversion of displacement spectra."""
    if config.noise_weighting:
        logging.info('Using noise weighting for inversion.')
    else:
        logging.info('Using frequency weighting for inversion.')
    if config.inv_algorithm == 'TNC':
        logging.info('Using truncated Newton algorithm for inversion.')
    elif config.inv_algorithm == 'LM':
        logging.info('Using Levenburg-Marquardt algorithm for inversion.')
    else:
        raise ValueError('Invalid choice of inversion algorithm.')

    sourcepar = dict()
    vs_m = config.vs*1000
    for station in set(x.stats.station for x in spec_st.traces):
        spec_st_sel = spec_st.select(station=station)
        for spec in spec_st_sel.traces:
            if spec.stats.channel[2] != 'H':
                continue
            dprint(station)

            # spectral amplitude is in Mw units
            amp = spec.data_mag

            # We calculate the initial value for Mw,
            # as an average of the first 5 values of amp
            Mw_0 = amp[0:5].mean()

            # We try to retrieve fc_0 from the configuration...
            fc_0 = config.fc_0
            # ...if it is not available, we calculate it
            if fc_0 is None:
                log_m0 = math.log10(mag_to_moment(Mw_0))
                log_beta = math.log10(vs_m)
                log_bsd = math.log10(config.bsd)
                log_fc = log_bsd - log_m0 + 3*log_beta - 0.935
                log_fc /= 3.
                fc_0 = math.pow(10, log_fc)
                logging.info('%s fc_0 autoset to: %.2f' % (spec.id, fc_0))

            # initial value for t_star
            t_star_0 = config.t_star_0

            hd = spec.stats.hypo_dist
            hd_m = hd*1000

            # azimuth computation
            coords = spec.stats.coords
            hypo = spec.stats.hypo
            stla = coords.latitude
            stlo = coords.longitude
            evla = hypo.latitude
            evlo = hypo.longitude
            geod = gps2dist_azimuth(evla, evlo, stla, stlo)
            az = geod[1]
            dprint('%s %s %f %f' % (station, spec.stats.instrtype, hd, az))
            loge = math.log10(math.e)
            coeff = math.pi*loge*hd_m/vs_m
            dprint('coeff= %f' % coeff)

            params_name = ('Mw', 'fc', 't_star')
            initial_values = InitialValues(Mw_0, fc_0, t_star_0)

            xdata = spec.get_freq()
            ydata = amp

            if config.noise_weighting:
                weight = select_trace(weight_st, spec.id, spec.stats.instrtype)
                # 'curve_fit' interprets 'yerr' as standard deviation vector
                # and calculates weights as 1/yerr^2 .
                # Therefore we build yerr as:
                yerr = 1./np.sqrt(weight)
            else:
                # Spectral weighting:
                #   config.weight for f<=f_weight
                #   1      for f> f_weight
                yerr = np.ones(len(ydata))
                yerr[xdata <= config.f_weight] = 1./math.sqrt(config.weight)
                weight = 1./np.power(yerr, 2)

            # Curve fitting using the Levenburg-Marquardt algorithm
            # or the truncated Newton algorithm (TNC), with bounds.
            try:
                if config.inv_algorithm == 'TNC':
                    bounds = Bounds(config, spec, initial_values)
                    logging.info('%s %s: initial values: %s' %
                                 (spec.id, spec.stats.instrtype,
                                  str(initial_values)))
                    logging.info('%s %s: bounds: %s' %
                                 (spec.id, spec.stats.instrtype, str(bounds)))
                    minimize_func = objective_func(xdata, ydata, weight)
                    res =\
                        minimize(minimize_func,
                                 x0=initial_values.get_params0(), method='TNC',
                                 callback=callback, bounds=bounds.get_bounds())
                    params_opt = res.x
                elif config.inv_algorithm == 'LM':
                    logging.info('%s %s: initial values: %s' %
                                 (spec.id, spec.stats.instrtype,
                                  str(initial_values)))
                    params_opt, params_cov =\
                        curve_fit(spectral_model,
                                  xdata, ydata,
                                  p0=initial_values.get_params0(),
                                  sigma=yerr)
            except RuntimeError:
                logging.warning('%s %s: unable to fit spectral model' %
                                (spec.id, spec.stats.instrtype))
                continue

            par = dict(zip(params_name, params_opt))
            par['Mo'] = mag_to_moment(par['Mw'])
            par['hyp_dist'] = hd
            par['az'] = az
            par['Ml'] = Ml  # FIXME: this is the network magnitude!
            statId = '%s %s' % (spec.id, spec.stats.instrtype)
            sourcepar[statId] = par

            spec_synth = spec.copy()
            spec_synth.stats.channel = spec.stats.channel[0:2] + 'S'
            spec_synth.stats.par = par
            spec_synth.data_mag = spectral_model(xdata, *params_opt)
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

    return sourcepar
