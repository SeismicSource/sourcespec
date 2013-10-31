# -*- coding: utf8 -*-
# ssp_inversion.py
#
# Spectral inversion routine for source_spec
# (c) 2012 Claudio Satriano <satriano@ipgp.fr>
# (c) 2013 Claudio Satriano <satriano@ipgp.fr>,
#          Emanuela Matrullo <matrullo@geologie.ens.fr>,
#          Agnes Chounet <chounet@ipgp.fr>
from __future__ import division
import logging
import math
import numpy as np
from scipy.optimize import curve_fit, minimize
from ssp_setup import dprint
from ssp_spectral_model import spectral_model, objective_func, callback
from ssp_util import mag_to_moment, select_trace
from obspy.core.util.geodetics import gps2DistAzimuth


class Bounds():

    def __init__(self, hd, config):
        self.hd = float(hd)
        self.config = config
        Mw_min_max = self.__nan_to_none__(config.Mw_min_max)
        fc_min_max = self.__nan_to_none__(config.fc_min_max)
        if any(math.isnan(x) for x in config.Qo_min_max):
            t_star_min_max = self.__nan_to_none__(config.t_star_min_max)
        else:
            t_star_min_max = self.__Qo_to_t_star__()
        self.bounds = [Mw_min_max,
                       fc_min_max,
                       t_star_min_max]

    def __nan_to_none__(self, val):
        ret = [None if math.isnan(x) else x for x in tuple(val)]
        return tuple(ret)

    def __Qo_to_t_star__(self):
        t_star_max, t_star_min =\
            self.hd/(self.config.vs*np.array(self.config.Qo_min_max))
        return self.__nan_to_none__((t_star_min, t_star_max))

    def __str__(self):
        s = 'Mw: %s, %s; ' % tuple(map(str, self.bounds[0]))
        s += 'fc: %s, %s; ' % tuple(map(str, self.bounds[1]))
        s += 't_star: %s, %s' % tuple(map(str, self.bounds[2]))
        return s

    def get_bounds(self):
        return self.bounds



def spectral_inversion(config, spec_st, weight_st, Ml):
    '''
    Inversion of displacement spectra
    '''

    if config.noise_weighting:
        logging.info('Using noise weighting for inversion.')
    else:
        logging.info('Using frequency weighting for inversion.')
    if config.inv_algorithm == 'TNC':
        logging.info('Using truncated Newton algorithm for inversion.')
    elif config.inv_algorithm == 'LM':
        logging.info('Using Levenburg-Marquardt algorithm for inversion.')
    else:
        raise ValueError, 'Invalid choice of inversion algorithm.'

    sourcepar = dict()
    vs_m = config.vs*1000
    for station in set(x.stats.station for x in spec_st.traces):
        spec_st_sel = spec_st.select(station=station)
        for spec in spec_st_sel.traces:
            if spec.stats.channel != 'H': continue
            dprint(station)

            # spectral amplitude is in Mw units
            amp = spec.data_mag

            # We calculate the initial value for Mw,
            # as an average of the first 5 values of amp
            Mw_0 = amp[0:5].mean()

            # We try to retrieve fc_0 from the configuration...
            fc_0 = config.fc_0
            # ...if it is not available, we calculate it
            if math.isnan(fc_0):
                l_m0 = Mw_0 * 1.5 + 9.1
                l_beta = math.log10(vs_m)
                l_bsd = math.log10(config.bsd)
                l_fc = l_bsd - l_m0 + 3*l_beta - 0.935
                l_fc /= 3.
                fc_0 = math.pow(10, l_fc)
                logging.info('%s fc_0 autoset to: %.2f' % (spec.id, fc_0))

            # initial value for t_star
            t_star_0 = config.t_star_0

            # print initial values of fc, M0 and t_star
            dprint('INITIAL CORNER FREQUENCY= %f' % fc_0)
            dprint('INITIAL MOMENT MAGNITUDE= %f' % Mw_0)
            dprint('INITIAL T_STAR= %f' % t_star_0)

            hd = spec.stats.hypo_dist
            hd_m = hd*1000

            # azimuth computation
            coords = spec.stats.coords
            hypo = spec.stats.hypo
            stla = coords.latitude
            stlo = coords.longitude
            evla = hypo.latitude
            evlo = hypo.longitude
            geod = gps2DistAzimuth(evla, evlo, stla, stlo)
            az   = geod[1]
            dprint('%s %s %f %f' % (station, spec.stats.instrtype, hd, az))
            loge = math.log10(math.e)
            coeff = math.pi*loge*hd_m/vs_m
            dprint('coeff= %f' % coeff)

            params_name = ('Mw', 'fc', 't_star')
            params_0 = np.array([Mw_0, fc_0, t_star_0])

            xdata = spec.get_freq()
            ydata = amp

            if config.noise_weighting:
                weight = select_trace(weight_st, spec.id, spec.stats.instrtype)
                # 'curve_fit' interprets 'yerr' as standard deviation vector and calculates
                # weights as 1/yerr^2 . Therefore we build yerr as:
                yerr = 1./np.sqrt(weight)
            else:
                # Spectral weighting:
                #   config.weight for f<=f_weight
                #   1      for f> f_weight
                yerr = np.ones(len(ydata))
                yerr[xdata<=config.f_weight] = 1./math.sqrt(config.weight)
                weight = 1./np.power(yerr, 2)

            # Curve fitting using the Levenburg-Marquardt algorithm
            # or the truncated Newton algorithm (TNC), with bounds.
            try:
                if config.inv_algorithm == 'TNC':
                    bounds = Bounds(hd, config)
                    logging.info('%s %s: bounds: %s' %
                            (spec.id, spec.stats.instrtype, str(bounds)))
                    minimize_func = objective_func(xdata, ydata, weight)
                    res = minimize(
                        minimize_func, params_0, method='TNC',
                        callback=callback, bounds=bounds.get_bounds())
                    params_opt = res.x
                elif config.inv_algorithm == 'LM':
                    params_opt, params_cov = curve_fit(
                        spectral_model, xdata, ydata, p0=params_0, sigma=yerr)
            except RuntimeError:
                    logging.warning('Unable to fit spectral model for station: %s' % station)

            par = dict(zip(params_name, params_opt))
            par['hyp_dist'] = hd
            par['az'] = az
            par['Ml'] = Ml #FIXME: this is the network magnitude!
            chanId = '%s.%s' % (station, spec.stats.instrtype)
            sourcepar[chanId] = par

            spec_synth = spec.copy()
            spec_synth.stats.channel = 'Synth'
            spec_synth.stats.par = par
            spec_synth.data_mag = spectral_model(xdata, *params_opt)
            spec_synth.data = mag_to_moment(spec_synth.data_mag)
            spec_st.append(spec_synth)

    # Filter stations with negative t_star
    # or with anomalous corner frequencies
    f1 = config.min_corner_freq
    f2 = config.max_corner_freq
    for statId in sourcepar.keys():
        par = sourcepar[statId]
        t_star = par['t_star']
        if t_star < 0:
            logging.warning('Ignoring station: %s t_star: %f' % (statId, t_star))
            sourcepar.pop(statId, None)
        fc = par['fc']
        if fc < f1 or fc > f2:
            logging.warning('Ignoring station: %s fc: %f' % (statId, fc))
            sourcepar.pop(statId, None)

    return sourcepar
