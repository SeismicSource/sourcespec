#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Main function for source_spec
# (c) 2012 Claudio Satriano <satriano@ipgp.fr>
# (c) 2013 Claudio Satriano <satriano@ipgp.fr>,
#          Emanuela Matrullo <matrullo@geologie.ens.fr>
from __future__ import division
import os
import logging
import math
import numpy as np
from scipy.optimize import curve_fit
from lib.ssp_setup import dprint, configure, setup_logging, ssp_exit
from lib.ssp_read_traces import read_traces
from lib.ssp_process_traces import process_traces
from lib.ssp_build_spectra import build_spectra
from lib.ssp_local_magnitude import local_magnitude
from lib.ssp_plot_spectra import plot_spectra
from lib.ssp_plot_spectranoise import plot_spectranoise
from lib.ssp_plot_weight import plot_weight
from lib.ssp_output import write_output
from lib.ssp_spectral_model import spectral_model
from lib.ssp_util import mag_to_moment
from obspy.core.util.geodetics import gps2DistAzimuth

def main():
    # Setup stage
    config = configure()
    setup_logging(config)
    st = read_traces(config)

    # Now that we (hopefully) have the evid
    # we rename the logfile to use the evid
    #TODO: improve this:
    evid = st.traces[0].stats.hypo.evid
    setup_logging(config, evid)

    # Deconvolve, filter, cut traces:
    proc_st, noise_st = process_traces(config, st)
    # Build spectra (amplitude in magnitude units)
    spec_st, specnoise_st, weight_st = build_spectra(config, proc_st, noise_st)

    #Ml = local_magnitude(config, proc_st)
    Ml = local_magnitude(config, st, deconvolve=True)

    # Inversion of displacement spectra
    # Spectral weighting:
    #   weight for f<=f_weight
    #   1      for f> f_weight
    #f_weight = config.f_weight
    #weight   = config.weight
    sourcepar = dict()
    vs_m = config.vs*1000
    for station in set(x.stats.station for x in spec_st.traces):
        spec_st_sel = spec_st.select(station=station)
        for spec in spec_st_sel.traces:
            if spec.stats.channel != 'H': continue
            dprint(station)

            spec_id = spec.id 
            #print spec_id, weight_st.select(id=spec_id)
            weight = weight_st.select(id=spec_id)[0] 

            # spectral amplitude is in Mw units
            amp = spec.data_mag

            # We calculate the initial value for Mw,
            # as an average of the first 5 values of amp
            Mw_0 = amp[0:5].mean()

            # We try to retrieve fc_0 from the configuration...
            try:
                fc_0 = config.fc_0
            # ...if it is not available, we calculate it
            except: 
                l_m0 = Mw_0 * 1.5 + 9.1
                l_beta = math.log10(vs_m)
                l_bsd = math.log10(config.bsd)
                l_fc = l_bsd - l_m0 + 3*l_beta - 0.935
                l_fc /= 3.
                fc_0 = math.pow(10, l_fc)

            # initial value for t_star
            t_star_0 = config.t_star_0
            
            # print initial values of fc, M0 and t_star
            dprint('INITIAL CORNER FREQUENCY= %f' % fc_0)
            dprint('INITIAL MOMENT MAGNITUDE= %f' % Mw_0)
            dprint('INITIAL T_STAR= %f' % t_star_0)

            hd   = spec.stats.hypo_dist
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
            #yerr = np.ones(len(ydata))
            # 'curve_fit' interprets 'yerr' as standard deviation vector and calculates
            # weights as 1/yerr^2 . Therefore we build yerr as:
            #yerr[xdata<=f_weight] = 1./math.sqrt(weight)
            yerr = 1./np.sqrt(weight)
            # Curve fitting using the Levenburg-Marquardt algorithm
            try:
                    params_opt, params_cov = curve_fit(spectral_model, xdata, ydata, p0=params_0, sigma=yerr)
            except RuntimeError:
                    logging.warning('Unable to fit spectral model for station: %s' % station)
            #print params_opt, params_cov
            par = dict(zip(params_name, params_opt))
            par['hyp_dist'] = hd
            par['az'] = az
            par['Ml'] = Ml #FIXME: this is the network magnitude!
            chanId = '%s.%s' % (station, spec.stats.instrtype)
            sourcepar[chanId] = par

            spec_synth = spec.copy()
            spec_synth.stats.channel = 'Synth'
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

    # Save output
    write_output(config, evid, sourcepar)

    # Plotting
    plot_spectra(config, spec_st)
    plot_spectranoise(config, specnoise_st)
    plot_weight(config, weight_st)

    ssp_exit()


if __name__ == '__main__':
    try:
        thismodule=os.path.basename(__file__).replace('.py','')
        this=__import__(thismodule)
        main()
    except KeyboardInterrupt:
        ssp_exit(1)
