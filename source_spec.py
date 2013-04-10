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
from obspy.core import Stream
from lib.ssp_setup import dprint, configure, setup_logging, ssp_exit
from lib.ssp_read_traces import read_traces
from lib.ssp_util import *
from lib.ssp_plot_spectra import *
from lib.ssp_output import *
from lib.ssp_spectral_model import spectral_model
from lib import spectrum
from copy import deepcopy, copy
from obspy.core.util.geodetics import gps2DistAzimuth
from obspy.signal import estimateMagnitude
from obspy.signal.konnoohmachismoothing import konnoOhmachiSmoothing

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

    magnitudes = []

    # Loop on stations for building spectra and the spec_st object 
    spec_st = Stream()
    for trace in st.traces:
        traceId = trace.getId()

        # copy trace for local magnitude estimation
        trace_cp = copy(trace)
        trace_cp.stats = deepcopy(trace.stats)


        # check if the trace has (significant) signal
        # since the count value is generally huge, we need to demean twice
        # to take into account for the rounding error
        trace.detrend(type='constant')
        trace.detrend(type='constant')
        rms2 = np.power(trace.data, 2).sum()
        rms_min = 1e-10 #TODO: parametrize?
        if rms2 <= rms_min:
            logging.warning('%s: Trace RMS smaller than %g: skipping trace' % (traceId, rms_min))
            continue

        # Remove instrument response
        if remove_instr_response(trace, config.correct_sensitivity_only,
                                 config.pre_filt) == None:
            logging.warning('%s: Unable to remove instrument response: skipping trace' % traceId)
            continue

        stats = trace.stats
        comp  = stats.channel
        # skip vertical components
        if comp[-1] == 'Z':
            continue
        station = stats.station
        dprint('%s %s' % (station, comp))

        # compute hypocentral distance hd
        hd = hypo_dist(trace)
        if hd == None:
            logging.warning('%s: Unable to compute hypocentral distance: skipping trace' % traceId)
            continue
        hd_m = hd*1000

        # S time window
        s_arrival_time = swave_arrival(trace, config.vs)
        t1 = s_arrival_time - config.pre_s_time
        t2 = t1 + config.s_win_length
        #trace_cut = trace.slice(t1, t2)
        trace_cut = copy(trace)
        trace_cut.stats = deepcopy(trace.stats)
        trace_cut.trim(starttime=t1, endtime=t2, pad=True, fill_value=0)

        npts = len(trace_cut.data)
        if npts == 0:
            logging.warning('%s: No data for the selected cut interval: skipping trace' % traceId)
            continue
        nzeros = len(np.where(trace_cut.data==0)[0])
        if nzeros > npts/4:
            logging.warning('%s: Too many gaps for the selected cut interval: skipping trace' % traceId)
            continue
        
        # TODO: parameterize
        # coefficient for converting displ spectrum
        # to seismic moment (Aki&Richards,1980)
        vs_m = config.vs*1000
        vs3  = pow(vs_m,3) 
        coeff = 4 * math.pi * vs3 * config.rho /config.rps/2.
        dprint('coeff= %f' % coeff)

        # compute S-wave (displacement) spectra from
        # accelerometers and velocimeters, uncorrected for attenuation,
        # corrected for instrumental constants, normalized by
        # hypocentral distance
        instrtype = stats.instrtype
        if instrtype == 'acc':
            nint = 2 #number of intergrations to perform
            # band-pass frequencies:
            # TODO: calculate from sampling rate?
            bp_freqmin = config.bp_freqmin_acc
            bp_freqmax = config.bp_freqmax_acc
            # cut frequencies:
            freq1 = config.freq1_acc
            freq2 = config.freq2_acc
        elif instrtype == 'shortp':
            nint = 1 #number of intergrations to perform
            # band-pass frequencies:
            # TODO: calculate from sampling rate?
            bp_freqmin = config.bp_freqmin_shortp
            bp_freqmax = config.bp_freqmax_shortp
            # cut frequencies:
            freq1 = config.freq1_shortp
            freq2 = config.freq2_shortp
        elif instrtype == 'broadb':
            nint = 1
            bp_freqmin = config.bp_freqmin_broadb
            bp_freqmax = config.bp_freqmax_broadb
            # cut frequencies:
            freq1 = config.freq1_broadb
            freq2 = config.freq2_broadb
        else:
            dprint('%s: Unknown instrument type: %s: skipping trace' % (traceId, instrtype))
            continue

        # remove the mean...
        trace_cut.detrend(type='constant')
        # ...and the linear trend...
        trace_cut.detrend(type='linear')
        trace_cut.filter(type='bandpass', freqmin=bp_freqmin, freqmax=bp_freqmax)
        # ...and taper
        cosine_taper(trace_cut.data, width=0.5)

        # normalization for the hypocentral distance
        trace_cut.data *= hd_m

        # calculate fft
        spec = spectrum.do_spectrum(trace_cut)
        spec.stats.instrtype = instrtype
        spec.stats.coords    = stats.coords
        spec.stats.hypo      = stats.hypo
        spec.stats.hypo_dist = stats.hypo_dist

        # Integrate in frequency domain (divide by the pulsation omega)
        for i in range(0,nint):
            spec.data /= (2 * math.pi * spec.get_freq())

        # smooth the abs of fft
        #data_smooth = smooth(spec.data, 6)
        #datanoise_smooth = smooth(specnoise.data, 6)
        data_smooth = konnoOhmachiSmoothing(spec.data, spec.get_freq(),40,normalize=True)

        # Uncomment these lines to see the effect of smoothing
        #import matplotlib.pyplot as plt
        #plt.figure()
        #plt.loglog(spec.get_freq(), spec.data, color='gray')
        #plt.loglog(spec.get_freq(), data_smooth)
        #plt.grid(True)
        #plt.show()
        #ssp_exit()

        spec.data = data_smooth
        # convert to seismic moment
        spec.data *= coeff

        # Cut the spectrum between freq1 and freq2
        spec_cut = spec.slice(freq1, freq2)
        spec_st.append(spec_cut)

        # use min/max amplitude for local magnitude estimation

        # remove the mean...
        trace_cp.detrend(type='constant')
        # ...and the linear trend...
        trace_cp.detrend(type='linear')
        # ...filter
        trace_cp.filter(type='bandpass',  freqmin=0.1, freqmax=20)
        #trace_cp.filter(type='bandpass',  freqmin=0.5, freqmax=20)
        # ...and taper
        cosine_taper(trace_cp.data, width=0.5)

        delta_amp = trace_cp.data.max() - trace_cp.data.min() 
        #delta_amp = max(abs(trace_cp.data))
        delta_t = trace_cp.data.argmax() - trace_cp.data.argmin()    
        delta_t = delta_t / trace_cp.stats.sampling_rate

        #estimate Magnitude 
        ml = estimateMagnitude(trace_cp.stats.paz, delta_amp, delta_t, hd)

        magnitudes.append(ml)
        #mlId = '%s %s.%s %.1f' % (evid, station, spec.stats.instrtype, ml)
        print '%s %s.%s %s %.1f' % (evid, station, spec.stats.instrtype, "Ml", ml)
    # end of loop on stations for building spectra 

    # average local magnitude
    Ml = np.mean(magnitudes)
    print '\nNetwork Local Magnitude: %.2f' % Ml

    # convert the spectral amplitudes to moment magnitude
    for spec in spec_st.traces:
        spec.data = (np.log10(spec.data) - 9.1 ) / 1.5

    # Add to spec_st the "horizontal" component, obtained from the
    # modulus of the N-S and E-W components.
    for station in set(x.stats.station for x in spec_st.traces):
        spec_st_sel = spec_st.select(station=station)
        for instrtype in set(x.stats.instrtype for x in spec_st_sel):
            spec_h = None
            for spec in spec_st_sel.traces:
                # this should never happen:
                if spec.stats.instrtype != instrtype:
                    continue
                if spec_h == None:
                    spec_h = spec.copy()
                    spec_h.stats.channel = 'H'
                else:
                    data_h = spec_h.data
                    data   = spec.data
                    data_h = np.power(data_h, 2) + np.power(data, 2)
                    data_h = np.sqrt(data_h / 2) #divide by the number of horizontal components
                    spec_h.data = data_h
            spec_st.append(spec_h)
            

    # Inversion of displacement spectra
    loge = math.log10(math.e)
    # Spectral weighting:
    #   weight for f<=f_weight
    #   1      for f> f_weight
    f_weight = config.f_weight
    weight   = config.weight
    sourcepar = dict()
    for station in set(x.stats.station for x in spec_st.traces):
        spec_st_sel = spec_st.select(station=station)
        for spec in spec_st_sel.traces:
            if spec.stats.channel != 'H': continue
            dprint(station)

            # spectral amplitude is in Mw units
            amp = spec.data

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
            coeff = math.pi*loge*hd_m/vs_m
            dprint('coeff= %f' % coeff)

            params_name = ('Mw', 'fc', 't_star')
            params_0 = np.array([Mw_0, fc_0, t_star_0])

            xdata = spec.get_freq()
            ydata = amp
            yerr = np.ones(len(ydata))
            # 'curve_fit' interprets 'yerr' as standard deviation vector and calculates
            # weights as 1/yerr^2 . Therefore we build yerr as:
            yerr[xdata<=f_weight] = 1./math.sqrt(weight)
            # Curve fitting using the Levenburg-Marquardt algorithm
            params_opt, params_cov = curve_fit(spectral_model, xdata, ydata, p0=params_0, sigma=yerr)
            try:
                    params_opt, params_cov = curve_fit(spectral_model, xdata, ydata, p0=params_0, sigma=yerr)
            except RuntimeError:
                    logging.warning('Unable to fit spectral model for station: %s' % station)
            #print params_opt, params_cov
            par = dict(zip(params_name, params_opt))
            par['hyp_dist'] = hd
            par['az'] = az
            par['Ml'] = Ml
            chanId = '%s.%s' % (station, spec.stats.instrtype)
            sourcepar[chanId] = par

            spec_synth = spec.copy()
            spec_synth.stats.channel = 'Synth'
            spec_synth.data = spectral_model(xdata, *params_opt)
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

    ssp_exit()


if __name__ == '__main__':
    try:
        thismodule=os.path.basename(__file__).replace('.py','')
        this=__import__(thismodule)
        main()
    except KeyboardInterrupt:
        ssp_exit(1)
