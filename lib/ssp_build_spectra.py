# -*- coding: utf-8 -*-
# ssp_process_traces.py
#
# Builds spectra for source_spec 
# (c) 2012 Claudio Satriano <satriano@ipgp.fr>
# (c) 2013 Claudio Satriano <satriano@ipgp.fr>,
#          Emanuela Matrullo <matrullo@geologie.ens.fr>
import logging
import numpy as np
import math
from scipy.integrate import cumtrapz
from copy import deepcopy, copy
from obspy.core import Stream
from ssp_setup import dprint
from ssp_util import swave_arrival, cosine_taper, moment_to_mag
from lib import spectrum
from obspy.signal.konnoohmachismoothing import konnoOhmachiSmoothing

def build_spectra(config, st):
    ''' Builds spectra and the spec_st object '''
    spec_st = Stream()
    for trace in st.traces:
        traceId = trace.getId()
        stats = trace.stats

        # S time window
        s_arrival_time = swave_arrival(trace, config.vs)
        t1 = s_arrival_time - config.pre_s_time
        t2 = t1 + config.s_win_length

        trace_cut = copy(trace)
        trace_cut.stats = deepcopy(trace.stats)

        # Do a preliminary trim, in order to check if there is enough
        # data within the selected time window
        trace_cut.trim(starttime=t1, endtime=t2, pad=True, fill_value=0)
        npts = len(trace_cut.data)
        if npts == 0:
            logging.warning('%s: No data for the selected cut interval: skipping trace' % traceId)
            continue
        nzeros = len(np.where(trace_cut.data==0)[0])
        if nzeros > npts/4:
            logging.warning('%s: Too many gaps for the selected cut interval: skipping trace' % traceId)
            continue
        
        # If the check is ok, recover the full trace
        # (it will be cutted later on)
        trace_cut = copy(trace)
        trace_cut.stats = deepcopy(trace.stats)

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

        if config.time_domain_int:
            for i in range(0,nint):
                # remove the mean...
                trace_cut.detrend(type='constant')
                # ...and the linear trend...
                trace_cut.detrend(type='linear')
                trace_cut.filter(type='bandpass', freqmin=bp_freqmin, freqmax=bp_freqmax)
                # integrate
                trace_cut.data = cumtrapz(trace_cut.data) * trace_cut.stats.delta
                trace_cut.stats.npts -= 1
        else:
            # remove the mean...
            trace_cut.detrend(type='constant')
            # ...and the linear trend...
            trace_cut.detrend(type='linear')
            trace_cut.filter(type='bandpass', freqmin=bp_freqmin, freqmax=bp_freqmax)

        # trim...
        trace_cut.trim(starttime=t1, endtime=t2, pad=True, fill_value=0)
        # ...and taper
        cosine_taper(trace_cut.data, width=config.taper_halfwidth)

        # normalization for the hypocentral distance
        trace_cut.data *= trace_cut.stats.hypo_dist * 1000

        # calculate fft
        spec = spectrum.do_spectrum(trace_cut)
        spec.stats.instrtype = instrtype
        spec.stats.coords    = stats.coords
        spec.stats.hypo      = stats.hypo
        spec.stats.hypo_dist = stats.hypo_dist

        if not config.time_domain_int:
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

        # TODO: parameterize
        # coefficient for converting displ spectrum
        # to seismic moment (Aki&Richards,1980)
        vs_m = config.vs*1000
        vs3  = pow(vs_m,3) 
        coeff = 4 * math.pi * vs3 * config.rho /config.rps/2.
        dprint('coeff= %f' % coeff)

        # convert to seismic moment
        spec.data *= coeff

        # Cut the spectrum between freq1 and freq2
        spec_cut = spec.slice(freq1, freq2)
        spec_st.append(spec_cut)

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
                    #data_h = np.sqrt(data_h / 2) #divide by the number of horizontal components
                    data_h = np.sqrt(data_h)
                    spec_h.data = data_h
            spec_st.append(spec_h)

    # convert the spectral amplitudes to moment magnitude
    for spec in spec_st:
        spec.data_mag = moment_to_mag(spec.data)

    return spec_st
