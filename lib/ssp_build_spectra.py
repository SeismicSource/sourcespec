# -*- coding: utf-8 -*-
# ssp_process_traces.py
#
# Builds spectra for source_spec 
# (c) 2012 Claudio Satriano <satriano@ipgp.fr>
# (c) 2013 Claudio Satriano <satriano@ipgp.fr>,
#          Emanuela Matrullo <matrullo@geologie.ens.fr>
from __future__ import division
import logging
import numpy as np
import math
from scipy.integrate import cumtrapz
from copy import deepcopy, copy
from obspy.core import Stream
from ssp_setup import dprint
from ssp_util import wave_arrival, cosine_taper,\
        moment_to_mag, select_trace
import spectrum
from obspy.signal.konnoohmachismoothing import konnoOhmachiSmoothing


def __process__(trace, bp_freqmin, bp_freqmax, integrate=False):
    # remove the mean...
    trace.detrend(type='constant')
    # ...and the linear trend...
    trace.detrend(type='linear')
    trace.filter(type='bandpass', freqmin=bp_freqmin, freqmax=bp_freqmax)
    if integrate:
        trace.data = cumtrapz(trace.data) * trace.stats.delta
        trace.stats.npts -= 1


def __compute_h__(spec_st, instrtype):
    spec_h = None
    for spec in spec_st.traces:
        # this avoids taking a component from co-located station:
        if spec.stats.instrtype != instrtype:
            continue
        if spec_h == None:
            spec_h = spec.copy()
            spec_h.data = np.power(spec_h.data, 2)
            spec_h.stats.channel = 'H'
        else:
            spec_h.data += np.power(spec.data, 2)
    spec_h.data = np.sqrt(spec_h.data)
    return spec_h
        

def build_spectra(config, st, noise_st=None):
    ''' 
    Builds spectra and the spec_st object.
    Computes S-wave (displacement) spectra from
    accelerometers and velocimeters, uncorrected for attenuation,
    corrected for instrumental constants, normalized by
    hypocentral distance
    '''
    spec_st = Stream()
    specnoise_st = Stream()
    weight_st = Stream()

    for trace in st.traces:
        traceId = trace.getId()
        stats = trace.stats

        # S time window
        s_arrival_time = wave_arrival(trace, config.vs, 'S')
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
            n_time_procs = nint
            n_freq_int = 0
            integrate = True
        else:
            n_time_procs = 1
            n_freq_int = nint
            integrate = False

        # Time domain processing
        for i in range(0, n_time_procs):
            __process__(trace_cut, bp_freqmin, bp_freqmax, integrate)
            if noise_st:
                trace_noise = select_trace(noise_st, traceId, instrtype)
                __process__(trace_noise, bp_freqmin, bp_freqmax, integrate)

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

        # same processing for noise, if requested
        if noise_st:
            cosine_taper(trace_noise.data, width=config.taper_halfwidth)
            trace_noise.data *= trace_noise.stats.hypo_dist * 1000
            specnoise = spectrum.do_spectrum(trace_noise)
            specnoise.stats.instrtype = instrtype
            specnoise.stats.coords    = stats.coords
            specnoise.stats.hypo      = stats.hypo
            specnoise.stats.hypo_dist = stats.hypo_dist

        # Integrate in frequency domain (divide by the pulsation omega)
        # (performed only if config.time_domain_int == False)
        for i in range(0, n_freq_int):
            spec.data /= (2 * math.pi * spec.get_freq())
            if noise_st:
                specnoise.data /= (2 * math.pi * specnoise.get_freq())

        # smooth the abs of fft
        spec.data = konnoOhmachiSmoothing(spec.data, spec.get_freq(),
                                          40, normalize=True)

        # TODO: parameterize
        # coefficient for converting displ spectrum
        # to seismic moment (Aki&Richards,1980)
        vs_m = config.vs*1000
        vs3  = pow(vs_m,3) 
        coeff = 4 * math.pi * vs3 * config.rho /config.rps/2.
        dprint('coeff= %f' % coeff)

        # convert to seismic moment
        spec.data *= coeff

        # same processing for noise, if requested
        if noise_st:
            specnoise.data = konnoOhmachiSmoothing(specnoise.data,
                                          specnoise.get_freq(),
                                          40, normalize=True)
            specnoise.data *= coeff

        # Cut the spectrum between freq1 and freq2
        spec_cut = spec.slice(freq1, freq2)
        # append to spec streams
        spec_st.append(spec_cut)

        if noise_st:
            specnoise_cut = specnoise.slice(freq1, freq2)
            specnoise_st.append(specnoise_cut)
    # end of loop on traces


    # Add to spec_st the "H" component, obtained from the
    # modulus of all the available components.
    # The same for noise, if requested. In this case we compute
    # weighting function as well.
    for station in set(x.stats.station for x in spec_st.traces):
        spec_st_sel = spec_st.select(station=station)
        if noise_st:
            specnoise_st_sel = specnoise_st.select(station=station)
        for instrtype in set(x.stats.instrtype for x in spec_st_sel):
            spec_h = __compute_h__(spec_st_sel, instrtype)
            spec_st.append(spec_h)

            # Compute "H" component for noise, if requested,
            # and weighting function.
            if noise_st:
                specnoise_h = __compute_h__(specnoise_st_sel, instrtype)
                specnoise_st.append(specnoise_h)

                # Weighting function is the ratio between "H" components
                # of signal and noise
                weight = copy(spec_h)
                weight.stats = deepcopy(spec_h.stats)
                weight.data = deepcopy(spec_h.data)
                weight.data /= specnoise_h.data
                # smooth weight
                weight.data = konnoOhmachiSmoothing(weight.data,
                                          weight.get_freq(),
                                          10, normalize=True)
                # normalization
                weight.data /= np.max(weight.data)
                weight_st.append(weight)

    # convert the spectral amplitudes to moment magnitude
    for spec in spec_st:
        spec.data_mag = moment_to_mag(spec.data)
    if noise_st:
        for specnoise in specnoise_st:
            specnoise.data_mag = moment_to_mag(specnoise.data)


    if noise_st:
        return spec_st, specnoise_st, weight_st
    else:
        return spec_st


