# -*- coding: utf-8 -*-
# ssp_process_traces.py
#
# (c) 2012 Claudio Satriano <satriano@ipgp.fr>
# (c) 2013-2014 Claudio Satriano <satriano@ipgp.fr>,
#               Emanuela Matrullo <matrullo@geologie.ens.fr>,
#               Agnes Chounet <chounet@ipgp.fr>
# (c) 2015 Claudio Satriano <satriano@ipgp.fr>
'''
Build spectral objects.
'''
from __future__ import division
import logging
import numpy as np
import math
from scipy.integrate import cumtrapz
from copy import deepcopy, copy
from obspy.core import Stream
from ssp_setup import dprint
from ssp_util import wave_arrival, cosine_taper,\
        moment_to_mag
from ssp_process_traces import filter_trace
from ssp_correction import station_correction
import spectrum
from obspy.signal.konnoohmachismoothing import konnoOhmachiSmoothing


def _time_integrate(config, trace):
    instrtype = trace.stats.instrtype
    if instrtype == 'acc':
        nint = 2
    elif instrtype == 'shortp':
        nint = 1
    elif instrtype == 'broadb':
        nint = 1
    else:
        raise ValueError
    trace.detrend(type='constant')
    trace.detrend(type='linear')
    for i in range(0, nint):
        trace.data = cumtrapz(trace.data) * trace.stats.delta
        trace.stats.npts -= 1
        filter_trace(config, trace)


def _frequency_integrate(config, spec):
    instrtype = spec.stats.instrtype
    if instrtype == 'acc':
        nint = 2
    elif instrtype == 'shortp':
        nint = 1
    elif instrtype == 'broadb':
        nint = 1
    else:
        raise ValueError
    for i in range(0, nint):
        spec.data /= (2 * math.pi * spec.get_freq())


def _cut_spectrum(config, spec):
    instrtype = spec.stats.instrtype
    if instrtype == 'acc':
        freq1 = config.freq1_acc
        freq2 = config.freq2_acc
    elif instrtype == 'shortp':
        freq1 = config.freq1_shortp
        freq2 = config.freq2_shortp
    elif instrtype == 'broadb':
        freq1 = config.freq1_broadb
        freq2 = config.freq2_broadb
    else:
        raise ValueError
    return spec.slice(freq1, freq2)


def _compute_h(spec_st, code):
    '''
    Computes the component 'H' from geometric
    mean of the stream components
    (which can also be all three components)
    '''
    spec_h = None
    for spec in spec_st.traces:
        # this avoids taking a component from co-located station:
        # ('code' is band+instrument code)
        if spec.stats.channel[0:2] != code:
            continue
        if spec_h is None:
            spec_h = spec.copy()
            spec_h.data = np.power(spec_h.data, 2)
            spec_h.stats.channel = code + 'H'
        else:
            spec_h.data += np.power(spec.data, 2)
    spec_h.data = np.sqrt(spec_h.data)
    return spec_h


def build_spectra(config, st, noise_weight=False):
    '''
    Builds spectra and the spec_st object.
    Computes S-wave (displacement) spectra from
    accelerometers and velocimeters, uncorrected for attenuation,
    corrected for instrumental constants, normalized by
    hypocentral distance.
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
        trace.stats.arrivals['S1'] = ('S1', t1 - trace.stats.starttime)
        trace.stats.arrivals['S2'] = ('S2', t2 - trace.stats.starttime)

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

        # Integrate in time domain, if required.
        # (otherwhise frequency-domain integration is
        # performed later)
        if config.time_domain_int:
            _time_integrate(config, trace_cut)

        if noise_weight:
            # Noise time window for weighting function:
            p_arrival_time = wave_arrival(trace, config.vp, 'P')
            noise_start_t = p_arrival_time - config.pre_p_time
            noise_end_t = noise_start_t + config.s_win_length
            trace.stats.arrivals['N1'] = ('N1', noise_start_t - trace.stats.starttime)
            trace.stats.arrivals['N2'] = ('N2', noise_end_t - trace.stats.starttime)
            # We inherit the same processing of trace_cut
            # (which, despite its name, has not been yet cut!)
            trace_noise = copy(trace_cut)
            trace_noise.stats = deepcopy(trace_cut.stats)

        # trim...
        trace_cut.trim(starttime=t1, endtime=t2, pad=True, fill_value=0)
        # ...taper...
        cosine_taper(trace_cut.data, width=config.taper_halfwidth)
        if not math.isnan(config.spectral_win_length):
            # ...and zero pad to spectral_win_length
            trace_cut.trim(
                    starttime=t1,
                    endtime=t1+config.spectral_win_length,
                    pad=True,
                    fill_value=0
                    )

        if noise_weight:
            # ...trim...
            trace_noise.trim(starttime=noise_start_t, endtime=noise_end_t, pad=True, fill_value=0)
            # ...taper...
            cosine_taper(trace_noise.data, width=config.taper_halfwidth)
            if not math.isnan(config.spectral_win_length):
                # ...and zero pad to spectral_win_length
                trace_noise.trim(
                        starttime=noise_start_t,
                        endtime=noise_start_t+config.spectral_win_length,
                        pad=True,
                        fill_value=0
                        )

        # normalization for the hypocentral distance
        trace_cut.data *= trace_cut.stats.hypo_dist * 1000

        # calculate fft
        spec = spectrum.do_spectrum(trace_cut)
        spec.stats.instrtype = stats.instrtype
        spec.stats.coords = stats.coords
        spec.stats.hypo = stats.hypo
        spec.stats.hypo_dist = stats.hypo_dist

        # same processing for noise, if requested
        if noise_weight:
            cosine_taper(trace_noise.data, width=config.taper_halfwidth)
            trace_noise.data *= trace_noise.stats.hypo_dist * 1000
            specnoise = spectrum.do_spectrum(trace_noise)
            specnoise.stats.instrtype = stats.instrtype
            specnoise.stats.coords = stats.coords
            specnoise.stats.hypo = stats.hypo
            specnoise.stats.hypo_dist = stats.hypo_dist

        # Integrate in frequency domain, if no time-domain
        # integration has been performed
        if not config.time_domain_int:
            _frequency_integrate(config, spec)
            if noise_weight:
                _frequency_integrate(config, specnoise)

        # smooth the abs of fft
        spec.data = konnoOhmachiSmoothing(spec.data, spec.get_freq(),
                                          40, normalize=True)

        # TODO: parameterize
        # coefficient for converting displ spectrum
        # to seismic moment (Aki&Richards,1980)
        vs_m = config.vs*1000
        vs3 = pow(vs_m,3)
        coeff = 4 * math.pi * vs3 * config.rho /config.rps/2.
        dprint('coeff= %f' % coeff)

        # convert to seismic moment
        spec.data *= coeff

        # same processing for noise, if requested
        if noise_weight:
            specnoise.data = konnoOhmachiSmoothing(specnoise.data,
                                          specnoise.get_freq(),
                                          40, normalize=True)
            specnoise.data *= coeff

        # Cut the spectrum
        spec_cut = _cut_spectrum(config, spec)
        # append to spec streams
        spec_st.append(spec_cut)

        if noise_weight:
            specnoise_cut = _cut_spectrum(config, specnoise)
            specnoise_st.append(specnoise_cut)
    # end of loop on traces

    # Add to spec_st the "H" component, obtained from the
    # modulus of all the available components.
    # The same for noise, if requested. In this case we compute
    # weighting function as well.
    for station in set(x.stats.station for x in spec_st.traces):
        spec_st_sel = spec_st.select(station=station)
        if noise_weight:
            specnoise_st_sel = specnoise_st.select(station=station)
        # 'code' is band+instrument code
        for code in set(x.stats.channel[0:2] for x in spec_st_sel):
            spec_h = _compute_h(spec_st_sel, code)
            spec_st.append(spec_h)

            # Compute "H" component for noise, if requested,
            # and weighting function.
            if noise_weight:
                specnoise_h = _compute_h(specnoise_st_sel, code)
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

    # optionally, apply station correction
    if config.options.correction:
        spec_st = station_correction(spec_st, config)

    if noise_weight:
        for specnoise in specnoise_st:
            specnoise.data_mag = moment_to_mag(specnoise.data)

    if noise_weight:
        return spec_st, specnoise_st, weight_st
    else:
        return spec_st
