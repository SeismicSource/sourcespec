# -*- coding: utf-8 -*-
"""
Build spectral objects.

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
import numpy as np
import math
from scipy.integrate import cumtrapz
from obspy.core import Stream
from obspy.signal.konnoohmachismoothing import (calculate_smoothing_matrix,
                                                apply_smoothing_matrix)

from sourcespec import spectrum
from sourcespec.ssp_setup import dprint
from sourcespec.ssp_util import cosine_taper, moment_to_mag
from sourcespec.ssp_process_traces import filter_trace
from sourcespec.ssp_correction import station_correction

smoothing_matrix = None
smoothing_matrix_weight = None


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
    """
    Compute the component 'H' from geometric mean of the stream components.

    (which can also be all three components)
    """
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
    if spec_h is not None:
        spec_h.data = np.sqrt(spec_h.data)
    return spec_h


def _check_data_len(config, trace):
    traceId = trace.get_id()

    trace_cut = trace.copy()
    t1 = trace.stats.arrivals['S1'][1]
    t2 = trace.stats.arrivals['S2'][1]
    trace_cut.trim(starttime=t1, endtime=t2, pad=True, fill_value=0)
    npts = len(trace_cut.data)
    if npts == 0:
        logging.warning('%s: No data for the selected cut interval: '
                        'skipping trace' % traceId)
        raise RuntimeError
    nzeros = len(np.where(trace_cut.data == 0)[0])
    if nzeros > npts/4:
        logging.warning('%s: Too many gaps for the selected cut '
                        'interval: skipping trace' % traceId)
        raise RuntimeError


def _cut_signal_noise(config, trace):
    trace_signal = trace.copy()
    trace_noise = trace.copy()

    # Integrate in time domain, if required.
    # (otherwhise frequency-domain integration is
    # performed later)
    if config.time_domain_int:
        _time_integrate(config, trace_signal)
        _time_integrate(config, trace_noise)

    # trim...
    t1 = trace.stats.arrivals['S1'][1]
    t2 = trace.stats.arrivals['S2'][1]
    trace_signal.trim(starttime=t1, endtime=t2, pad=True, fill_value=0)
    # Noise time window for weighting function:
    noise_t1 = trace.stats.arrivals['N1'][1]
    noise_t2 = trace.stats.arrivals['N2'][1]
    trace_noise.trim(starttime=noise_t1, endtime=noise_t2, pad=True,
                     fill_value=0)
    # ...taper...
    cosine_taper(trace_signal.data, width=config.taper_halfwidth)
    cosine_taper(trace_noise.data, width=config.taper_halfwidth)
    if config.spectral_win_length is not None:
        # ...and zero pad to spectral_win_length
        trace_signal.trim(starttime=t1,
                          endtime=t1+config.spectral_win_length,
                          pad=True,
                          fill_value=0)
        trace_noise.trim(starttime=noise_t1,
                         endtime=noise_t1+config.spectral_win_length,
                         pad=True,
                         fill_value=0)

    return trace_signal, trace_noise


def _check_noise_level(trace_signal, trace_noise):
    traceId = trace_signal.get_id()
    trace_signal_rms = ((trace_signal.data**2).sum())**0.5
    trace_noise_rms = ((trace_noise.data**2).sum())**0.5
    if trace_noise_rms/trace_signal_rms < 1e-6:
        logging.warning('%s: Noise level is too low or zero: '
                        'ignoring for noise weighting' % traceId)
        raise RuntimeError


def _build_spectrum(config, trace):
    stats = trace.stats

    # normalization for the hypocentral distance
    trace.data *= trace.stats.hypo_dist * 1000

    # calculate fft
    spec = spectrum.do_spectrum(trace)
    spec.stats.instrtype = stats.instrtype
    spec.stats.coords = stats.coords
    spec.stats.hypo = stats.hypo
    spec.stats.hypo_dist = stats.hypo_dist
    spec.stats.epi_dist = stats.epi_dist

    # Integrate in frequency domain, if no time-domain
    # integration has been performed
    if not config.time_domain_int:
        _frequency_integrate(config, spec)

    # smooth the abs of fft
    global smoothing_matrix
    compute_smoothing_matrix = False
    if smoothing_matrix is None:
        logging.info('Computing Konno-Ohmachi smoothing matrix for spectra. '
                     'This might take a while but is done only once.')
        compute_smoothing_matrix = True
    elif spec.data.shape[0] != smoothing_matrix.shape[0]:
        logging.info('Re-computing Konno-Ohmachi smoothing matrix for spectra '
                     'with different shape.')
        compute_smoothing_matrix = True
    if compute_smoothing_matrix:
        smoothing_matrix = calculate_smoothing_matrix(
            spec.get_freq(), bandwidth=40, normalize=True)
    spec.data = apply_smoothing_matrix(spec.data, smoothing_matrix)

    # TODO: parameterize
    # coefficient for converting displ spectrum
    # to seismic moment (Aki&Richards,1980)
    vs_m = config.vs*1000
    vs3 = pow(vs_m, 3)
    coeff = 4 * math.pi * vs3 * config.rho / config.rps / 2.
    dprint('coeff= %f' % coeff)

    # convert to seismic moment
    spec.data *= coeff

    # Cut the spectrum
    spec_cut = _cut_spectrum(config, spec)
    return spec_cut


def _build_weight(spec, specnoise):
    weight = spec.copy()
    if specnoise is not None:
        weight.data /= specnoise.data
        # save data to raw_data
        weight.data_raw = weight.data.copy()
        # The inversion is done in magnitude units,
        # so let's take log10 of weight
        weight.data = np.log10(weight.data)
        # Make sure weight is positive,
        # i.e put weight to zero when S/N < 1
        weight.data[weight.data < 0] = 0
        # smooth weight
        global smoothing_matrix_weight
        compute_smoothing_matrix = False
        if smoothing_matrix_weight is None:
            logging.info('Computing Konno-Ohmachi smoothing matrix for '
                         'weights. This might take a while but is done '
                         'only once.')
            compute_smoothing_matrix = True
        elif weight.data.shape[0] != smoothing_matrix_weight.shape[0]:
            logging.info('Re-computing Konno-Ohmachi smoothing matrix for '
                         'weights with different shape.')
            compute_smoothing_matrix = True
        if compute_smoothing_matrix:
            smoothing_matrix_weight = calculate_smoothing_matrix(
                spec.get_freq(), bandwidth=20, normalize=True)
        weight.data = apply_smoothing_matrix(
            weight.data, smoothing_matrix_weight)
        # normalization
        weight.data /= np.max(weight.data)
    else:
        logging.warning('%s: No available noise window: '
                        % weight.get_id()[0:-1] +
                        'a uniform weight will be applied')
        weight.data = np.ones(len(spec.data))
    return weight


def _build_H_and_weight(spec_st, specnoise_st):
    """
    Add to spec_st the "H" component.

    H component is obtained from the modulus of all the available components.

    The same for noise, if requested. In this case we compute
    weighting function as well.
    """
    if specnoise_st:
        noise_weight = True
    else:
        noise_weight = False
    weight_st = Stream()
    for station in set(x.stats.station for x in spec_st.sort()):
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
                if specnoise_h is not None:
                    specnoise_st.append(specnoise_h)

                # Weighting function is the ratio between "H" components
                # of signal and noise
                weight = _build_weight(spec_h, specnoise_h)
                weight_st.append(weight)
    return weight_st


def build_spectra(config, st, noise_weight=False):
    """
    Build spectra and the spec_st object.

    Computes S-wave (displacement) spectra from
    accelerometers and velocimeters, uncorrected for attenuation,
    corrected for instrumental constants, normalized by
    hypocentral distance.
    """
    spec_st = Stream()
    specnoise_st = Stream()

    for trace in st.sort():
        try:
            _check_data_len(config, trace)
            trace_signal, trace_noise = _cut_signal_noise(config, trace)
            _check_noise_level(trace_signal, trace_noise)
        except (ValueError, RuntimeError):
            continue
        spec = _build_spectrum(config, trace_signal)
        if noise_weight:
            specnoise = _build_spectrum(config, trace_noise)
            weight = _build_weight(spec, specnoise)
            if config.spectral_sn_freq_range is not None:
                sn_fmin, sn_fmax = config.spectral_sn_freq_range
                freqs = weight.get_freq()
                idx = np.where((sn_fmin <= freqs)*(freqs <= sn_fmax))
            else:
                idx = range(len(weight.data_raw))
            spectral_ssn =\
                weight.data_raw[idx].sum()/len(weight.data_raw[idx])
            logging.info('%s: spectral S/N: %.2f' %
                         (spec.get_id(), spectral_ssn))
            if config.spectral_sn_min:
                ssnmin = config.spectral_sn_min
                if spectral_ssn < ssnmin:
                    logging.warning('%s: spectral S/N smaller than %.2f: '
                                    'skipping spectrum' %
                                    (spec.get_id(), ssnmin))
                    continue
            specnoise_st.append(specnoise)
        spec_st.append(spec)

    # build H component and weight_st
    weight_st = _build_H_and_weight(spec_st, specnoise_st)

    # convert the spectral amplitudes to moment magnitude
    for spec in spec_st:
        spec.data_mag = moment_to_mag(spec.data)

    # optionally, apply station correction
    if config.options.correction:
        spec_st = station_correction(spec_st, config)

    if noise_weight:
        for specnoise in specnoise_st:
            specnoise.data_mag = moment_to_mag(specnoise.data)
        return spec_st, specnoise_st, weight_st
    else:
        return spec_st
