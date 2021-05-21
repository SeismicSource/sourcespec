# -*- coding: utf-8 -*-
"""
Build spectral objects.

:copyright:
    2012 Claudio Satriano <satriano@ipgp.fr>

    2013-2014 Claudio Satriano <satriano@ipgp.fr>,
              Emanuela Matrullo <matrullo@geologie.ens.fr>,
              Agnes Chounet <chounet@ipgp.fr>

    2015-2021 Claudio Satriano <satriano@ipgp.fr>
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
from scipy.interpolate import interp1d
from obspy.core import Stream
from sourcespec import spectrum
from sourcespec.ssp_util import smooth, cosine_taper, moment_to_mag, get_vel
from sourcespec.ssp_process_traces import filter_trace
from sourcespec.ssp_correction import station_correction
logger = logging.getLogger(__name__.split('.')[-1])


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
    # see if there is a station-specfic frequency range
    station = spec.stats.station
    try:
        freq1 = float(config['freq1_' + station])
        freq2 = float(config['freq2_' + station])
    except KeyError:
        try:
            instrtype = spec.stats.instrtype
            freq1 = float(config['freq1_' + instrtype])
            freq2 = float(config['freq2_' + instrtype])
        except KeyError:
            logger.warning('%s: Unknown instrument type: %s: '
                           'skipping spectrum' % (spec.id, instrtype))
            raise ValueError
    return spec.slice(freq1, freq2)


def _compute_h(spec_st, code, wave_type='S'):
    """
    Compute the component 'H' from geometric mean of the stream components.

    (which can also be all three components)
    """
    spec_h = None
    for spec in spec_st.traces:
        # this avoids taking a component from co-located station:
        # ('code' is band+instrument code)
        channel = spec.stats.channel
        if channel[:-1] != code:
            continue
        # avoid reusing a previous 'H' channel
        if channel[-1] == 'H':
            continue
        # only use transverse component for SH
        if wave_type == 'SH' and channel[-1] != 'T':
            continue
        # only use radial and, optionally, vertical component for SV
        if wave_type == 'SV' and channel[-1] == 'R':
            continue
        if spec_h is None:
            spec_h = spec.copy()
            spec_h.data = np.power(spec_h.data, 2)
            spec_h.data_log = np.power(spec_h.data_log, 2)
            spec_h.stats.channel = code + 'H'
        else:
            spec_h.data += np.power(spec.data, 2)
            spec_h.data_log += np.power(spec.data_log, 2)
    if spec_h is not None:
        spec_h.data = np.sqrt(spec_h.data)
        spec_h.data_log = np.sqrt(spec_h.data_log)
    return spec_h


def _check_data_len(config, trace):
    traceId = trace.get_id()

    trace_cut = trace.copy()
    t1 = trace.stats.arrivals['S1'][1]
    t2 = trace.stats.arrivals['S2'][1]
    trace_cut.trim(starttime=t1, endtime=t2, pad=True, fill_value=0)
    npts = len(trace_cut.data)
    if npts == 0:
        logger.warning('%s: No data for the selected cut interval: '
                       'skipping trace' % traceId)
        raise RuntimeError
    nzeros = len(np.where(trace_cut.data == 0)[0])
    if nzeros > npts/4:
        logger.warning('%s: Too many gaps for the selected cut '
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
    # Be sure that both traces have same length:
    if len(trace_signal) != len(trace_noise):
        npts = min(len(trace_signal), len(trace_noise))
        trace_signal.data = trace_signal.data[:npts]
        trace_noise.data = trace_noise.data[:npts]

    return trace_signal, trace_noise


def _check_noise_level(trace_signal, trace_noise):
    traceId = trace_signal.get_id()
    trace_signal_rms = ((trace_signal.data**2).sum())**0.5
    trace_noise_rms = ((trace_noise.data**2).sum())**0.5
    if trace_noise_rms/trace_signal_rms < 1e-6:
        logger.warning('%s: Noise level is too low or zero: '
                       'ignoring for noise weighting' % traceId)
        raise RuntimeError


def _displacement_to_moment(stats, config):
    """
    Return the coefficient for converting displacement to seismic moment.

    From Aki&Richards,1980
    """
    vs_hypo = config.hypo.vs
    vs_station = get_vel(
        stats.coords.longitude, stats.coords.latitude, -stats.coords.elevation,
        'S', config.NLL_model_dir)
    if vs_station is None:
        vs_station = config.vs
    vs_hypo *= 1000.
    vs_station *= 1000.
    vs3 = vs_hypo**(5./2) * vs_station**(1./2)
    return 4 * math.pi * vs3 * config.rho / (2 * config.rps)


def _smooth_spectrum(spec, npts=5):
    """Smooth spectrum in a log_freq space."""
    freq = spec.get_freq()
    f = interp1d(freq, spec.data, fill_value='extrapolate')
    freq_log = np.logspace(np.log10(freq[0]), np.log10(freq[-1]))
    spec.freq_log = freq_log
    spec.data_log = f(freq_log)
    spec.data_log = smooth(spec.data_log, window_len=npts)
    f = interp1d(freq_log, spec.data_log, fill_value='extrapolate')
    spec.data = f(freq)


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
    spec.stats.ignore = stats.ignore
    # Integrate in frequency domain, if no time-domain
    # integration has been performed
    if not config.time_domain_int:
        _frequency_integrate(config, spec)
    # cut the spectrum
    spec = _cut_spectrum(config, spec)
    # convert to seismic moment
    coeff = _displacement_to_moment(stats, config)
    spec.data *= coeff
    # smooth
    _smooth_spectrum(spec, npts=5)
    return spec


def _build_weight(spec, specnoise):
    weight = spec.copy()
    if specnoise is not None:
        weight.data /= specnoise.data
        # save data to raw_data
        weight.data_raw = weight.data.copy()
        # The inversion is done in magnitude units,
        # so let's take log10 of weight
        weight.data = np.log10(weight.data)
        # Make sure weight is positive
        weight.data[weight.data <= 0] = 0.001
        _smooth_spectrum(weight, npts=11)
        weight.data /= np.max(weight.data)
    else:
        logger.warning('%s: No available noise window: '
                       % weight.get_id()[0:-1] +
                       'a uniform weight will be applied')
        weight.data = np.ones(len(spec.data))
    # interpolate to log-frequencies
    f = interp1d(weight.get_freq(), weight.data, fill_value='extrapolate')
    weight.data_log = f(weight.freq_log)
    return weight


def _build_H_and_weight(spec_st, specnoise_st, wave_type='S'):
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
    stalist = set(sp.id[:-1] for sp in spec_st if not sp.stats.ignore)
    for specid in stalist:
        network, station, location, code = specid.split('.')
        spec_st_sel = spec_st.select(
            network=network, station=station, location=location)
        spec_st_sel = Stream(sp for sp in spec_st_sel if not sp.stats.ignore)
        if noise_weight:
            specnoise_st_sel = specnoise_st.select(
                network=network, station=station, location=location)
            specnoise_st_sel = Stream(
                sp for sp in specnoise_st_sel if not sp.stats.ignore)
        # 'code' is band+instrument code
        for code in set(x.stats.channel[:-1] for x in spec_st_sel):
            spec_h = _compute_h(spec_st_sel, code, wave_type)
            if spec_h is None:
                continue
            spec_st.append(spec_h)

            # Compute "H" component for noise, if requested,
            # and weighting function.
            if noise_weight:
                specnoise_h = _compute_h(specnoise_st_sel, code, wave_type)
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

    # sort by sampling rate: this limits the number of times on which
    # konno-ohmachi smoothing matrix is recomputed
    for trace in st.sort(keys=['sampling_rate', 'station']):
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
            spectral_snratio =\
                weight.data_raw[idx].sum()/len(weight.data_raw[idx])
            spec.stats.spectral_snratio = spectral_snratio
            logger.info('%s: spectral S/N: %.2f' %
                        (spec.get_id(), spectral_snratio))
            if config.spectral_sn_min:
                ssnmin = config.spectral_sn_min
                if spectral_snratio < ssnmin:
                    logger.warning('%s: spectral S/N smaller than %.2f: '
                                   'ignoring spectrum' %
                                   (spec.get_id(), ssnmin))
                    trace.stats.ignore = True
                    spec.stats.ignore = True
                    specnoise.stats.ignore = True
            specnoise_st.append(specnoise)
        spec_st.append(spec)

    # build H component and weight_st
    weight_st = _build_H_and_weight(spec_st, specnoise_st, config.wave_type)

    # convert the spectral amplitudes to moment magnitude
    for spec in spec_st:
        spec.data_mag = moment_to_mag(spec.data)
        spec.data_log_mag = moment_to_mag(spec.data_log)

    # optionally, apply station correction
    if config.options.correction:
        spec_st = station_correction(spec_st, config)

    if noise_weight:
        for specnoise in specnoise_st:
            specnoise.data_mag = moment_to_mag(specnoise.data)
        return spec_st, specnoise_st, weight_st
    else:
        return spec_st
