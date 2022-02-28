# -*- coding: utf-8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Build spectral objects.

:copyright:
    2012 Claudio Satriano <satriano@ipgp.fr>

    2013-2014 Claudio Satriano <satriano@ipgp.fr>,
              Emanuela Matrullo <matrullo@geologie.ens.fr>,
              Agnes Chounet <chounet@ipgp.fr>

    2015-2022 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
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
from sourcespec.ssp_setup import ssp_exit
from sourcespec.ssp_util import smooth, cosine_taper, moment_to_mag, get_vel
from sourcespec.ssp_process_traces import filter_trace
from sourcespec.ssp_correction import station_correction
from sourcespec.ssp_radiation_pattern import get_radiation_pattern_coefficient
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
            msg = '{}: Unknown instrument type: {}: skipping spectrum'
            msg = msg.format(spec.id, instrtype)
            raise RuntimeError(msg)
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
        if wave_type == 'SV' and channel[-1] == 'T':
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
        msg = '{}: No data for the selected cut interval: skipping trace'
        msg = msg.format(traceId)
        raise RuntimeError(msg)
    nzeros = len(np.where(trace_cut.data == 0)[0])
    if nzeros > npts/4:
        msg = '{}: S-wave window is zero for more than 25%: skipping trace'
        msg = msg.format(traceId)
        raise RuntimeError(msg)


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
        msg =\
            '{}: Noise level is too low or zero: ignoring for noise weighting'
        msg = msg.format(traceId)
        raise RuntimeError(msg)


def _displacement_to_moment(stats, config):
    """
    Return the coefficient for converting displacement to seismic moment.

    From Aki&Richards,1980
    """
    vs_hypo = config.hypo.vs
    vs_station = get_vel(
        stats.coords.longitude, stats.coords.latitude, -stats.coords.elevation,
        'S', config)
    vs_hypo *= 1000.
    vs_station *= 1000.
    vs3 = vs_hypo**(5./2) * vs_station**(1./2)
    rps = get_radiation_pattern_coefficient(stats, config)
    coeff = 4 * math.pi * vs3 * config.rho / (2 * rps)
    return coeff


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
    # store coeff to correct back data in displacement units
    # for radiated_energy()
    spec.coeff = coeff
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
        _smooth_spectrum(weight, npts=11)
        weight.data /= np.max(weight.data)
        # slightly taper weight at low frequencies, to avoid overestimating
        # weight at low frequencies, in cases where noise is underestimated
        cosine_taper(weight.data, spec.stats.delta/4, left_taper=True)
        # Make sure weight is positive
        weight.data[weight.data <= 0] = 0.001
    else:
        msg = '{}: No available noise window: a uniform weight will be applied'
        msg = msg.format(weight.get_id()[0:-1])
        logger.warning(msg)
        weight.data = np.ones(len(spec.data))
    # interpolate to log-frequencies
    f = interp1d(weight.get_freq(), weight.data, fill_value='extrapolate')
    weight.data_log = f(weight.freq_log)
    return weight


def _build_weight_st(spec_st, specnoise_st):
    """Build the weight spectrum."""
    weight_st = Stream()
    spec_ids = set(sp.id[:-1] for sp in spec_st if not sp.stats.ignore)
    for specid in spec_ids:
        try:
            spec_h = _select_spectra(spec_st, specid + 'H')[0]
            specnoise_h = _select_spectra(specnoise_st, specid + 'H')[0]
        except Exception:
            continue
        weight = _build_weight(spec_h, specnoise_h)
        weight_st.append(weight)
    return weight_st


def _select_spectra(spec_st, specid):
    """Select spectra from stream, based on specid."""
    network, station, location, channel = specid.split('.')
    channel = channel + '?'*(3-len(channel))
    spec_st_sel = spec_st.select(
        network=network, station=station, location=location, channel=channel)
    spec_st_sel = Stream(sp for sp in spec_st_sel if not sp.stats.ignore)
    return spec_st_sel


def _build_H(spec_st, specnoise_st=None, wave_type='S'):
    """
    Add to spec_st and specnoise_st the "H" component.

    H component is obtained from the modulus of all the available components.
    """
    spec_ids = set(sp.id[:-1] for sp in spec_st if not sp.stats.ignore)
    for specid in spec_ids:
        spec_st_sel = _select_spectra(spec_st, specid)
        specnoise_st_sel = _select_spectra(specnoise_st, specid)
        # 'code' is band+instrument code
        for code in set(x.stats.channel[:-1] for x in spec_st_sel):
            spec_h = _compute_h(spec_st_sel, code, wave_type)
            if spec_h is not None:
                spec_st.append(spec_h)
            specnoise_h = _compute_h(specnoise_st_sel, code, wave_type)
            if specnoise_h is not None:
                specnoise_st.append(specnoise_h)


def _check_spectral_sn_ratio(config, spec, specnoise):
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
    logger.info(
        '{}: spectral S/N: {:.2f}'.format(
            spec.get_id(), spectral_snratio))
    if config.spectral_sn_min:
        ssnmin = config.spectral_sn_min
        if spectral_snratio < ssnmin:
            msg = '{}: spectral S/N smaller than {:.2f}: '
            msg += 'ignoring spectrum'
            msg = msg.format(spec.get_id(), ssnmin)
            raise ValueError(msg)


def build_spectra(config, st):
    """
    Build spectra and the spec_st object.

    Computes S-wave (displacement) spectra from
    accelerometers and velocimeters, uncorrected for attenuation,
    corrected for instrumental constants, normalized by
    hypocentral distance.
    """
    logger.info('Building spectra...')
    spec_st = Stream()
    specnoise_st = Stream()

    # sort by trace id
    for trace in sorted(st, key=lambda tr: tr.id):
        try:
            _check_data_len(config, trace)
            trace_signal, trace_noise = _cut_signal_noise(config, trace)
            _check_noise_level(trace_signal, trace_noise)
            spec = _build_spectrum(config, trace_signal)
            specnoise = _build_spectrum(config, trace_noise)
            _check_spectral_sn_ratio(config, spec, specnoise)
        except RuntimeError as msg:
            # RuntimeError is for skipped spectra
            logger.warning(msg)
            continue
        except ValueError as msg:
            # ValueError is for ignored spectra, which are still stored
            logger.warning(msg)
            trace.stats.ignore = True
            spec.stats.ignore = True
            specnoise.stats.ignore = True
        spec_st.append(spec)
        specnoise_st.append(specnoise)

    if not spec_st:
        logger.error('No spectra left! Exiting.')
        ssp_exit()

    # build H component
    _build_H(spec_st, specnoise_st, config.wave_type)

    # convert the spectral amplitudes to moment magnitude
    for spec in spec_st:
        spec.data_mag = moment_to_mag(spec.data)
        spec.data_log_mag = moment_to_mag(spec.data_log)

    # apply station correction if a residual file is specified in config
    spec_st = station_correction(spec_st, config)

    # build the weight spectrum
    weight_st = _build_weight_st(spec_st, specnoise_st)

    logger.info('Building spectra: done')
    if config.weighting == 'noise':
        for specnoise in specnoise_st:
            specnoise.data_mag = moment_to_mag(specnoise.data)
        return spec_st, specnoise_st, weight_st
    else:
        return spec_st
