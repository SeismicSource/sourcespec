# -*- coding: utf-8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Build spectral objects.

:copyright:
    2012 Claudio Satriano <satriano@ipgp.fr>

    2013-2014 Claudio Satriano <satriano@ipgp.fr>,
              Emanuela Matrullo <matrullo@geologie.ens.fr>,
              Agnes Chounet <chounet@ipgp.fr>

    2015-2026 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
import logging
import math
import fnmatch
import numbers
from collections.abc import Mapping
import numpy as np
# cumtrapz was deprecated in scipy 1.12.0 and removed in 1.14.0
try:
    from scipy.integrate import cumtrapz
except ImportError:
    from scipy.integrate import cumulative_trapezoid as cumtrapz
from obspy.core import Stream
from sourcespec.spectrum import Spectrum, SpectrumStream
from sourcespec.ssp_setup import ssp_exit
from sourcespec.ssp_util import (
    smooth, cosine_taper, moment_to_mag, MediumProperties)
# NOTE: 'smooth' is a time-domain window convolution; it assumes finite values.
from sourcespec.ssp_geom_spreading import (
    geom_spread_r_power_n, geom_spread_r_power_n_segmented,
    geom_spread_boatwright, geom_spread_teleseismic)
from sourcespec.ssp_process_traces import filter_trace
from sourcespec.ssp_correction import station_correction
from sourcespec.ssp_radiation_pattern import get_radiation_pattern_coefficient
logger = logging.getLogger(__name__.rsplit('.', maxsplit=1)[-1])


class _IgnoredException(Exception):
    """
    Base class for ignored exceptions.

    Not meant to be used directly, but only as a base class for
    other ignored exceptions.
    """
    def __init__(self, message, reason):
        # Call the base class constructor with the parameters it needs
        super().__init__(message)
        # Custom code...
        self.reason = reason


class SpectrumIgnored(_IgnoredException):
    """Spectrum ignored exception"""


class TraceIgnored(_IgnoredException):
    """Trace ignored exception"""


def _get_nint(config, trace):
    if config.trace_units == 'auto':
        instrtype = trace.stats.instrtype
    else:
        instrtype = config.trace_units
    if instrtype == 'acc':
        nint = 2
    elif instrtype in ['broadb', 'shortp', 'vel']:
        nint = 1
    elif instrtype == 'disp':
        nint = 0
    else:
        raise ValueError
    return nint


def _time_integrate(config, trace):
    nint = _get_nint(config, trace)
    trace.detrend(type='constant')
    trace.detrend(type='linear')
    for _ in range(nint):
        trace.data = cumtrapz(trace.data) * trace.stats.delta
        trace.stats.npts -= 1
        filter_trace(config, trace)


def _frequency_integrate(config, spec):
    nint = _get_nint(config, spec)
    for _ in range(nint):
        spec.data /= (2 * math.pi * spec.freq)


def _cut_spectrum(config, spec):
    # see if there is a station-specific frequency range
    station = spec.stats.station
    try:
        freq1 = float(config[f'freq1_{station}'])
        freq2 = float(config[f'freq2_{station}'])
    except KeyError:
        try:
            instrtype = spec.stats.instrtype
            freq1 = float(config[f'freq1_{instrtype}'])
            freq2 = float(config[f'freq2_{instrtype}'])
        except KeyError as e:
            raise RuntimeError(
                f'{spec.id}: Unknown instrument type: {instrtype}: '
                'skipping spectrum'
            ) from e
    return spec.slice(freq1, freq2)


def _compute_spec_h(spec_st, code, vertical_channel_codes=None, wave_type='S'):
    """
    Compute the component 'H' from geometric mean of the stream components.

    (which can also be all three components)
    """
    if vertical_channel_codes is None:
        vertical_channel_codes = ['Z']
    spec_h = None
    nspecs = 0
    spectral_snratio_mean = spectral_snratio_max = 0
    ignore_list = []
    ignore_reason_list = []
    for spec in spec_st:
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
        # do not use transverse component for SV
        if wave_type == 'SV' and channel[-1] == 'T':
            continue
        if spec_h is None:
            spec_h = spec.copy()
            spec_h.stats.channel = f'{code}H'
            # remove unneeded stats
            for attr in (
                'spectral_snratio',
                'ignore',
                'ignore_reason'
            ):
                spec_h.stats.pop(attr, None)
            spec_h.data = np.zeros_like(spec.data)
            spec_h.data_logspaced = np.zeros_like(spec.data_logspaced)
        ignore_list.append(spec.stats.ignore)
        ignore_reason_list.append(getattr(spec.stats, 'ignore_reason', None))
        if not spec.stats.ignore:
            spec_h.data += np.power(spec.data, 2)
            spec_h.data_logspaced += np.power(spec.data_logspaced, 2)
            spec_h.stats.spectral_snratio_fmin = min(
                getattr(spec.stats, 'spectral_snratio_fmin', np.nan),
                getattr(spec_h.stats, 'spectral_snratio_fmin', np.nan)
            )
            spec_h.stats.spectral_snratio_fmax = max(
                getattr(spec.stats, 'spectral_snratio_fmax', np.nan),
                getattr(spec_h.stats, 'spectral_snratio_fmax', np.nan)
            )
        spectral_snratio = getattr(
            spec.stats, 'spectral_snratio', 0)
        spectral_snratio_mean += spectral_snratio
        spectral_snratio_max = max(spectral_snratio_max, spectral_snratio)
        nspecs += 1
    # compute mean, avoid possible division by zero
    spectral_snratio_mean /= max(nspecs, 1)
    # finalize spec_h data and metadata
    spec_h.data = np.sqrt(spec_h.data)
    spec_h.data_logspaced = np.sqrt(spec_h.data_logspaced)
    spec_h.stats.spectral_snratio_mean = spectral_snratio_mean
    spec_h.stats.spectral_snratio_max = spectral_snratio_max
    spec_h.stats.ignore = all(ignore_list)
    if spec_h.stats.ignore:
        spec_h.data = np.ones_like(spec_h.data) * np.nan
        spec_h.data_logspaced = np.ones_like(spec_h.data_logspaced) * np.nan
        ignore_reason_list = {
            ir for ir in ignore_reason_list if ir is not None}
        spec_h.stats.ignore_reason = ', '.join(ignore_reason_list)
    return spec_h


def _check_data_len(config, trace):
    traceId = trace.get_id()
    trace_cut = trace.copy()
    if config.wave_type[0] == 'S':
        t1 = trace.stats.arrivals['S1'][1]
        t2 = trace.stats.arrivals['S2'][1]
    elif config.wave_type[0] == 'P':
        t1 = trace.stats.arrivals['P1'][1]
        t2 = trace.stats.arrivals['P2'][1]
    trace_cut.trim(starttime=t1, endtime=t2, pad=True, fill_value=0)
    npts = len(trace_cut.data)
    if npts == 0:
        msg = (
            f'{traceId}: No data for the selected cut interval: '
            'skipping trace'
        )
        raise TraceIgnored(msg, 'no data')
    nzeros = len(np.where(trace_cut.data == 0)[0])
    if nzeros > npts / 4:
        msg = (
            f'{traceId}: Signal window is zero for more than 25%: '
            'skipping trace'
        )
        raise TraceIgnored(msg, 'mostly zero signal window')


def _cut_signal_noise(config, trace):
    trace_signal = trace.copy()
    trace_noise = trace.copy()

    # Integrate in time domain, if required.
    # (otherwise frequency-domain integration is performed later)
    if config.time_domain_int:
        _time_integrate(config, trace_signal)
        _time_integrate(config, trace_noise)

    # trim...
    if config.wave_type[0] == 'S':
        t1 = trace.stats.arrivals['S1'][1]
        t2 = trace.stats.arrivals['S2'][1]
    elif config.wave_type[0] == 'P':
        t1 = trace.stats.arrivals['P1'][1]
        t2 = trace.stats.arrivals['P2'][1]
    trace_signal.trim(starttime=t1, endtime=t2, pad=True, fill_value=0)
    trace_signal.stats.type = 'signal'
    # Noise time window for weighting function:
    noise_t1 = trace.stats.arrivals['N1'][1]
    noise_t2 = trace.stats.arrivals['N2'][1]
    trace_noise.trim(
        starttime=noise_t1, endtime=noise_t2, pad=True, fill_value=0)
    trace_noise.stats.type = 'noise'
    # ...taper...
    cosine_taper(trace_signal.data, width=config.taper_halfwidth)
    cosine_taper(trace_noise.data, width=config.taper_halfwidth)

    # Be sure that both traces have same length:
    npts = min(len(trace_signal), len(trace_noise))
    if len(trace_signal) <= len(trace_noise):
        # truncate noise to signal length
        trace_noise.data = trace_noise.data[:npts]
    else:
        tr_noise_id = trace_noise.get_id()[:-1]
        signal_win_length = len(trace_signal) * trace_signal.stats.delta
        noise_win_length = len(trace_noise) * trace_noise.stats.delta
        if config.weighting == 'noise' and not config.force_noise_zero_padding:
            msg = (
                f'{tr_noise_id}: truncating signal window to noise length: '
                f'{signal_win_length:.1f} -> {noise_win_length:.1f} s'
            )
            trace_signal.data = trace_signal.data[:npts]
            # recompute signal window end time, so that it's plotted correctly
            _recompute_time_window(
                trace, config.wave_type[0], npts, keep='start')
        else:
            msg = (
                f'{tr_noise_id}: zero-padding noise window to signal length: '
                f'{noise_win_length:.1f} -> {signal_win_length:.1f} s'
            )
            # Notes:
            # 1. no risk of ringing, as noise has been tapered
            # 2. we use np.pad instead of obspy trim method
            #    to avoid potential rounding errors with start/end times
            # 3. we just pad at the end for simplicity (no changes
            #    in trace headers required) and it is identical to
            #    symmetric padding (tested)
            pad_len = len(trace_signal) - len(trace_noise)
            trace_noise.data = np.pad(trace_noise.data, (0, pad_len))
            # recompute noise window start time, so that it's plotted correctly
            _recompute_time_window(trace, 'N', len(trace_signal), keep='end')
        logger.warning(msg)
    signal_win_length = len(trace_signal) * trace_signal.stats.delta
    win_length_min = config.win_length_min or 0
    if signal_win_length < win_length_min:
        msg = (
            f'{trace.get_id()}: Signal window is too short: '
            f'{signal_win_length:.1f} s < {win_length_min:.1f}: skipping trace'
        )
        raise TraceIgnored(msg, 'window too short')
    return trace_signal, trace_noise


def _recompute_time_window(trace, wave_type, npts, keep='start'):
    """Recompute start or end time of signal or noise window,
    based on new number of points"""
    length = npts * trace.stats.delta
    if keep == 'end':
        label, _ = trace.stats.arrivals[f'{wave_type}1']
        t1 = trace.stats.arrivals[f'{wave_type}2'][1] - length
        trace.stats.arrivals[f'{wave_type}1'] = (label, t1)
    elif keep == 'start':
        label, _ = trace.stats.arrivals[f'{wave_type}2']
        t2 = trace.stats.arrivals[f'{wave_type}1'][1] + length
        trace.stats.arrivals[f'{wave_type}2'] = (label, t2)
    else:
        raise ValueError('keep must be "start" or "end"')


def _check_noise_level(trace_signal, trace_noise, config):
    traceId = trace_signal.get_id()
    trace_signal_rms = ((trace_signal.data**2).sum())**0.5
    # Scale trace_noise_rms to length of signal window,
    # based on length of non-zero noise window
    try:
        scale_factor = float(len(trace_signal)) / len(trace_noise.data != 0)
    except ZeroDivisionError:
        scale_factor = 1
    trace_noise_rms = ((trace_noise.data**2 * scale_factor).sum())**0.5
    if (
        trace_noise_rms / trace_signal_rms < 1e-6 and
        config.weighting == 'noise'
    ):
        # Skip trace if noise level is too low and if noise weighting is used
        msg = (
            f'{traceId}: Noise level is too low or zero: '
            'station will be skipped'
        )
        raise TraceIgnored(msg, 'low noise level')


def _geometrical_spreading_coefficient(config, spec):
    """
    Return the geometrical spreading coefficient for the given spectrum.
    """
    hypo_dist_in_km = spec.stats.hypo_dist
    epi_dist_in_km = spec.stats.epi_dist
    # set geometrical spreading distance to a very large value if it is None
    geom_spread_min_dist = config.geom_spread_min_teleseismic_distance or 1e99
    geom_spread_model =\
        'teleseismic' if epi_dist_in_km >= geom_spread_min_dist\
        else config.geom_spread_model
    logger.info(f'{spec.id}: geometrical spreading model: {geom_spread_model}')
    if geom_spread_model == 'teleseismic':
        angular_distance = spec.stats.gcarc
        source_depth_in_km = spec.stats.event.hypocenter.depth.value_in_km
        station_depth_in_km = -spec.stats.coords.elevation
        phase = config.wave_type[0]
        return geom_spread_teleseismic(
            angular_distance, source_depth_in_km, station_depth_in_km, phase)
    if config.geom_spread_model == 'r_power_n':
        exponent = config.geom_spread_n_exponent
        return geom_spread_r_power_n(hypo_dist_in_km, exponent)
    if config.geom_spread_model == 'r_power_n_segmented':
        exponents = config.geom_spread_n_exponents
        hinge_distances = config.geom_spread_n_distances
        if len(hinge_distances) != len(exponents):
            raise ValueError(
                f'The number of exponents must be equal to the number of '
                f'hinge distances. You provided {len(exponents)} exponents '
                f'and {len(hinge_distances)} hinge distances'
            )
        return geom_spread_r_power_n_segmented(hypo_dist_in_km, exponents,
                                               hinge_distances)
    if config.geom_spread_model == 'boatwright':
        cutoff_dist_in_km = config.geom_spread_cutoff_distance
        return geom_spread_boatwright(
            hypo_dist_in_km, cutoff_dist_in_km, spec.freq)
    raise ValueError(
        f'Unknown geometrical spreading model: {config.geom_spread_model}')


def _get_free_surface_amplification(stats, config):
    """
    Return the free-surface amplification factor for the given station.
    """
    traceid = '.'.join((
        stats.network, stats.station, stats.location, stats.channel))
    fsa = config.free_surface_amplification
    # If scalar, just use it for all
    if (
        isinstance(fsa, (np.generic, numbers.Real))
        and not isinstance(fsa, bool)
    ):
        return float(fsa)
    # If tuple of tuples, convert to list for matching
    if isinstance(fsa, tuple):
        fsa_items = list(fsa)
    else:
        raise ValueError(
            'free_surface_amplification must be numeric or tuple of tuples')

    # Sort patterns by specificity (most specific first):
    # - More literal characters -> more specific
    # - Fewer '*' (greedy wildcards) -> more specific
    # - Fewer '?' (single-char wildcards) -> more specific
    # - Longer patterns -> more specific (tie-breaker)
    def _get_specificity(pattern):
        # Ensure '*' alone is always last
        if pattern == '*':
            return (1, 0, 0, 0, 0)
        num_literals = sum(c not in ('*', '?') for c in pattern)
        num_star = pattern.count('*')
        num_qmark = pattern.count('?')
        return (0, -num_literals, num_star, num_qmark, -len(pattern))

    sorted_items = sorted(
        fsa_items, key=lambda item: _get_specificity(item[0]))
    # Find the first matching pattern
    for pattern, value in sorted_items:
        if fnmatch.fnmatch(traceid, pattern):
            return float(value)
    raise ValueError(
        f'No free-surface amplification factor found for station {traceid}'
    )


# store log messages to avoid duplicates
PROPERTY_LOG_MESSAGES = []


def _displacement_to_moment(stats, config):
    """
    Return the coefficient for converting displacement to seismic moment.

    From Aki&Richards,1980
    """
    phase = config.wave_type[0]
    lon = stats.coords.longitude
    lat = stats.coords.latitude
    depth = -stats.coords.elevation
    medium_properties = MediumProperties(lon, lat, depth, config)
    depth_string = medium_properties.to_string('station depth', depth)
    v_name = f'v{phase.lower()}'
    v_source = config.event.hypocenter[v_name]
    v_source_string = medium_properties.to_string(f'{v_name}_source', v_source)
    v_station = medium_properties.get(mproperty=v_name, where='stations')
    stats.v_station = v_station
    stats.v_station_type = phase
    v_station_string = medium_properties.to_string(
        f'{v_name}_station', v_station)
    rho_source = config.event.hypocenter.rho
    rho_source_string = medium_properties.to_string('rho_source', rho_source)
    rho_station = medium_properties.get(mproperty='rho', where='stations')
    stats.rho_station = rho_station
    rho_station_string = medium_properties.to_string(
        'rho_station', rho_station)
    specid = '.'.join((
        stats.network, stats.station, stats.location, stats.channel))
    msg = (
        f'{specid}: {depth_string}, '
        f'{v_source_string}, {v_station_string}'
    )
    if msg not in PROPERTY_LOG_MESSAGES:
        logger.info(msg)
        PROPERTY_LOG_MESSAGES.append(msg)
    msg = (
        f'{specid}: {depth_string}, '
        f'{rho_source_string}, {rho_station_string}'
    )
    if msg not in PROPERTY_LOG_MESSAGES:
        logger.info(msg)
        PROPERTY_LOG_MESSAGES.append(msg)
    v_source *= 1000.
    v_station *= 1000.
    v3 = v_source**(5. / 2) * v_station**(1. / 2)
    rho = rho_source**0.5 * rho_station**0.5
    fsa = _get_free_surface_amplification(stats, config)
    stats.fsa = fsa
    msg = f'{specid}: free-surface amplification factor: {fsa:.2f}'
    if msg not in PROPERTY_LOG_MESSAGES:
        logger.info(msg)
        PROPERTY_LOG_MESSAGES.append(msg)
    return 4 * math.pi * v3 * rho / (fsa * stats.radiation_pattern)


def _pre_smoothing_snr_mask(config, spec, specnoise):
    """
    Suppress low-frequency bins dominated by noise *before* log-frequency smoothing.
    This prevents the smoothing kernel from bleeding high noise amplitudes
    into the usable band.
    Policy:
      - Compute per-bin SNR on linear spectra (moment units) after
        geometrical spreading & conversion (same stage as in _build_spectrum,
        just before smoothing).
      - If 'spectral_snr_mask_threshold' <= 0 or missing -> no-op.
      - mode 'left_of_first' (default): find the first frequency f_cross where
        SNR >= threshold and clamp all lower-frequency amplitudes to at most the
        amplitude at f_cross. Optionally apply a gentle ramp (half-cosine) over
        'spectral_snr_mask_ramp_decades' to avoid a hard step.
      - mode 'binary': clamp every bin with SNR < threshold to at most the
        minimum unmasked amplitude (conservative).
    Notes:
      - We *do not* insert NaNs: the smoother needs finite values (ssp_util.smooth).
      - This runs before _smooth_spectrum(), i.e., before 'make_freq_logspaced()'.
    """
    th = getattr(config, 'spectral_snr_mask_threshold', None)
    if not th or th <= 0:
        return
    mode = getattr(config, 'spectral_snr_mask_mode', 'left_of_first')
    ramp_dec = float(getattr(config, 'spectral_snr_mask_ramp_decades', 0.0) or 0.0)

    # spec and specnoise share the same linear frequency grid at this stage
    freq = spec.freq
    s = spec.data
    n = specnoise.data
    if s.size == 0 or n.size == 0 or freq.size == 0:
        return
    eps = np.finfo(float).eps
    sn = s / np.maximum(n, eps)

    mask_info = {
        'threshold': float(th),
        'mode': mode,
        'ramp_decades': float(ramp_dec),
        'applied': False
    }

    freq_range_cfg = getattr(config, 'spectral_snr_mask_freq_range', None)
    freq_range = None
    range_idx = None
    if freq_range_cfg:
        try:
            if isinstance(freq_range_cfg, (list, tuple)):
                values = list(freq_range_cfg)
            else:
                cleaned = str(freq_range_cfg).replace('[', '').replace(']', '')
                values = [
                    float(val.strip()) for val in cleaned.split(',') if val.strip()
                ]
            if not values:
                raise ValueError
            if len(values) == 1:
                f_range_min = float(values[0])
                f_range_max = 1.0
            else:
                f_range_min = float(values[0])
                f_range_max = float(values[1])
        except (TypeError, ValueError):
            logger.warning(
                'Invalid spectral SNR mask frequency range %r: ignoring',
                freq_range_cfg)
        else:
            if f_range_max < f_range_min:
                f_range_min, f_range_max = f_range_max, f_range_min
            f_range_min = max(f_range_min, freq[0])
            f_range_max = min(max(f_range_min, f_range_max), 1.0)
            if f_range_min < f_range_max:
                # reduce to available frequency support
                idx = np.where((freq >= f_range_min) & (freq <= f_range_max))[0]
                if idx.size:
                    range_idx = idx
                    # store effective range for logging/plots
                    freq_range = (float(freq[idx[0]]), float(freq[idx[-1]]))
                    mask_info.update({
                        'freq_range_min': freq_range[0],
                        'freq_range_max': freq_range[1]
                    })
                else:
                    logger.info(
                        '%s: spectral SNR mask frequency range [%g, %g] Hz '
                        'does not intersect spectrum: ignoring',
                        spec.get_id(), f_range_min, f_range_max)
            else:
                logger.warning(
                    'Invalid spectral SNR mask frequency range %r: ignoring',
                    freq_range_cfg)

    if mode == 'left_of_first':
        range_forced = range_idx is not None and np.any(sn[range_idx] < th)
        if range_forced:
            idx = range_idx[-1]
        else:
            # Index of first crossing where SNR >= th; if none, do nothing.
            idx = np.argmax(sn >= th)
            if sn[idx] < th:
                # never crosses threshold
                spec.stats.snr_mask_info = mask_info
                spec_id = spec.get_id()
                logger.info(
                    f'{spec_id}: pre-smoothing SNR mask enabled '
                    f'(threshold={th:g}), but no crossing found: skipped')
                return
        A0 = s[idx]
        f0 = freq[idx]
        # Clamp everything left of f0 to at most A0
        left = np.where(freq < f0)[0]
        if left.size:
            s[left] = np.minimum(s[left], A0)
        # Optional soft ramp (half-cosine) across a band of width 'ramp_dec' decades
        if ramp_dec > 0.0:
            f1 = f0 / (10.0 ** ramp_dec)
            if freq_range is not None:
                f1 = max(f1, freq_range[0])
            band = np.where((freq >= f1) & (freq < f0))[0]
            if band.size:
                # w goes 0->1 across [f1, f0]; blend clamped (A0) with original
                x = (freq[band] - f1) / max(f0 - f1, eps)
                w = 0.5 * (1.0 - np.cos(np.pi * x))
                s[band] = np.minimum(s[band], A0 * w + s[band] * (1.0 - w))
        spec.data = s
        mask_info.update({
            'applied': True,
            'f_cross': float(f0),
        })
        if range_forced:
            mask_info['range_forced'] = True
        if ramp_dec > 0.0:
            mask_info['f_ramp_start'] = float(f1)
        spec.stats.snr_mask_info = mask_info
        spec_id = spec.get_id()
        logger.info(f'{spec_id}: pre-smoothing SNR mask applied '
                    f'(threshold={th:g}, mode={mode}, f_cross={f0:.4f} Hz, '
                    f'ramp_decades={ramp_dec:.3f})')
        if range_forced and freq_range is not None:
            logger.info(
                '%s: pre-smoothing SNR mask forced by RSS range %.4f-%.4f Hz',
                spec_id, freq_range[0], freq_range[1])
    else:
        # 'binary' fallback: clamp any low-SNR bin.
        mask = sn < th
        if np.any(mask):
            # Use the minimum amplitude among unmasked bins as a conservative cap.
            unmasked = np.where(~mask)[0]
            A0 = np.min(s[unmasked]) if unmasked.size else np.min(s)
            s[mask] = np.minimum(s[mask], A0)
            spec.data = s
            mask_info.update({
                'applied': True,
                'freq_min': float(np.min(freq[mask])),
                'freq_max': float(np.max(freq[mask])),
                'masked_bins': int(np.count_nonzero(mask)),
            })
        spec.stats.snr_mask_info = mask_info


def _apply_mask_info_to_component(spec, mask_info):
    """Apply RSS-derived mask information to a single component spectrum."""
    if not mask_info:
        return
    # Store a copy of the mask info so later consumers (plots) can reuse it.
    spec.stats.snr_mask_info = dict(mask_info)
    if not mask_info.get('applied'):
        return
    mode = mask_info.get('mode', 'left_of_first')
    freq = spec.freq
    data = spec.data
    if mode == 'left_of_first':
        f_cross = mask_info.get('f_cross')
        if not f_cross:
            return
        ramp_dec = float(mask_info.get('ramp_decades', 0.0) or 0.0)
        eps = np.finfo(float).eps
        idx_candidates = np.where(freq >= f_cross)[0]
        if not idx_candidates.size:
            return
        idx = idx_candidates[0]
        A0 = data[idx]
        left = np.where(freq < f_cross)[0]
        if left.size:
            data[left] = np.minimum(data[left], A0)
        if ramp_dec > 0.0:
            f_ramp_start = mask_info.get('f_ramp_start')
            if f_ramp_start is None:
                f_ramp_start = f_cross / (10.0 ** ramp_dec)
            if f_ramp_start < f_cross:
                band = np.where((freq >= f_ramp_start) & (freq < f_cross))[0]
                if band.size:
                    x = (freq[band] - f_ramp_start) / max(f_cross - f_ramp_start, eps)
                    w = 0.5 * (1.0 - np.cos(np.pi * x))
                    data[band] = np.minimum(data[band], A0 * w + data[band] * (1.0 - w))
        spec.data = data
    else:
        freq_min = mask_info.get('freq_min')
        freq_max = mask_info.get('freq_max')
        if freq_min is None or freq_max is None:
            return
        mask = (freq >= freq_min) & (freq <= freq_max)
        if np.any(mask):
            unmasked = np.where(~mask)[0]
            A0 = np.min(data[unmasked]) if unmasked.size else np.min(data)
            data[mask] = np.minimum(data[mask], A0)
            spec.data = data


def _apply_rss_snr_mask(config, spec_st, specnoise_st):
    """Apply the spectral SNR mask using the root sum of squares spectrum."""
    if not spec_st:
        return
    spec_ids = {sp.id[:-1] for sp in spec_st}
    vertical_codes = config.vertical_channel_codes
    wave_type = config.wave_type
    for specid in spec_ids:
        spec_sel = _select_spectra(spec_st, specid)
        if specnoise_st:
            specnoise_sel = _select_spectra(specnoise_st, specid)
        else:
            specnoise_sel = SpectrumStream()
        codes = {sp.stats.channel[:-1] for sp in spec_sel}
        for code in codes:
            components = SpectrumStream(
                sp for sp in spec_sel
                if sp.stats.channel[:-1] == code and
                sp.stats.channel[-1] not in ['H', 'h', 'S', 's', 't']
            )
            if not components:
                continue
            noises = SpectrumStream(
                sp for sp in specnoise_sel
                if sp.stats.channel[:-1] == code and
                sp.stats.channel[-1] not in ['H', 'h', 'S', 's', 't']
            )
            spec_h = _compute_spec_h(
                components, code, vertical_codes, wave_type)
            if spec_h is None or getattr(spec_h.stats, 'ignore', False):
                continue
            specnoise_h = _compute_spec_h(
                noises, code, vertical_codes, wave_type)
            if specnoise_h is None:
                continue
            _pre_smoothing_snr_mask(config, spec_h, specnoise_h)
            mask_info = getattr(spec_h.stats, 'snr_mask_info', None)
            if not mask_info:
                continue
            mask_info = dict(mask_info)
            mask_info['source_channel'] = spec_h.stats.channel
            for component in components:
                if getattr(component.stats, 'ignore', False):
                    continue
                _apply_mask_info_to_component(component, mask_info)


def _smooth_spectrum(spec, smooth_width_decades=0.2):
    """
    Smooth spectrum in a log10-freq space.

    This function also adds to the original spectrum the log-spaced
    frequencies and data, which are used in the inversion.
    """
    # 1. Generate log10-spaced frequencies
    spec.make_freq_logspaced()
    # 2. Interpolate data to log10-spaced frequencies
    spec.make_logspaced_from_linear()
    # 3. Smooth log10-spaced data points
    log_df = spec.stats.delta_logspaced
    npts = max(1, int(round(smooth_width_decades / log_df)))
    spec.data_logspaced = smooth(spec.data_logspaced, window_len=npts)
    # 4. Reinterpolate to linear frequencies
    spec.make_linear_from_logspaced()
    # 5. Optimize the sampling rate of log spectrum,
    #    based on the width of the smoothing window
    log_df = smooth_width_decades / 5
    # note: make_freq_logspaced() will reinterpolate the log-spaced data
    # to the new log-spaced frequencies
    spec.make_freq_logspaced(log_df)


def _build_spectrum(config, trace, do_smooth=True):
    try:
        spec = Spectrum(obspy_trace=trace)
    except ValueError as e:
        raise RuntimeError(
            f'{trace.id}: Error building spectrum: skipping spectrum\n{str(e)}'
        ) from e
    spec.stats.instrtype = trace.stats.instrtype
    spec.stats.coords = trace.stats.coords
    spec.stats.event = trace.stats.event
    spec.stats.hypo_dist = trace.stats.hypo_dist
    spec.stats.epi_dist = trace.stats.epi_dist
    spec.stats.gcarc = trace.stats.gcarc
    spec.stats.azimuth = trace.stats.azimuth
    spec.stats.travel_times = trace.stats.travel_times
    spec.stats.takeoff_angles = trace.stats.takeoff_angles
    spec.stats.ignore = trace.stats.ignore
    # Integrate in frequency domain, if no time-domain
    # integration has been performed
    if not config.time_domain_int:
        _frequency_integrate(config, spec)
    # cut the spectrum
    spec = _cut_spectrum(config, spec)
    # correct geometrical spreading
    try:
        geom_spread = _geometrical_spreading_coefficient(config, spec)
    except Exception as e:
        raise RuntimeError(
            f'{spec.id}: Error computing geometrical spreading: '
            f'skipping spectrum\n{str(e)}'
        ) from e
    spec.data *= geom_spread
    # store the radiation pattern coefficient in the spectrum stats
    spec.stats.radiation_pattern =\
        get_radiation_pattern_coefficient(spec.stats, config)
    # convert to seismic moment
    try:
        coeff = _displacement_to_moment(spec.stats, config)
        spec.data *= coeff
        # store coeff to correct back data in displacement units
        # for radiated_energy()
        spec.stats.coeff = coeff
    except ValueError as e:
        raise RuntimeError(
            f'{spec.id}: Error converting to seismic moment: '
            f'skipping spectrum\n{str(e)}'
        ) from e
    # smooth spectrum. This also creates log-spaced frequencies and data
    if do_smooth:
        try:
            _smooth_spectrum(spec, config.spectral_smooth_width_decades)
        except ValueError as e:
            raise RuntimeError(
                f'{spec.id}: Error smoothing spectrum: '
                f'skipping spectrum\n{str(e)}'
            ) from e
    return spec


def _build_uniform_weight(spec):
    weight = spec.copy()
    # Clear data_mag and data_mag_logspaced,
    # since magnitude values are not relevant for weights
    weight.data_mag = []
    weight.data_mag_logspaced = []
    weight.snratio = None
    weight.data = np.ones_like(weight.data)
    weight.data_logspaced = np.ones_like(weight.data_logspaced)
    return weight


def _build_weight_from_frequency(config, spec):
    weight = spec.copy()
    # Clear data_mag and data_mag_logspaced,
    # since magnitude values are not relevant for weights
    weight.data_mag = []
    weight.data_mag_logspaced = []
    freq = weight.freq
    weight.data = np.ones_like(weight.data)
    weight.data[freq <= config.f_weight] = config.weight
    weight.data /= np.nanmax(weight.data)
    freq_logspaced = weight.freq_logspaced
    weight.data_logspaced = np.ones_like(weight.data_logspaced)
    weight.data_logspaced[freq_logspaced <= config.f_weight] = config.weight
    weight.data_logspaced /= np.nanmax(weight.data_logspaced)
    return weight


def _build_weight_from_inv_frequency(spec, power=0.25):
    """
    Build spectral weights from inverse frequency (raised to a power < 1)
    """
    if power >= 1:
        raise ValueError('pow must be < 1')
    # Note: weight.data is used for plotting,
    #       weight.data_logspaced for actual weighting
    weight = spec.copy()
    # Clear data_mag and data_mag_logspaced,
    # since magnitude values are not relevant for weights
    weight.data_mag = []
    weight.data_mag_logspaced = []
    freq = weight.freq
    weight.data *= 0
    # Limit non-zero weights to fmin/fmax from spectral_snratio if available
    snr_fmin = getattr(spec.stats, 'spectral_snratio_fmin', None)
    snr_fmax = getattr(spec.stats, 'spectral_snratio_fmax', None)
    # Check if snr_fmin and snr_fmax are valid (not NaN)
    if snr_fmin is not None and not np.isnan(snr_fmin):
        i0 = np.where(freq >= snr_fmin)[0][0]
    else:
        i0 = 0
    if snr_fmax is not None and not np.isnan(snr_fmax):
        i1 = np.where(freq <= snr_fmax)[0][-1]
    else:
        i1 = len(freq) - 1
    # Build weights as if frequencies always start from 0.25 Hz
    # to obtain similar curves regardless of fmin
    # and to avoid too much weight for very low frequencies
    weight.data[i0: i1 + 1] = 1. / (freq[i0: i1 + 1] - freq[i0] + 0.25)**power
    weight.data /= np.nanmax(weight.data)
    freq_logspaced = weight.freq_logspaced
    weight.data_logspaced *= 0
    i0 = np.where(freq_logspaced >= snr_fmin)[0][0] if snr_fmin else 0
    i1 = np.where(freq_logspaced <= snr_fmax)[0][-1] if snr_fmax\
        else len(freq_logspaced) - 1
    weight.data_logspaced[i0: i1 + 1] =\
        1. / (freq_logspaced[i0: i1 + 1] - freq_logspaced[i0] + 0.25)**power
    weight.data_logspaced /= np.nanmax(weight.data_logspaced)
    return weight


def _build_weight_from_ratio(spec, specnoise, smooth_width_decades):
    weight = spec.copy()
    # Clear data_mag and data_mag_logspaced,
    # since magnitude values are not relevant for weights
    weight.data_mag = []
    weight.data_mag_logspaced = []
    try:
        weight.data /= specnoise.data
    except ValueError as e:
        logger.error(
            f'{spec.id}: Error building weight from noise: {str(e)}'
        )
        ssp_exit(1)
    # save signal-to-noise ratio before log10, smoothing, and normalization
    weight.snratio = weight.data.copy()
    # The inversion is done in magnitude units,
    # so let's take log10 of weight
    weight.data = np.log10(weight.data)
    # a small value used to replace negative and NaN values
    _small_value = 1e-9
    weight.data[weight.data < 0] = _small_value
    # replace NaN values with a small value
    weight.data[np.isnan(weight.data)] = _small_value
    # Weight spectrum is smoothed once more
    _smooth_spectrum(weight, smooth_width_decades)
    weight.data /= np.nanmax(weight.data)
    # slightly taper weight at low frequencies, to avoid overestimating
    # weight at low frequencies, in cases where noise is underestimated
    cosine_taper(
        weight.data,
        min(0.25, weight.stats.delta / 4),
        left_taper=True)
    # Replace NaN values with a small value (in case the taper introduces NaNs)
    weight.data[np.isnan(weight.data)] = _small_value
    # Replace weights <= 0.2 with a small value
    weight.data[weight.data <= 0.2] = _small_value
    return weight


def _build_weight_from_noise(config, spec, specnoise):
    if specnoise is None or np.all(specnoise.data == 0):
        spec_id = spec.get_id()[:-1]
        logger.warning(
            f'{spec_id}: No available noise window: '
            'a uniform weight will be applied')
        weight = _build_uniform_weight(spec)
    else:
        weight = _build_weight_from_ratio(
            spec, specnoise, config.spectral_smooth_width_decades)
    # interpolate to log-frequencies
    weight.make_logspaced_from_linear()
    weight.data_logspaced /= np.nanmax(weight.data_logspaced)
    # Make sure weight is positive
    weight.data_logspaced[weight.data_logspaced <= 0] = 0.001
    return weight


def _build_weight_spectral_stream(config, spec_st, specnoise_st):
    """Build a stream of weights from a stream of spectra and a stream of
    noise spectra."""
    weight_st = SpectrumStream()
    spec_ids = {sp.id[:-1] for sp in spec_st}
    for specid in spec_ids:
        try:
            spec_h = _select_spectra(spec_st, f'{specid}H')[0]
            specnoise_h = _select_spectra(specnoise_st, f'{specid}H')[0]
        except Exception:
            continue
        if config.weighting == 'noise':
            weight = _build_weight_from_noise(config, spec_h, specnoise_h)
        elif config.weighting == 'frequency':
            weight = _build_weight_from_frequency(config, spec_h)
        elif config.weighting == 'inv_frequency':
            weight = _build_weight_from_inv_frequency(spec_h)
        elif config.weighting == 'no_weight':
            weight = _build_uniform_weight(spec_h)
        weight_st.append(weight)
    return weight_st


def _select_spectra(spec_st, specid):
    """Select spectra from stream, based on specid."""
    network, station, location, channel = specid.split('.')
    channel = channel + '?' * (3 - len(channel))
    spec_st_sel = spec_st.select(
        network=network, station=station, location=location, channel=channel)
    return SpectrumStream(spec_st_sel)


def _build_H(spec_st, specnoise_st=None, vertical_channel_codes=None,
             wave_type='S'):
    """
    Add to spec_st and specnoise_st the "H" component.

    H component is obtained from the modulus of all the available components.
    """
    if vertical_channel_codes is None:
        vertical_channel_codes = ['Z']
    spec_ids = {sp.id[:-1] for sp in spec_st}
    for specid in spec_ids:
        spec_st_sel = _select_spectra(spec_st, specid)
        specnoise_st_sel = _select_spectra(specnoise_st, specid)
        # 'code' is band+instrument code
        for code in {x.stats.channel[:-1] for x in spec_st_sel}:
            spec_h = _compute_spec_h(
                spec_st_sel, code, vertical_channel_codes, wave_type)
            if spec_h is not None:
                spec_st.append(spec_h)
            specnoise_h = _compute_spec_h(
                specnoise_st_sel, code, vertical_channel_codes, wave_type)
            if specnoise_h is not None:
                specnoise_st.append(specnoise_h)


def _check_spectral_sn_ratio(config, spec, specnoise):
    spec_id = spec.get_id()
    weight = _build_weight_from_noise(config, spec, specnoise)
    freqs = weight.freq
    # if no noise window is available, snratio is not computed
    if weight.snratio is None:
        spec.stats.spectral_snratio = None
        return
    if config.spectral_sn_freq_range is not None:
        sn_fmin, sn_fmax = config.spectral_sn_freq_range[:2]
        valid_freqs_idx = np.where((sn_fmin <= freqs) * (freqs <= sn_fmax))
        valid_snratio = weight.snratio[valid_freqs_idx]
    else:
        valid_snratio = weight.snratio
    if len(valid_snratio) == 0:
        spec.stats.spectral_snratio = np.nan
        msg = (
            f'{spec_id}: no valid frequency to compute spectral S/N: '
            'ignoring spectrum'
        )
        reason = 'no valid frequency'
        raise SpectrumIgnored(msg, reason)
    spectral_snratio = valid_snratio.mean()
    spec.stats.spectral_snratio = spectral_snratio
    # Save frequency range where SNR > 3 so it can be used for building weights
    # Note: not sure if we could use config.spectral_sn_min here instead of 3
    snr_valid_freqs = freqs[weight.snratio >= 3]
    try:
        spec.stats.spectral_snratio_fmin = snr_valid_freqs[0]
        spec.stats.spectral_snratio_fmax = snr_valid_freqs[-1]
    except IndexError:
        spec.stats.spectral_snratio_fmin = np.nan
        spec.stats.spectral_snratio_fmax = np.nan
    logger.info(f'{spec_id}: average spectral S/N: {spectral_snratio:.2f}')
    ssnmin = config.spectral_sn_min or -np.inf
    if spectral_snratio < ssnmin:
        msg = (
            f'{spec_id}: spectral S/N smaller than {ssnmin:.2f}: '
            'ignoring spectrum')
        reason = 'low spectral S/N'
        raise SpectrumIgnored(msg, reason)


def _ignore_spectrum(msg, spec, specnoise):
    """Ignore spectrum. Set ignore flag and reason."""
    logger.warning(msg)
    spec.stats.ignore = True
    spec.stats.ignore_reason = msg.reason
    specnoise.stats.ignore = True
    specnoise.stats.ignore_reason = msg.reason


def _ignore_trace(msg, trace):
    """Ignore trace. Set ignore flag and reason."""
    # NOTE: no logger.warning here, because it is already done in
    # _ignore_spectrum()
    trace.stats.ignore = True
    trace.stats.ignore_reason = msg.reason


def _build_signal_and_noise_streams(config, st):
    """Build signal and noise streams."""
    # remove traces with ignore flag
    traces = Stream([tr for tr in st if not tr.stats.ignore])
    # Compute and store in config the number of input stations
    # (i.e., trace id without component code). This value will be stored
    # later in the sspec_output object.
    n_input_stations = len({tr.id[:-1] for tr in traces})
    config.n_input_stations = n_input_stations
    signal_st = Stream()
    noise_st = Stream()
    for trace in sorted(traces, key=lambda tr: tr.id):
        try:
            _check_data_len(config, trace)
            trace_signal, trace_noise = _cut_signal_noise(config, trace)
            _check_noise_level(trace_signal, trace_noise, config)
            signal_st.append(trace_signal)
            noise_st.append(trace_noise)
        except RuntimeError as msg:
            # this should not happen, since the function above
            # raise a TraceIgnored exception
            logger.warning(msg)
            continue
        except TraceIgnored as msg:
            _ignore_trace(msg, trace)
            logger.warning(msg)
            continue
    return signal_st, noise_st


def _trim_components(config, signal_st, noise_st, st):
    """
    Trim components of the same instrument to the same number of samples.

    Recompute time window of the signal and noise traces for correct plotting.
    """
    for traceid in sorted({tr.id[:-1] for tr in signal_st}):
        st_sel = signal_st.select(id=f'{traceid}*') +\
            noise_st.select(id=f'{traceid}*')
        all_npts = {tr.stats.npts for tr in st_sel}
        if len(all_npts) == 1:
            continue
        logger.warning(
            f'{traceid}: components have different window lengths. '
            'Trimming signal and noise windows to the shortest one')
        npts = min(all_npts)
        for tr in st_sel:
            if tr.stats.type == 'signal':
                tr.data = tr.data[:npts]
            elif tr.stats.type == 'noise':
                tr.data = tr.data[-npts:]
        for tr in st.select(id=f'{traceid}*'):
            _recompute_time_window(tr, config.wave_type[0], npts, keep='start')
            _recompute_time_window(tr, 'N', npts, keep='end')


def _build_signal_and_noise_spectral_streams(
        config, signal_st, noise_st, original_st):
    """
    Build signal and noise spectral streams.

    Note: original_st is only used to keep track of ignored traces.
    """
    spec_st = SpectrumStream()
    specnoise_st = SpectrumStream()
    spec_items = []

    th = getattr(config, 'spectral_snr_mask_threshold', None)
    if isinstance(th, str):
        try:
            th = float(th)
        except ValueError:
            logger.warning(
                'Invalid spectral SNR mask threshold %r: disabling mask',
                th)
            th = None
        else:
            setattr(config, 'spectral_snr_mask_threshold', th)
    use_mask = th is not None and th > 0

    for trace_signal in sorted(signal_st, key=lambda tr: tr.id):
        trace_noise = noise_st.select(id=trace_signal.id)[0]
        try:
            spec = _build_spectrum(
                config, trace_signal, do_smooth=not use_mask)
            specnoise = _build_spectrum(
                config, trace_noise, do_smooth=not use_mask)
        except RuntimeError as msg:
            logger.warning(msg)
            continue
        spec_items.append((spec, specnoise, trace_signal))
        spec_st.append(spec)
        specnoise_st.append(specnoise)

    if not spec_st:
        logger.error('No spectra left! Exiting.')
        ssp_exit()

    if use_mask:
        # Apply the mask on root-sum-of-squares spectra before smoothing.
        _apply_rss_snr_mask(config, spec_st, specnoise_st)
        for spec in spec_st:
            _smooth_spectrum(spec, config.spectral_smooth_width_decades)
        for specnoise in specnoise_st:
            _smooth_spectrum(specnoise, config.spectral_smooth_width_decades)

    for spec, specnoise, trace_signal in spec_items:
        try:
            _check_spectral_sn_ratio(config, spec, specnoise)
        except SpectrumIgnored as msg:
            _ignore_spectrum(msg, spec, specnoise)
            trace_original = original_st.select(id=trace_signal.id)[0]
            _ignore_trace(msg, trace_original)
        except RuntimeError as msg:
            logger.warning(msg)

    # build H component
    _build_H(
        spec_st, specnoise_st, config.vertical_channel_codes, config.wave_type)
    # convert the spectral amplitudes to moment magnitude
    for spec in spec_st:
        spec.data_mag = moment_to_mag(spec.data)
        spec.data_mag_logspaced = moment_to_mag(spec.data_logspaced)
    for specnoise in specnoise_st:
        specnoise.data_mag = moment_to_mag(specnoise.data)
    # apply station correction if a residual file is specified in config
    station_correction(spec_st, config)
    return spec_st, specnoise_st


def _zero_pad(config, trace):
    """Zero-pad trace to spectral_win_length"""
    spec_win_len = config.spectral_win_length
    if spec_win_len is None:
        return
    wtype = config.wave_type[0]
    if trace.stats.type == 'signal':
        t1 = trace.stats.arrivals[f'{wtype}1'][1]
    elif trace.stats.type == 'noise':
        t1 = trace.stats.arrivals['N1'][1]
    trace.trim(starttime=t1, endtime=t1 + spec_win_len, pad=True, fill_value=0)


def build_spectra(config, st):
    """
    Build spectra and the ``spec_st`` object.

    Computes P- or S-wave (displacement) spectra from
    accelerometers and velocimeters, uncorrected for anelastic attenuation,
    corrected for instrumental constants, normalized by geometrical spreading.
    """
    wave_type = config.wave_type
    logger.info(f'Building {wave_type}-wave spectra...')
    signal_st, noise_st = _build_signal_and_noise_streams(config, st)
    _trim_components(config, signal_st, noise_st, st)
    for trace in signal_st + noise_st:
        _zero_pad(config, trace)
    spec_st, specnoise_st = _build_signal_and_noise_spectral_streams(
        config, signal_st, noise_st, st)
    weight_st = _build_weight_spectral_stream(config, spec_st, specnoise_st)
    logger.info(f'Building {wave_type}-wave spectra: done')
    logger.info('---------------------------------------------------')
    return spec_st, specnoise_st, weight_st
