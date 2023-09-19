# -*- coding: utf-8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Trace processing for sourcespec.

:copyright:
    2012 Claudio Satriano <satriano@ipgp.fr>

    2013-2014 Claudio Satriano <satriano@ipgp.fr>,
              Emanuela Matrullo <matrullo@geologie.ens.fr>

    2015-2023 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
import logging
import re
import numpy as np
from scipy.signal import savgol_filter
from obspy.core import Stream
from obspy.core.util import AttribDict
from sourcespec.ssp_setup import ssp_exit
from sourcespec.ssp_util import (
    remove_instr_response, station_to_event_position)
from sourcespec.ssp_wave_arrival import add_arrival_to_trace
from sourcespec.clipping_detection import (
    compute_clipping_score, clipping_peaks)
logger = logging.getLogger(__name__.rsplit('.', maxsplit=1)[-1])


def _get_bandpass_frequencies(config, trace):
    """Get frequencies for bandpass filter."""
    # see if there is a station-specific filter
    station = trace.stats.station
    try:
        bp_freqmin = float(config[f'bp_freqmin_{station}'])
        bp_freqmax = float(config[f'bp_freqmax_{station}'])
    except KeyError:
        instrtype = trace.stats.instrtype
        try:
            bp_freqmin = float(config[f'bp_freqmin_{instrtype}'])
            bp_freqmax = float(config[f'bp_freqmax_{instrtype}'])
        except KeyError as e:
            raise ValueError(
                f'{trace.id}: Unknown instrument type: {instrtype}: '
                'skipping trace'
            ) from e
    return bp_freqmin, bp_freqmax


def filter_trace(config, trace):
    """Filter trace."""
    bp_freqmin, bp_freqmax = _get_bandpass_frequencies(config, trace)
    nyquist = 1. / (2. * trace.stats.delta)
    if bp_freqmax >= nyquist:
        bp_freqmax = nyquist * 0.999
        logger.warning(
            f'{trace.id}: maximum frequency for bandpass filtering '
            f'is larger or equal to Nyquist. Setting it to {bp_freqmax} Hz')
    filter_options = {
        'type': 'bandpass',
        'freqmin': bp_freqmin,
        'freqmax': bp_freqmax
    }
    trace.filter(**filter_options)
    # save filter info to trace stats
    trace.stats.filter = AttribDict(filter_options)


def _check_signal_level(config, trace):
    rms2 = np.power(trace.data, 2).sum()
    rms = np.sqrt(rms2)
    rms_min = config.rmsmin
    if rms <= rms_min:
        raise RuntimeError(
            f'{trace.stats.info}: Trace RMS smaller than {rms_min:g}: '
            'skipping trace')


def _check_clipping(config, trace):
    trace.stats.clipped = False
    if config.clipping_detection_algorithm == 'none':
        return
    # cut the trace between the beginning of noise window
    # and the end of the signal window
    t1 = trace.stats.arrivals['N1'][1]
    if config.wave_type[0] == 'S':
        t2 = trace.stats.arrivals['S2'][1]
    elif config.wave_type[0] == 'P':
        t2 = trace.stats.arrivals['P2'][1]
    tr = trace.copy().trim(t1, t2)
    if config.clipping_detection_algorithm == 'clipping_score':
        score = compute_clipping_score(
            tr, config.remove_baseline, config.clipping_debug_plot)
        logger.info(f'{tr.stats.info}: clipping score: {score:.1f}%')
        if score > config.clipping_score_threshold:
            trace.stats.clipped = True
            trace.stats.ignore = True
            trace.stats.ignore_reason = f'clipping: {score:.1f}%'
    elif config.clipping_detection_algorithm == 'clipping_peaks':
        trace_clipped, properties = clipping_peaks(
            tr,
            config.clipping_peaks_sensitivity,
            config.clipping_peaks_percentile,
            config.clipping_debug_plot)
        logger.info(
            f'{tr.stats.info}: total peaks: {properties["npeaks"]}, '
            f'clipped peaks: {properties["npeaks_clipped"]}')
        if trace_clipped:
            trace.stats.clipped = True
            trace.stats.ignore = True
            trace.stats.ignore_reason = 'clipped'
    if trace.stats.clipped:
        logger.warning(
            f'{tr.stats.info}: Trace is clipped or significantly distorted: '
            'skipping trace')


def _check_sn_ratio(config, trace):
    trace_noise = _get_detrended_trace_copy(trace)
    t1 = trace_noise.stats.arrivals['N1'][1]
    t2 = trace_noise.stats.arrivals['N2'][1]
    trace_noise.trim(starttime=t1, endtime=t2, pad=True, fill_value=0)
    trace_signal = _get_detrended_trace_copy(trace)
    if config.wave_type[0] == 'S':
        t1 = trace_signal.stats.arrivals['S1'][1]
        t2 = trace_signal.stats.arrivals['S2'][1]
    elif config.wave_type[0] == 'P':
        t1 = trace_signal.stats.arrivals['P1'][1]
        t2 = trace_signal.stats.arrivals['P2'][1]
    trace_signal.trim(starttime=t1, endtime=t2, pad=True, fill_value=0)
    rmsnoise2 = np.power(trace_noise.data, 2).sum()
    rmsnoise = np.sqrt(rmsnoise2)
    rmsS2 = np.power(trace_signal.data, 2).sum()
    rmsS = np.sqrt(rmsS2)
    if rmsnoise == 0:
        if config.weighting == 'noise':
            raise RuntimeError(
                f'{trace.stats.info}: empty noise window: skipping trace')
        logger.warning(f'{trace.stats.info}: empty noise window!')
        rmsnoise = 1.
    sn_ratio = rmsS / rmsnoise
    logger.info(f'{trace.stats.info}: S/N: {sn_ratio:.1f}')
    trace.stats.sn_ratio = sn_ratio
    # stop here if trace is already ignored
    if trace.stats.ignore:
        return
    snratio_min = config.sn_min
    if sn_ratio < snratio_min:
        logger.warning(
            f'{trace.stats.info}: S/N smaller than {snratio_min:g}: '
            'skipping trace')
        trace.stats.ignore = True
        trace.stats.ignore_reason = 'low S/N'


def _get_detrended_trace_copy(trace):
    # noise time window for s/n ratio
    tr_copy = trace.copy()
    # remove the mean...
    tr_copy.detrend(type='constant')
    # ...and the linear trend...
    tr_copy.detrend(type='linear')
    return tr_copy


def _remove_baseline(config, trace):
    """
    Get the signal baseline using a Savitzky-Golay filter and subtract it
    from the trace.
    """
    if not config.remove_baseline:
        return
    npts = len(trace.data)
    wlen = npts // 10
    if wlen % 2 == 0:
        wlen += 1
    baseline = savgol_filter(trace.data, wlen, 3)
    trace.data -= baseline


def _process_trace(config, trace):
    # copy trace for manipulation
    trace_process = trace.copy()
    comp = trace_process.stats.channel
    if config.ignore_vertical and comp[-1] in config.vertical_channel_codes:
        if config.wave_type == 'P':
            logger.warning(
                f'{trace.stats.info}: cannot ignore vertical trace, '
                'since "wave_type" is set to "P"')
        else:
            logger.info(
                f'{trace.stats.info}: ignoring vertical trace as requested')
            trace_process.stats.ignore = True
            trace_process.stats.ignore_reason = 'vertical'
    # check if the trace has (significant) signal
    _check_signal_level(config, trace_process)
    # check if trace is clipped
    _check_clipping(config, trace_process)
    # Remove instrument response
    bp_freqmin, bp_freqmax = _get_bandpass_frequencies(config, trace)
    if config.correct_instrumental_response:
        try:
            pre_filt = (
                bp_freqmin, bp_freqmin * 1.1,
                bp_freqmax * 0.9, bp_freqmax)
            remove_instr_response(trace_process, pre_filt)
        except Exception as e:
            raise RuntimeError(
                f'{e}\n'
                f'{trace.stats.info}: Unable to remove instrument response: '
                'skipping trace') from e
    _remove_baseline(config, trace_process)
    filter_trace(config, trace_process)
    # Check if the trace has significant signal to noise ratio
    _check_sn_ratio(config, trace_process)
    return trace_process


def _add_station_to_event_position(trace):
    """
    Add to ``trace.stats`` station-to-event distance (hypocentral and
    epicentral), great-circle distance, azimuth and backazimuth.
    Raise RuntimeError if unable to compute distances.
    """
    try:
        station_to_event_position(trace)
    except Exception as e:
        raise RuntimeError(
            f'{trace.id}: Unable to compute hypocentral distance: {e}. '
            'Skipping trace'
        ) from e


def _check_epicentral_distance(config, trace):
    """
    Reject traces with hypocentral distance outside the range specified
    in the configuration file.
    """
    if config.epi_dist_ranges is None:
        return
    # transform integers to true integers, for better string representation
    edr = [int(r) if r.is_integer() else r for r in config.epi_dist_ranges]
    # if edr has an odd number of elements, add one last element
    if len(edr) % 2 == 1:
        edr.append(999999)
    # reshape edr to a list of pairs
    edr = [[edr[i], edr[i+1]] for i in range(0, len(edr), 2)]
    # string representation of edr
    edr_str = ', '.join(f'{ed[0]}-{ed[1]} km' for ed in edr)
    epi_dist = trace.stats.epi_dist
    # check if epi_dist is in one of the ranges
    if not any(ed[0] <= epi_dist <= ed[1] for ed in edr):
        raise RuntimeError(
            f'{trace.id}: Epicentral distance ({epi_dist:.1f} km) '
            f'not in the selected range ({edr_str}): skipping trace')


def _add_arrivals(config, trace):
    """Add to trace P and S arrival times, travel times and angles."""
    for phase in 'P', 'S':
        try:
            add_arrival_to_trace(trace, phase, config)
        except Exception as e:
            for line in str(e).splitlines():
                logger.warning(line)
            raise RuntimeError(
                f'{trace.id}: Unable to get {phase} arrival time: '
                'skipping trace'
            ) from e


def _define_signal_and_noise_windows(config, trace):
    """Define signal and noise windows for spectral analysis."""
    p_arrival_time = trace.stats.arrivals['P'][1]
    if config.wave_type[0] == 'P' and p_arrival_time < trace.stats.starttime:
        raise RuntimeError(f'{trace.id}: P-window incomplete: skipping trace')
    s_arrival_time = trace.stats.arrivals['S'][1]
    if config.wave_type[0] == 'S' and s_arrival_time < trace.stats.starttime:
        raise RuntimeError(f'{trace.id}: S-window incomplete: skipping trace')
    # Signal window for spectral analysis (S phase)
    s_minus_p = s_arrival_time - p_arrival_time
    s_pre_time = config.signal_pre_time
    if s_minus_p / 2 < s_pre_time:
        # use (Ts-Tp)/2 if it is smaller than signal_pre_time
        # (for short-distance records with short S-P interval)
        s_pre_time = s_minus_p / 2
        logger.warning(
            f'{trace.id}: signal_pre_time is larger than (Ts-Tp)/2. '
            'Using (Ts-Tp)/2 instead')
    t1 = s_arrival_time - s_pre_time
    t1 = max(trace.stats.starttime, t1)
    t2 = t1 + config.win_length
    trace.stats.arrivals['S1'] = ('S1', t1)
    trace.stats.arrivals['S2'] = ('S2', t2)
    # Signal window for spectral analysis (P phase)
    t1 = p_arrival_time - config.signal_pre_time
    t1 = max(trace.stats.starttime, t1)
    t2 = t1 + min(config.win_length, s_minus_p + s_pre_time)
    trace.stats.arrivals['P1'] = ('P1', t1)
    trace.stats.arrivals['P2'] = ('P2', t2)
    # Noise window for spectral analysis
    t1 = max(trace.stats.starttime, p_arrival_time - config.noise_pre_time)
    t2 = t1 + config.win_length
    if t2 >= p_arrival_time:
        logger.warning(f'{trace.id}: noise window ends after P-wave arrival')
        # Note: maybe we should also take into account signal_pre_time here
        t2 = p_arrival_time
        t1 = min(t1, t2)
    trace.stats.arrivals['N1'] = ('N1', t1)
    trace.stats.arrivals['N2'] = ('N2', t2)


def _check_signal_window(config, st):
    """
    Check if the signal window has sufficient amount of signal
    (i.e., not too many gaps).

    This is done on the stream, before merging.
    """
    traceid = st[0].id
    st_cut = st.copy()
    if config.wave_type[0] == 'S':
        t1 = st[0].stats.arrivals['S1'][1]
        t2 = st[0].stats.arrivals['S2'][1]
    elif config.wave_type[0] == 'P':
        t1 = st[0].stats.arrivals['P1'][1]
        t2 = st[0].stats.arrivals['P2'][1]
    st_cut.trim(starttime=t1, endtime=t2)
    if not st_cut:
        raise RuntimeError(
            f'{traceid}: no signal for the selected '
            f'{config.wave_type[0]}-wave cut interval: skipping trace >\n'
            f'> Cut interval: {t1} - {t2}')
    gaps_olaps = st_cut.get_gaps()
    gaps = [g for g in gaps_olaps if g[6] >= 0]
    overlaps = [g for g in gaps_olaps if g[6] < 0]
    duration = st_cut[-1].stats.endtime - st_cut[0].stats.starttime
    gap_duration = sum(g[6] for g in gaps)
    if gap_duration > duration / 4:
        raise RuntimeError(
            f'{traceid}: too many gaps for the selected cut interval: '
            'skipping trace')
    overlap_duration = -1 * sum(g[6] for g in overlaps)
    if overlap_duration > 0:
        logger.info(
            f'{traceid}: signal window has {overlap_duration:.3f} seconds '
            'of overlaps.')


def _merge_stream(config, st):
    """
    Check for gaps and overlaps; remove mean; merge stream.
    """
    traceid = st[0].id
    # First, compute gap/overlap statistics for the whole trace.
    gaps_olaps = st.get_gaps()
    gaps = [g for g in gaps_olaps if g[6] >= 0]
    overlaps = [g for g in gaps_olaps if g[6] < 0]
    gap_duration = sum(g[6] for g in gaps)
    if gap_duration > 0:
        logger.info(
            f'{traceid}: trace has {gap_duration:.3f} seconds of gaps.')
        gap_max = config.gap_max
        if gap_max is not None and gap_duration > gap_max:
            raise RuntimeError(
                f'{traceid}: Gap duration larger than gap_max '
                f'({gap_max:.1f} s): skipping trace')
    overlap_duration = -1 * sum(g[6] for g in overlaps)
    if overlap_duration > 0:
        logger.info(
            f'{traceid}: trace has {overlap_duration:.3f} seconds of '
            'overlaps.')
        overlap_max = config.overlap_max
        if overlap_max is not None and overlap_duration > overlap_max:
            raise RuntimeError(
                f'{traceid}: overlap duration larger than overlap_max '
                f'({overlap_max:.1f} s): skipping trace')
    # Finally, demean (only if trace has not be already preprocessed)
    if config.trace_units == 'auto':
        # Since the count value is generally huge, we need to demean twice
        # to take into account for the rounding error
        st.detrend(type='constant')
        st.detrend(type='constant')
    # Merge stream to remove gaps and overlaps
    try:
        st.merge(fill_value=0)
        # st.merge raises a generic Exception if traces have
        # different sampling rates
    except Exception as e:
        raise RuntimeError(
            f'{traceid}: unable to fill gaps: skipping trace'
        ) from e
    return st[0]


def _skip_ignored(config, traceid):
    """Skip traces ignored from config."""
    network, station, location, channel = traceid.split('.')
    # build a list of all possible ids, from station only
    # to full net.sta.loc.chan
    ss = [
        station,
        '.'.join((network, station)),
        '.'.join((network, station, location)),
        '.'.join((network, station, location, channel)),
    ]
    if config.use_traceids is not None:
        combined = (
            "(" + ")|(".join(config.use_traceids) + ")"
        ).replace('.', r'\.')
        if not any(re.match(combined, s) for s in ss):
            raise RuntimeError(f'{traceid}: ignored from config file')
    if config.ignore_traceids is not None:
        combined = (
            "(" + ")|(".join(config.ignore_traceids) + ")"
        ).replace('.', r'\.')
        if any(re.match(combined, s) for s in ss):
            raise RuntimeError(f'{traceid}: ignored from config file')


def process_traces(config, st):
    """Remove mean, deconvolve and ignore unwanted components."""
    logger.info('Processing traces...')
    out_st = Stream()
    for traceid in sorted({tr.id for tr in st}):
        try:
            _skip_ignored(config, traceid)
            # We still use a stream, since the trace can have gaps or overlaps
            st_sel = st.select(id=traceid)
            for _trace in st_sel:
                _add_station_to_event_position(_trace)
                _check_epicentral_distance(config, _trace)
                _add_arrivals(config, _trace)
                _define_signal_and_noise_windows(config, _trace)
            _check_signal_window(config, st_sel)
            trace = _merge_stream(config, st_sel)
            trace.stats.ignore = False
            trace_process = _process_trace(config, trace)
            out_st.append(trace_process)
        except (ValueError, RuntimeError) as msg:
            logger.warning(msg)
            continue

    if len(out_st) == 0:
        logger.error('No traces left! Exiting.')
        ssp_exit()

    # Rotate traces, if SH or SV is requested
    if config.wave_type in ['SH', 'SV']:
        for traceid in sorted({tr.id[:-1] for tr in out_st}):
            net, sta, loc, chan = traceid.split('.')
            st_sel = out_st.select(
                network=net, station=sta, location=loc, channel=f'{chan}?'
            )
            t0 = max(tr.stats.starttime for tr in st_sel)
            t1 = min(tr.stats.endtime for tr in st_sel)
            st_sel.trim(t0, t1)
            st_sel.rotate('NE->RT')

    logger.info('Processing traces: done')
    logger.info('---------------------------------------------------')
    return out_st
