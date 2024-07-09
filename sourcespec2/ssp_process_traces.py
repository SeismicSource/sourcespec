# -*- coding: utf-8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Trace processing for sourcespec.

:copyright:
    2012 Claudio Satriano <satriano@ipgp.fr>

    2013-2014 Claudio Satriano <satriano@ipgp.fr>,
              Emanuela Matrullo <matrullo@geologie.ens.fr>

    2015-2026 Claudio Satriano <satriano@ipgp.fr>
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
from .config import config
from .ssp_setup import ssp_exit
from .ssp_util import (
    remove_instr_response, station_to_event_position)
from .ssp_wave_arrival import add_arrival_to_trace
from .ssp_wave_picking import refine_trace_picks
from .clipping_detection import (
    check_min_amplitude, compute_clipping_score, clipping_peaks)
logger = logging.getLogger(__name__.rsplit('.', maxsplit=1)[-1])


def _skip_trace_and_raise(trace, reason, short_reason=None):
    """Skip a trace with a reason and raise a RuntimeError."""
    trace.stats.ignore = True
    trace.stats.ignore_reason = short_reason or reason
    raise RuntimeError(f'{trace.id}: {reason}: skipping trace')


def _skip_stream_and_raise(st, reason, short_reason=None):
    """Skip a stream with a reason and raise a RuntimeError."""
    for trace in st:
        trace.stats.ignore = True
        trace.stats.ignore_reason = short_reason or reason
    raise RuntimeError(f'{st[0].id}: {reason}: skipping trace')


def _get_bandpass_frequencies(trace):
    """
    Get frequencies for bandpass filter.

    :param trace: ObsPy Trace object.
    :type trace: :class:`obspy.core.trace.Trace`

    :return: Tuple with minimum and maximum frequencies.
    :rtype: tuple

    :raises: ValueError if instrument type is unknown.
    """
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
        except KeyError:
            _skip_trace_and_raise(
                trace,
                reason=f'Unknown instrument type: {instrtype}',
                short_reason='unknown instr type'
            )
    return bp_freqmin, bp_freqmax


def filter_trace(trace):
    """
    Filter trace in place.

    Save filter info to trace stats.

    :param trace: ObsPy Trace object.
    :type trace: :class:`obspy.core.trace.Trace`
    """
    bp_freqmin, bp_freqmax = _get_bandpass_frequencies(trace)
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


def _check_signal_level(trace):
    """
    Check if the trace has significant signal level.

    :param trace: ObsPy Trace object.
    :type trace: :class:`obspy.core

    :raises: RuntimeError if trace RMS is smaller than config.rmsmin.
    """
    rms2 = np.power(trace.data, 2).sum()
    rms = np.sqrt(rms2)
    rms_min = config.rmsmin
    if rms <= rms_min:
        _skip_trace_and_raise(
            trace,
            reason=f'RMS smaller than {rms_min:g}',
            short_reason='low RMS'
        )


def _check_clipping(trace):
    """
    Check if the trace is clipped.

    :param trace: ObsPy Trace object.
    :type trace: :class:`obspy.core.trace.Trace`
    """
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
    else:
        # this should never happen
        raise ValueError(f'Unknown wave type: {config.wave_type[0]}')
    tr = trace.copy().trim(t1, t2)
    min_ampl_test, tr_max, min_thresh = check_min_amplitude(
        tr, config.clipping_min_amplitude_ratio)
    if not min_ampl_test:
        logger.info(
            f'{trace.id}: max amplitude ({tr_max:.1f}) below minimum '
            f'threshold ({min_thresh:.1f}): skipping clipping check'
        )
        return
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


def _check_sn_ratio(trace):
    """
    Check if the trace has significant signal to noise ratio.

    :param trace: ObsPy Trace object.
    :type trace: :class:`obspy.core.trace.Trace`

    :raises: RuntimeError if no noise window is available.
    """
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
            _skip_trace_and_raise(
                trace,
                reason='empty noise window',
                short_reason='empty noise window'
            )
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
    """
    Get a copy of the trace with the mean and linear trend removed.

    :param trace: ObsPy Trace object.
    :type trace: :class:`obspy.core.trace.Trace`

    :return: Trace with mean and linear trend removed.
    :rtype: :class:`obspy.core.trace.Trace`
    """
    # noise time window for s/n ratio
    tr_copy = trace.copy()
    # remove the mean...
    tr_copy.detrend(type='constant')
    # ...and the linear trend...
    tr_copy.detrend(type='linear')
    return tr_copy


def _remove_baseline(trace):
    """
    Get the signal baseline using a Savitzky-Golay filter and subtract it
    from the trace.

    :param trace: ObsPy Trace object.
    :type trace: :class:`obspy.core.trace.Trace`
    """
    if not config.remove_baseline:
        return
    npts = len(trace.data)
    wlen = npts // 10
    if wlen % 2 == 0:
        wlen += 1
    baseline = savgol_filter(trace.data, wlen, 3)
    trace.data -= baseline


def _process_trace(trace):
    """
    Process a trace.

    :param trace: ObsPy Trace object.
    :type trace: :class:`obspy.core.trace.Trace`

    :return: Processed trace.
    :rtype: :class:`obspy.core.trace.Trace`

    :raises: RuntimeError if unable to remove instrument response.
    """
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
    _check_signal_level(trace_process)
    # check if trace is clipped
    _check_clipping(trace_process)
    # Remove instrument response
    bp_freqmin, bp_freqmax = _get_bandpass_frequencies(trace)
    if config.correct_instrumental_response:
        try:
            pre_filt = (
                bp_freqmin, bp_freqmin * 1.1,
                bp_freqmax * 0.9, bp_freqmax)
            remove_instr_response(trace_process, pre_filt)
        except Exception as e:
            _skip_trace_and_raise(
                trace,
                reason=f'unable to remove instrument response: {e}',
                short_reason='no instr response'
            )
    _remove_baseline(trace_process)
    filter_trace(trace_process)
    # Check if the trace has significant signal to noise ratio
    _check_sn_ratio(trace_process)
    trace_process.stats.processed = True
    return trace_process


def _add_station_to_event_position(trace):
    """
    Add to ``trace.stats`` station-to-event distance (hypocentral and
    epicentral), great-circle distance, azimuth and backazimuth.
    Raise RuntimeError if unable to compute distances.

    :param trace: ObsPy Trace object.
    :type trace: :class:`obspy.core.trace.Trace`

    :raises: RuntimeError if unable to compute hypocentral distance.
    """
    try:
        station_to_event_position(trace)
    except Exception as e:
        _skip_trace_and_raise(
            trace,
            reason=f'unable to compute hypocentral distance: {e}',
            short_reason='no hypocentral distance'
        )


def _check_epicentral_distance(trace):
    """
    Check if the epicentral distance is within the selected range.

    :param trace: ObsPy Trace object.
    :type trace: :class:`obspy.core.trace.Trace`

    :raises: RuntimeError if epicentral distance is outside the range.
    """
    if config.epi_dist_ranges is None:
        return
    # transform integers to true integers, for better string representation
    edr = [
        int(r) if float(r).is_integer() else r for r in config.epi_dist_ranges]
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
        _skip_trace_and_raise(
            trace,
            reason=(
                f'Epicentral distance ({epi_dist:.1f} km) '
                f'not in the selected range ({edr_str})'
            ),
            short_reason='epi dist outside range'
        )


def _add_arrivals(trace):
    """
    Add to trace P and S arrival times, travel times and angles.

    :param trace: ObsPy Trace object.
    :type trace: :class:`obspy.core.trace.Trace`

    :raises: RuntimeError if unable to get arrival times.
    """
    for phase in 'P', 'S':
        try:
            add_arrival_to_trace(trace, phase)
        except Exception as e:
            for line in str(e).splitlines():
                logger.warning(line)
            _skip_trace_and_raise(
                trace,
                reason=f'unable to get {phase} arrival times',
                short_reason=f'no {phase} arrival times',
            )
    if config.refine_theoretical_arrivals:
        freqmin = config.autopick_freqmin
        debug = config.autopick_debug_plot
        refine_trace_picks(trace, freqmin, debug)


def _define_signal_and_noise_windows(trace):
    """
    Define signal and noise windows for spectral analysis.

    :param trace: ObsPy Trace object.
    :type trace: :class:`obspy.core.trace.Trace`

    :raises: RuntimeError if P or S window is incomplete.
    """
    p_arrival_time = trace.stats.arrivals['P'][1]
    if config.wave_type[0] == 'P' and p_arrival_time < trace.stats.starttime:
        _skip_trace_and_raise(trace, 'P-window incomplete')
    s_arrival_time = trace.stats.arrivals['S'][1]
    if config.wave_type[0] == 'S' and s_arrival_time < trace.stats.starttime:
        _skip_trace_and_raise(trace, 'S-window incomplete')
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
    tt_s = trace.stats.travel_times['S']
    # manage case where variable_win_length_factor is None:
    variable_win_length_factor = config.variable_win_length_factor or 0
    win_length_s = max(config.win_length, variable_win_length_factor * tt_s)
    t2 = t1 + win_length_s
    t2 = min(trace.stats.endtime, t2)
    win_length_s = t2 - t1
    trace.stats.arrivals['S1'] = ('S1', t1)
    trace.stats.arrivals['S2'] = ('S2', t2)
    # Signal window for spectral analysis (P phase)
    t1 = p_arrival_time - config.signal_pre_time
    t1 = max(trace.stats.starttime, t1)
    tt_p = trace.stats.travel_times['P']
    win_length_p = max(config.win_length, variable_win_length_factor * tt_p)
    t2 = t1 + min(win_length_p, s_minus_p + s_pre_time)
    t2 = min(trace.stats.endtime, t2)
    win_length_p = t2 - t1
    trace.stats.arrivals['P1'] = ('P1', t1)
    trace.stats.arrivals['P2'] = ('P2', t2)
    # Noise window for spectral analysis
    if config.wave_type[0] == 'P':
        win_length = win_length_p
    elif config.wave_type[0] == 'S':
        win_length = win_length_s
    if config.variable_win_length_factor:
        logger.info(f'{trace.id}: window length {win_length:.3f} seconds')
    if config.noise_pre_time is None:
        noise_pre_time = win_length + config.signal_pre_time
        logger.info(
            f'{trace.id}: noise_pre_time autoset to '
            f'{noise_pre_time:.3f} seconds')
    else:
        noise_pre_time = config.noise_pre_time
    t1 = max(
        trace.stats.starttime,
        p_arrival_time - noise_pre_time
    )
    t2 = t1 + win_length
    if t2 > (p_arrival_time - config.signal_pre_time):
        logger.warning(f'{trace.id}: noise window ends after P-window start')
        t2 = p_arrival_time - config.signal_pre_time
        t1 = min(t1, t2)
    trace.stats.arrivals['N1'] = ('N1', t1)
    trace.stats.arrivals['N2'] = ('N2', t2)


def _check_signal_window(st):
    """
    Check if the signal window has sufficient amount of signal
    (i.e., not too many gaps).

    This is done on the stream, before merging.

    :param st: ObsPy Stream object.
    :type st: :class:`obspy.core.stream.Stream`

    :raises: RuntimeError if the cut interval has no signal.
    :raises: RuntimeError if too many gaps in the signal window.
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
        _skip_stream_and_raise(
            st,
            reason=(
                'no signal for the selected '
                f'{config.wave_type[0]}-wave cut interval ({t1} - {t2})'
            ),
            short_reason='no signal'
        )
    gaps_olaps = st_cut.get_gaps()
    gaps = [g for g in gaps_olaps if g[6] >= 0]
    overlaps = [g for g in gaps_olaps if g[6] < 0]
    duration = st_cut[-1].stats.endtime - st_cut[0].stats.starttime
    gap_duration = sum(g[6] for g in gaps)
    if gap_duration > duration / 4:
        _skip_stream_and_raise(
            st,
            reason='too many gaps for the selected cut interval',
            short_reason='too many gaps'
        )
    overlap_duration = -1 * sum(g[6] for g in overlaps)
    if overlap_duration > 0:
        logger.info(
            f'{traceid}: signal window has {overlap_duration:.3f} seconds '
            'of overlaps.')


def _merge_stream(st):
    """
    Check for gaps and overlaps; remove mean; merge stream.

    :param st: ObsPy Stream object.
    :type st: :class:`obspy.core.stream.Stream`

    :return: An ObsPy Trace object.
    :rtype: :class:`obspy.core.trace.Trace`

    :raises: RuntimeError if gap duration is larger than config.gap_max.
    :raises: RuntimeError if overlap duration is larger than
        config.overlap_max.
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
            _skip_stream_and_raise(
                st,
                reason=f'Gap duration larger than {gap_max:.1f} s',
                short_reason='gap too long'
            )
    overlap_duration = -1 * sum(g[6] for g in overlaps)
    if overlap_duration > 0:
        logger.info(
            f'{traceid}: trace has {overlap_duration:.3f} seconds of '
            'overlaps.')
        overlap_max = config.overlap_max
        if overlap_max is not None and overlap_duration > overlap_max:
            _skip_stream_and_raise(
                st,
                reason=f'Overlap duration larger than {overlap_max:.1f} s',
                short_reason='overlap too long'
            )
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
    except Exception:
        _skip_stream_and_raise(st, 'unable to fill gaps')
    return st[0]


def _glob_to_regex(pattern):
    """
    Convert a glob-style pattern to a regex pattern.
    """
    # Escape regex-special characters except for ? and *
    pattern = pattern.strip()
    pattern = re.escape(pattern)            # Escapes all regex chars
    pattern = pattern.replace(r'\?', '.')   # Convert escaped ? to .
    pattern = pattern.replace(r'\*', '.*')  # Convert escaped * to .*
    return pattern


def _skip_ignored(st):
    """
    Skip traces ignored from config.

    :param st: ObsPy Stream object.
    :type st: :class:`obspy.core.stream.Stream`

    :raises: RuntimeError if traces are ignored from config file.
    """
    traceid = st[0].id
    network, station, location, channel = traceid.split('.')

    # Build all possible IDs: station → full net.sta.loc.chan
    ss = [
        station,
        '.'.join((network, station)),
        '.'.join((network, station, location)),
        '.'.join((network, station, location, channel)),
    ]

    def _matches(patterns, strings):
        """Return True if any string matches any glob pattern."""
        regex = r'^(' + '|'.join(_glob_to_regex(p) for p in patterns) + r')$'
        return any(re.match(regex, s) for s in strings)

    if (
        config.use_traceids is not None and
        not _matches(config.use_traceids, ss)
    ):
        _skip_stream_and_raise(
            st,
            reason='ignored from config file',
            short_reason='ignored from config'
        )

    if (
        config.ignore_traceids is not None and
        _matches(config.ignore_traceids, ss)
    ):
        _skip_stream_and_raise(
            st,
            reason='ignored from config file',
            short_reason='ignored from config'
        )


def process_traces(st):
    """
    Remove mean, deconvolve and ignore unwanted components.

    :param st: Stream object with traces to process.
    :type st: :class:`obspy.core.stream.Stream`

    :return: Stream object with processed traces.
    :rtype: :class:`obspy.core.stream.Stream`
    """
    logger.info('Processing traces...')
    out_st = Stream()
    for traceid in sorted({tr.id for tr in st}):
        try:
            # We still use a stream, since the trace can have gaps or overlaps
            st_sel = st.select(id=traceid)
            # Add event-related metadata to trace.stats
            for _trace in st_sel:
                _add_station_to_event_position(_trace)
                _check_epicentral_distance(_trace)
                _add_arrivals(_trace)
                _define_signal_and_noise_windows(_trace)
            # We skip traces ignored from config file here, so that we have
            # the metadata needed for the raw plot
            _skip_ignored(st_sel)
            _check_signal_window(st_sel)
            trace = _merge_stream(st_sel)
            trace.stats.ignore = False
            trace_process = _process_trace(trace)
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
