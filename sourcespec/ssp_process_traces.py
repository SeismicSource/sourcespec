# -*- coding: utf-8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Trace processing for sourcespec.

:copyright:
    2012 Claudio Satriano <satriano@ipgp.fr>

    2013-2014 Claudio Satriano <satriano@ipgp.fr>,
              Emanuela Matrullo <matrullo@geologie.ens.fr>

    2015-2022 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
import logging
import numpy as np
import re
from obspy.core import Stream
from sourcespec.ssp_setup import ssp_exit
from sourcespec.ssp_util import remove_instr_response, hypo_dist
from sourcespec.ssp_wave_arrival import add_arrivals_to_trace
logger = logging.getLogger(__name__.split('.')[-1])


def _get_bandpass_frequencies(config, trace):
    """Get frequencies for bandpass filter."""
    # see if there is a station-specfic filter
    station = trace.stats.station
    try:
        bp_freqmin = float(config['bp_freqmin_' + station])
        bp_freqmax = float(config['bp_freqmax_' + station])
    except KeyError:
        instrtype = trace.stats.instrtype
        try:
            bp_freqmin = float(config['bp_freqmin_' + instrtype])
            bp_freqmax = float(config['bp_freqmax_' + instrtype])
        except KeyError:
            msg = '{}: Unknown instrument type: {}: skipping trace'
            msg = msg.format(trace.id, instrtype)
            raise ValueError(msg)
    return bp_freqmin, bp_freqmax


def filter_trace(config, trace):
    bp_freqmin, bp_freqmax = _get_bandpass_frequencies(config, trace)
    # remove the mean...
    trace.detrend(type='constant')
    # ...and the linear trend...
    trace.detrend(type='linear')
    nyquist = 1./(2. * trace.stats.delta)
    if bp_freqmax >= nyquist:
        bp_freqmax = nyquist * 0.999
        msg = '{}: maximum frequency for bandpass filtering '
        msg += 'is larger or equal to Nyquist. Setting it to {} Hz'
        msg = msg.format(trace.id, bp_freqmax)
        logger.warning(msg)
    trace.filter(type='bandpass', freqmin=bp_freqmin, freqmax=bp_freqmax)


def _check_signal_level(config, trace):
    rms2 = np.power(trace.data, 2).sum()
    rms = np.sqrt(rms2)
    rms_min = config.rmsmin
    if rms <= rms_min:
        msg = '{} {}: Trace RMS smaller than {:g}: skipping trace'
        msg = msg.format(trace.id, trace.stats.instrtype, rms_min)
        raise RuntimeError(msg)


def _check_clipping(config, trace):
    t1 = (trace.stats.arrivals['P'][1] - config.noise_pre_time)
    t2 = (trace.stats.arrivals['S'][1] + config.win_length)
    tr = trace.copy().trim(t1, t2).detrend('demean')
    npts = len(tr.data)
    # Compute data histogram with a number of bins equal to 10% of data points
    hist, bins = np.histogram(np.abs(tr.data), bins=int(npts*0.1))
    # Clipped samples are in the last bin
    nclips = hist[-1]
    clip_max_percent = config.clip_max_percent
    if nclips/npts > clip_max_percent/100.:
        msg = (
            '{} {}: Trace is clipped for more than {:.2f}% '
            'skipping trace'.format(
                tr.id, tr.stats.instrtype, clip_max_percent)
        )
        raise RuntimeError(msg)


def _check_sn_ratio(config, trace):
    # noise time window for s/n ratio
    trace_noise = trace.copy()
    # remove the mean...
    trace_noise.detrend(type='constant')
    # ...and the linear trend...
    trace_noise.detrend(type='linear')
    t1 = trace_noise.stats.arrivals['N1'][1]
    t2 = trace_noise.stats.arrivals['N2'][1]
    trace_noise.trim(starttime=t1, endtime=t2, pad=True, fill_value=0)
    # signal window for s/n ratio
    trace_signal = trace.copy()
    # remove the mean...
    trace_signal.detrend(type='constant')
    # ...and the linear trend...
    trace_signal.detrend(type='linear')
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
        msg = '{} {}: empty noise window: skipping trace'
        msg = msg.format(trace.id, trace.stats.instrtype)
        raise RuntimeError(msg)
    sn_ratio = rmsS/rmsnoise
    logger.info('{} {}: S/N: {:.1f}'.format(
        trace.id, trace.stats.instrtype, sn_ratio))
    trace.stats.sn_ratio = sn_ratio
    snratio_min = config.sn_min
    if sn_ratio < snratio_min:
        msg = '{} {}: S/N smaller than {:g}: skipping trace'
        msg = msg.format(trace.id, trace.stats.instrtype, snratio_min)
        logger.warning(msg)
        trace.stats.ignore = True


def _process_trace(config, trace):
    # copy trace for manipulation
    trace_process = trace.copy()
    comp = trace_process.stats.channel
    instrtype = trace_process.stats.instrtype
    if config.ignore_vertical and comp[-1] in ['Z', '1']:
        raise RuntimeError
    # check if the trace has (significant) signal
    _check_signal_level(config, trace_process)
    # check if trace is clipped
    _check_clipping(config, trace_process)
    # Remove instrument response
    if not config.options.no_response:
        correct = config.correct_instrumental_response
        bp_freqmin, bp_freqmax = _get_bandpass_frequencies(config, trace)
        pre_filt = (bp_freqmin, bp_freqmin*1.1, bp_freqmax*0.9, bp_freqmax)
        if remove_instr_response(trace_process, correct, pre_filt) is None:
            msg = '{} {}: Unable to remove instrument response: '
            msg += 'skipping trace'
            msg = msg.format(trace_process.id, instrtype)
            raise RuntimeError(msg)
    filter_trace(config, trace_process)
    # Check if the trace has significant signal to noise ratio
    _check_sn_ratio(config, trace_process)
    return trace_process


def _merge_stream(config, st):
    traceid = st[0].id
    # First, compute gap/overlap statistics for the whole trace.
    gaps_olaps = st.get_gaps()
    gaps = [g for g in gaps_olaps if g[6] >= 0]
    overlaps = [g for g in gaps_olaps if g[6] < 0]
    gap_duration = sum(g[6] for g in gaps)
    if gap_duration > 0:
        msg = '{}: trace has {:.3f} seconds of gaps.'
        msg = msg.format(traceid, gap_duration)
        logger.info(msg)
        gap_max = config.gap_max
        if gap_max is not None and gap_duration > gap_max:
            msg = '{}: Gap duration larger than gap_max ({:.1f} s): '
            msg += 'skipping trace'
            msg = msg.format(traceid, gap_max)
            raise RuntimeError(msg)
    overlap_duration = -1 * sum(g[6] for g in overlaps)
    if overlap_duration > 0:
        msg = '{}: trace has {:.3f} seconds of overlaps.'
        msg = msg.format(traceid, overlap_duration)
        logger.info(msg)
        overlap_max = config.overlap_max
        if overlap_max is not None and overlap_duration > overlap_max:
            msg = '{}: Overlap duration larger than overlap_max ({:.1f} s): '
            msg += 'skipping trace'
            msg = msg.format(traceid, overlap_max)
            raise RuntimeError(msg)
    # Then, compute the same statisics for the signal window.
    st_cut = st.copy()
    if config.wave_type[0] == 'S':
        t1 = st[0].stats.arrivals['S1'][1]
        t2 = st[0].stats.arrivals['S2'][1]
    elif config.wave_type[0] == 'P':
        t1 = st[0].stats.arrivals['P1'][1]
        t2 = st[0].stats.arrivals['P2'][1]
    st_cut.trim(starttime=t1, endtime=t2)
    if not st_cut:
        msg = '{}: No signal for the selected {}-wave cut interval: '
        msg += 'skipping trace >\n'
        msg += '> Cut interval: {} - {}'
        msg = msg.format(traceid, config.wave_type[0], t1, t2)
        raise RuntimeError(msg)
    gaps_olaps = st_cut.get_gaps()
    gaps = [g for g in gaps_olaps if g[6] >= 0]
    overlaps = [g for g in gaps_olaps if g[6] < 0]
    duration = st_cut[-1].stats.endtime - st_cut[0].stats.starttime
    gap_duration = sum(g[6] for g in gaps)
    if gap_duration > duration/4:
        msg = '{}: Too many gaps for the selected cut interval: skipping trace'
        msg = msg.format(traceid)
        raise RuntimeError(msg)
    overlap_duration = -1 * sum(g[6] for g in overlaps)
    if overlap_duration > 0:
        msg = '{}: Signal window has {:.3f} seconds of overlaps.'
        msg = msg.format(traceid, overlap_duration)
        logger.info(msg)
    # Finally, demean and remove gaps and overlaps.
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
        msg = '{}: unable to fill gaps: skipping trace'.format(traceid)
        raise RuntimeError(msg)
    return st[0]


def _add_hypo_dist_and_arrivals(config, st):
    for trace in st:
        if hypo_dist(trace) is None:
            msg = '{}: Unable to compute hypocentral distance: skipping trace'
            msg = msg.format(trace.id)
            raise RuntimeError(msg)
        if config.max_epi_dist is not None and \
                trace.stats.epi_dist > config.max_epi_dist:
            msg = '{}: Epicentral distance ({:.1f} km) '
            msg += 'larger than max_epi_dist ({:.1f} km): skipping trace'
            msg = msg.format(
                trace.id, trace.stats.epi_dist, config.max_epi_dist)
            raise RuntimeError(msg)
        add_arrivals_to_trace(trace, config)
        try:
            p_arrival_time = trace.stats.arrivals['P'][1]
        except KeyError:
            msg = '{}: Unable to get P arrival time: skipping trace'
            msg = msg.format(trace.id)
            raise RuntimeError(msg)
        try:
            s_arrival_time = trace.stats.arrivals['S'][1]
        except KeyError:
            msg = '{}: Unable to get S arrival time: skipping trace'
            msg = msg.format(trace.id)
            raise RuntimeError(msg)
        # Signal window for spectral analysis (S phase)
        t1 = s_arrival_time - config.signal_pre_time
        t2 = t1 + config.win_length
        trace.stats.arrivals['S1'] = ('S1', t1)
        trace.stats.arrivals['S2'] = ('S2', t2)
        # Signal window for spectral analysis (P phase)
        t1 = p_arrival_time - config.signal_pre_time
        t2 = t1 + min(config.win_length, s_arrival_time-p_arrival_time)
        trace.stats.arrivals['P1'] = ('P1', t1)
        trace.stats.arrivals['P2'] = ('P2', t2)
        # Noise window for spectral analysis
        t1 = p_arrival_time - config.noise_pre_time
        t2 = t1 + config.win_length
        if t2 >= s_arrival_time:
            logger.warning(
                '{}: noise window ends after S-wave arrival'.format(trace.id))
        elif t2 >= p_arrival_time:
            logger.warning(
                '{}: noise window ends after P-wave arrival'.format(trace.id))
        trace.stats.arrivals['N1'] = ('N1', t1)
        trace.stats.arrivals['N2'] = ('N2', t2)


def _skip_ignored(config, id):
    """Skip traces ignored from config."""
    network, station, location, channel = id.split('.')
    # build a list of all possible ids, from station only
    # to full net.sta.loc.chan
    ss = [station, ]
    ss.append('.'.join((network, station)))
    ss.append('.'.join((network, station, location)))
    ss.append('.'.join((network, station, location, channel)))
    if config.use_traceids is not None:
        combined = (
            "(" + ")|(".join(config.use_traceids) + ")"
            ).replace('.', r'\.')
        if not any(re.match(combined, s) for s in ss):
            msg = '{}: ignored from config file'.format(id)
            raise RuntimeError(msg)
    if config.ignore_traceids is not None:
        combined = (
            "(" + ")|(".join(config.ignore_traceids) + ")"
            ).replace('.', r'\.')
        if any(re.match(combined, s) for s in ss):
            msg = '{}: ignored from config file'.format(id)
            raise RuntimeError(msg)


def process_traces(config, st):
    """Remove mean, deconvolve and ignore unwanted components."""
    logger.info('Processing traces...')
    out_st = Stream()
    for id in sorted(set(tr.id for tr in st)):
        # We still use a stream, since the trace can have gaps or overlaps
        st_sel = st.select(id=id)
        try:
            _skip_ignored(config, id)
            _add_hypo_dist_and_arrivals(config, st_sel)
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
        for id in sorted(set(tr.id[:-1] for tr in out_st)):
            net, sta, loc, chan = id.split('.')
            st_sel = out_st.select(network=net, station=sta,
                                   location=loc, channel=chan+'?')
            t0 = max(tr.stats.starttime for tr in st_sel)
            t1 = min(tr.stats.endtime for tr in st_sel)
            st_sel.trim(t0, t1)
            st_sel.rotate('NE->RT')

    logger.info('Processing traces: done')
    return out_st
