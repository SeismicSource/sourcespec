# -*- coding: utf-8 -*-
"""
Trace processing for sourcespec.

:copyright:
    2012 Claudio Satriano <satriano@ipgp.fr>

    2013-2014 Claudio Satriano <satriano@ipgp.fr>,
              Emanuela Matrullo <matrullo@geologie.ens.fr>

    2015-2021 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement, Version 2.1
    (http://www.cecill.info/index.en.html)
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import logging
import numpy as np
import re
from obspy.core import Stream
from sourcespec.ssp_setup import ssp_exit
from sourcespec.ssp_util import remove_instr_response, hypo_dist
from sourcespec.ssp_wave_arrival import wave_arrival
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
            logger.warning('%s: Unknown instrument type: %s: '
                           'skipping trace' % (trace.id, instrtype))
            raise ValueError
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
        msg = '%s: maximum frequency for bandpass filtering ' % trace.id
        msg += 'is larger or equal to Nyquist. '
        msg += 'Setting it to %s Hz' % bp_freqmax
        logger.warning(msg)
    trace.filter(type='bandpass', freqmin=bp_freqmin, freqmax=bp_freqmax)


def _check_signal_level(config, trace):
    rms2 = np.power(trace.data, 2).sum()
    rms = np.sqrt(rms2)
    rms_min = config.rmsmin
    if rms <= rms_min:
        logger.warning('%s %s: Trace RMS smaller than %g: '
                       'skipping trace' % (
                        trace.id, trace.stats.instrtype, rms_min))
        raise RuntimeError


def _check_clipping(config, trace):
    clip_tolerance = config.clip_tolerance
    clip_max = (1 - clip_tolerance/100.) * trace.data.max()
    clip_min = (1 - clip_tolerance/100.) * trace.data.min()
    nclips = (trace.data >= clip_max).sum()
    nclips += (trace.data <= clip_min).sum()
    clip_nmax = config.clip_nmax
    if float(nclips)/trace.stats.npts > clip_nmax/100.:
        logger.warning('%s %s: Trace is clipped for more than %.2f%% '
                       'with %.2f%% tolerance: skipping trace' %
                       (trace.id, trace.stats.instrtype, clip_nmax,
                        clip_tolerance))
        raise RuntimeError


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

    # S window for s/n ratio
    trace_cutS = trace.copy()
    # remove the mean...
    trace_cutS.detrend(type='constant')
    # ...and the linear trend...
    trace_cutS.detrend(type='linear')
    t1 = trace_cutS.stats.arrivals['S1'][1]
    t2 = trace_cutS.stats.arrivals['S2'][1]
    trace_cutS.trim(starttime=t1, endtime=t2, pad=True, fill_value=0)

    rmsnoise2 = np.power(trace_noise.data, 2).sum()
    rmsnoise = np.sqrt(rmsnoise2)
    rmsS2 = np.power(trace_cutS.data, 2).sum()
    rmsS = np.sqrt(rmsS2)

    if rmsnoise == 0:
        logger.warning('%s %s: empty noise window: skipping trace' %
                       (trace.id, trace.stats.instrtype))
        raise RuntimeError
    sn_ratio = rmsS/rmsnoise
    logger.info('%s %s: S/N: %.1f' % (
        trace.id, trace.stats.instrtype, sn_ratio))
    trace.stats.sn_ratio = sn_ratio

    snratio_min = config.sn_min
    if sn_ratio < snratio_min:
        logger.warning('%s %s: S/N smaller than %g: ignoring trace' %
                       (trace.id, trace.stats.instrtype, snratio_min))
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
            logger.warning('%s %s: Unable to remove instrument response: '
                           'skipping trace' % (trace_process.id, instrtype))
            raise RuntimeError

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
        logger.info('%s: trace has %.3f seconds of gaps.' %
                    (traceid, gap_duration))
        gap_max = config.gap_max
        if gap_max is not None and gap_duration > gap_max:
            logger.warning('%s: Gap duration larger than gap_max (%.1f): '
                           'skipping trace' % (traceid, gap_max))
            raise RuntimeError
    overlap_duration = -1 * sum(g[6] for g in overlaps)
    if overlap_duration > 0:
        logger.info('%s: trace has %.3f seconds of overlaps.' %
                    (traceid, overlap_duration))
        overlap_max = config.overlap_max
        if overlap_max is not None and overlap_duration > overlap_max:
            logger.warning('%s: Overlap duration larger than '
                           'overlap_max (%.1f): skipping trace' %
                           (traceid, overlap_max))
            raise RuntimeError
    # Then, compute the same statisics for the S-wave window.
    st_cut = st.copy()
    t1 = st[0].stats.arrivals['S1'][1]
    t2 = st[0].stats.arrivals['S2'][1]
    st_cut.trim(starttime=t1, endtime=t2)
    if not st_cut:
        logger.warning('%s: Not enough signal for the selected cut '
                       'interval: skipping trace' % traceid)
        raise RuntimeError
    gaps_olaps = st_cut.get_gaps()
    gaps = [g for g in gaps_olaps if g[6] >= 0]
    overlaps = [g for g in gaps_olaps if g[6] < 0]
    duration = st_cut[-1].stats.endtime - st_cut[0].stats.starttime
    gap_duration = sum(g[6] for g in gaps)
    if gap_duration > duration/4:
        logger.warning('%s: Too many gaps for the selected cut '
                       'interval: skipping trace' % traceid)
        raise RuntimeError
    overlap_duration = -1 * sum(g[6] for g in overlaps)
    if overlap_duration > 0:
        logger.info('%s: S-wave window has %.3f seconds of overlaps.' %
                    (traceid, overlap_duration))
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
        raise RuntimeError
    return st[0]


def _add_hypo_dist_and_arrivals(config, st):
    for trace in st:
        if hypo_dist(trace) is None:
            logger.warning('%s: Unable to compute hypocentral distance: '
                           'skipping trace' % trace.id)
            raise RuntimeError
        if config.max_epi_dist is not None and \
                trace.stats.epi_dist > config.max_epi_dist:
            logger.warning('%s: Epicentral distance (%.1f) '
                           'larger than max_epi_dist (%.1f): skipping trace' %
                           (trace.id, trace.stats.epi_dist,
                            config.max_epi_dist))
            raise RuntimeError

        p_arrival_time = wave_arrival(trace, 'P', config.p_arrival_tolerance,
                                      config.vp_tt, config.NLL_time_dir)
        s_arrival_time = wave_arrival(trace, 'S', config.s_arrival_tolerance,
                                      config.vs_tt, config.NLL_time_dir)
        if p_arrival_time is None or s_arrival_time is None:
            logger.warning('%s: Unable to get arrival times: '
                           'skipping trace' % trace.id)
            raise RuntimeError
        # Signal window for spectral analysis
        t1 = s_arrival_time - config.pre_s_time
        t2 = t1 + config.win_length
        trace.stats.arrivals['S1'] = ('S1', t1)
        trace.stats.arrivals['S2'] = ('S2', t2)
        # Noise window for spectral analysis
        t1 = p_arrival_time - config.pre_p_time
        t2 = t1 + config.win_length
        trace.stats.arrivals['N1'] = ('N1', t1)
        trace.stats.arrivals['N2'] = ('N2', t2)


def process_traces(config, st):
    """Remove mean, deconvolve and ignore unwanted components."""
    out_st = Stream()
    for id in sorted(set(tr.id for tr in st)):
        # We still use a stream, since the trace can have
        # gaps or overlaps
        st_sel = st.select(id=id)
        network, station, location, channel = id.split('.')
        # build a list of all possible ids, from station only
        # to full net.sta.loc.chan
        ss = [station, ]
        ss.append('.'.join((network, station)))
        ss.append('.'.join((network, station, location)))
        ss.append('.'.join((network, station, location, channel)))
        if config.use_stations is not None:
            combined = (
                "(" + ")|(".join(config.use_stations) + ")"
                ).replace('.', '\.')
            if not any(re.match(combined, s) for s in ss):
                logger.warning('%s: ignored from config file' % id)
                continue
        if config.ignore_stations is not None:
            combined = (
                "(" + ")|(".join(config.ignore_stations) + ")"
                ).replace('.', '\.')
            if any(re.match(combined, s) for s in ss):
                logger.warning('%s: ignored from config file' % id)
                continue
        try:
            _add_hypo_dist_and_arrivals(config, st_sel)
            trace = _merge_stream(config, st_sel)
            trace.stats.ignore = False
            trace_process = _process_trace(config, trace)
            out_st.append(trace_process)
        except (ValueError, RuntimeError):
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

    return out_st
