# -*- coding: utf-8 -*-
"""
Trace processing for sourcespec.

:copyright:
    2012 Claudio Satriano <satriano@ipgp.fr>

    2013-2014 Claudio Satriano <satriano@ipgp.fr>,
              Emanuela Matrullo <matrullo@geologie.ens.fr>

    2015-2016 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement, Version 2.1
    (http://www.cecill.info/index.en.html)
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import logging
import numpy as np
from obspy.core import Stream
from sourcespec.ssp_setup import ssp_exit
from sourcespec.ssp_util import remove_instr_response, hypo_dist, wave_arrival


def filter_trace(config, trace):
    instrtype = trace.stats.instrtype
    if instrtype == 'acc':
        # band-pass frequencies:
        # TODO: calculate from sampling rate?
        bp_freqmin = config.bp_freqmin_acc
        bp_freqmax = config.bp_freqmax_acc
    elif instrtype == 'shortp':
        # band-pass frequencies:
        # TODO: calculate from sampling rate?
        bp_freqmin = config.bp_freqmin_shortp
        bp_freqmax = config.bp_freqmax_shortp
    elif instrtype == 'broadb':
        # band-pass frequencies:
        bp_freqmin = config.bp_freqmin_broadb
        bp_freqmax = config.bp_freqmax_broadb
    else:
        raise ValueError
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
        logging.warning(msg)
    trace.filter(type='bandpass', freqmin=bp_freqmin, freqmax=bp_freqmax)


def _has_not_signal(config, trace):
    rms2 = np.power(trace.data, 2).sum()
    rms = np.sqrt(rms2)
    rms_min = config.rmsmin
    if rms <= rms_min:
        logging.warning('%s %s: Trace RMS smaller than %g: '
                        'skipping trace' % (trace.id, trace.stats.instrtype,
                                            rms_min))
        return True
    else:
        return False


def _is_clipped(config, trace):
    clip_tolerance = config.clip_tolerance
    clip_max = (1 - clip_tolerance/100.) * trace.data.max()
    clip_min = (1 - clip_tolerance/100.) * trace.data.min()
    nclips = (trace.data > clip_max).sum()
    nclips += (trace.data < clip_min).sum()
    clip_nmax = config.clip_nmax
    if float(nclips)/trace.stats.npts > clip_nmax/100.:
        logging.warning('%s %s: Trace is clipped for more than %.2f%% '
                        'with %.2f%% tolerance: skipping trace' %
                        (trace.id, trace.stats.instrtype, clip_nmax,
                         clip_tolerance))
        return True
    else:
        return False


def _has_low_sn_ratio(config, trace):
    # noise time window for s/n ratio
    p_arrival_time = trace.stats.arrivals['P'][1]
    noise_start_t = p_arrival_time - config.pre_p_time
    noise_end_t = noise_start_t + config.noise_win_length

    trace_noise = trace.copy()
    # remove the mean...
    trace_noise.detrend(type='constant')
    # ...and the linear trend...
    trace_noise.detrend(type='linear')
    trace_noise.trim(starttime=noise_start_t, endtime=noise_end_t,
                     pad=True, fill_value=0)

    # S window for s/n ratio
    s_arrival_time = trace.stats.arrivals['S'][1]
    s_start_t = s_arrival_time - config.pre_s_time
    s_end_t = s_start_t + config.noise_win_length
    trace_cutS = trace.copy()
    # remove the mean...
    trace_cutS.detrend(type='constant')
    # ...and the linear trend...
    trace_cutS.detrend(type='linear')
    trace_cutS.trim(starttime=s_start_t, endtime=s_end_t,
                    pad=True, fill_value=0)

    rmsnoise2 = np.power(trace_noise.data, 2).sum()
    rmsnoise = np.sqrt(rmsnoise2)
    rmsS2 = np.power(trace_cutS.data, 2).sum()
    rmsS = np.sqrt(rmsS2)

    sn_ratio = rmsS/rmsnoise
    logging.info('%s %s: S/N: %.1f' % (trace.id, trace.stats.instrtype,
                                       sn_ratio))

    snratio_min = config.sn_min
    if sn_ratio < snratio_min:
        logging.warning('%s %s: S/N smaller than %g: skipping trace' %
                        (trace.id, trace.stats.instrtype, snratio_min))
        return True
    else:
        return False


def _process_trace(config, trace):
    # compute hypocentral distance
    if hypo_dist(trace) is None:
        logging.warning('%s: Unable to compute hypocentral distance: '
                        'skipping trace' % trace.id)
        return None

    # retrieve arrival times
    p_arrival_time = wave_arrival(trace, config.vp, 'P',
                                  config.p_arrival_tolerance)
    s_arrival_time = wave_arrival(trace, config.vs, 'S',
                                  config.s_arrival_tolerance)
    if (p_arrival_time is None or s_arrival_time is None):
        logging.warning('%s: Unable to get arrival times: '
                        'skipping trace' % trace.id)
        return None

    # copy trace for manipulation
    trace_process = trace.copy()

    comp = trace_process.stats.channel
    instrtype = trace_process.stats.instrtype
    if config.ignore_vertical and comp[-1] in ['Z', '1']:
        return None

    # check if the trace has (significant) signal
    if _has_not_signal(config, trace_process):
        return None

    # check if trace is clipped
    if _is_clipped(config, trace_process):
        return None

    # Remove instrument response
    if remove_instr_response(trace_process,
                             config.correct_instrumental_response,
                             config.pre_filt) is None:
        logging.warning('%s %s: Unable to remove instrument response: '
                        'skipping trace' % (trace_process.id, instrtype))
        return None

    try:
        filter_trace(config, trace_process)
    except ValueError:
        logging.warning('%s: Unknown instrument type: %s: '
                        'skipping trace' % (trace_process.id, instrtype))
        return None

    # Check if the trace has significant signal to noise ratio
    if _has_low_sn_ratio(config, trace_process):
        return None

    return trace_process


def process_traces(config, st):
    """Remove mean, deconvolve and ignore unwanted components."""
    # Since the count value is generally huge, we need to demean twice
    # to take into account for the rounding error
    st.detrend(type='constant')
    st.detrend(type='constant')
    # Merge stream to remove gaps and overlaps
    # Gaps are detected later as zero values
    st.merge(fill_value=0)
    out_st = Stream()
    for trace in st:
        trace_process = _process_trace(config, trace)
        if trace_process is not None:
            out_st.append(trace_process)

    if len(out_st) == 0:
        logging.error('No traces left! Exiting.')
        ssp_exit()

    return out_st
