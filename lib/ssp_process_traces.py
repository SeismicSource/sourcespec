# -*- coding: utf-8 -*-
# ssp_process_traces.py
#
# (c) 2012 Claudio Satriano <satriano@ipgp.fr>
# (c) 2013-2014 Claudio Satriano <satriano@ipgp.fr>,
#               Emanuela Matrullo <matrullo@geologie.ens.fr>
# (c) 2015-2016 Claudio Satriano <satriano@ipgp.fr>
"""Trace processing for source_spec."""
from __future__ import division
import logging
import numpy as np
from copy import deepcopy, copy
from obspy.core import Stream
from ssp_setup import dprint, ssp_exit
from ssp_util import remove_instr_response, hypo_dist, wave_arrival


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


def process_traces(config, st):
    """Remove mean, deconvolve and ignore unwanted components."""
    out_st = Stream()
    for orig_trace in st.traces:
        # compute hypocentral distance
        hd = hypo_dist(orig_trace)
        if hd is None:
            logging.warning('%s: Unable to compute hypocentral distance: '
                            'skipping trace' % orig_trace.id)
            continue
        orig_trace.stats.hypo_dist = hd

        # retrieve arrival times
        p_arrival_time = wave_arrival(orig_trace, config.vp, 'P')
        s_arrival_time = wave_arrival(orig_trace, config.vs, 'S')
        if (p_arrival_time is None or s_arrival_time is None):
            logging.warning('%s: Unable to get arrival times: '
                            'skipping trace' % orig_trace.id)
            continue

        # copy trace for manipulation
        trace = copy(orig_trace)
        trace.stats = deepcopy(orig_trace.stats)

        stats = trace.stats
        comp = stats.channel
        instrtype = stats.instrtype
        if config.ignore_vertical and comp[-1] in ['Z', '1']:
            continue
        station = stats.station
        dprint('%s %s' % (station, comp))

        # check if the trace has (significant) signal
        # since the count value is generally huge, we need to demean twice
        # to take into account for the rounding error
        trace.detrend(type='constant')
        trace.detrend(type='constant')
        rms2 = np.power(trace.data, 2).sum()
        rms = np.sqrt(rms2)
        rms_min = config.rmsmin
        if rms <= rms_min:
            logging.warning('%s %s: Trace RMS smaller than %g: '
                            'skipping trace' % (trace.id, instrtype, rms_min))
            continue

        # check if trace is clipped
        clip_tollerance = config.clip_tollerance
        clip_max = (1 - clip_tollerance/100.) * trace.data.max()
        clip_min = (1 - clip_tollerance/100.) * trace.data.min()
        nclips = (trace.data > clip_max).sum()
        nclips += (trace.data < clip_min).sum()
        clip_nmax = config.clip_nmax
        if float(nclips)/trace.stats.npts > clip_nmax/100.:
            logging.warning('%s %s: Trace is clipped for more than %.2f%% '
                            'with %.2f%% tollerance: skipping trace' %
                            (trace.id, instrtype, clip_nmax, clip_tollerance))
            continue

        # Remove instrument response
        if remove_instr_response(trace, config.correct_instrumental_response,
                                 config.pre_filt) is None:
            logging.warning('%s %s: Unable to remove instrument response: '
                            'skipping trace' % (trace.id, instrtype))
            continue

        try:
            filter_trace(config, trace)
        except ValueError:
            logging.warning('%s: Unknown instrument type: %s: '
                            'skipping trace' % (trace.id, instrtype))
            continue

        # Check if the trace has (significant) signal to noise ratio
        # ***start signal/noise ratio

        # noise time window for s/n ratio
        #noise_start_t = trace.stats['starttime']
        #noise_end_t = p_arrival_time
        #noise_win_length = noise_end_t - noise_start_t
        noise_start_t = p_arrival_time - config.pre_p_time
        noise_end_t = noise_start_t + config.noise_win_length

        trace_noise = copy(trace)
        trace_noise.stats = deepcopy(trace.stats)
        # remove the mean...
        trace_noise.detrend(type='constant')
        # ...and the linear trend...
        trace_noise.detrend(type='linear')
        trace_noise.trim(starttime=noise_start_t, endtime=noise_end_t,
                         pad=True, fill_value=0)
        #trace_noise.plot()

        # S window for s/n ratio
        s_start_t = s_arrival_time - config.pre_s_time
        #print traceId, s_arrival_time, config.pre_s_time
        #s_end_t = s_start_t + noise_win_length
        s_end_t = s_start_t + config.noise_win_length
        trace_cutS = copy(trace)
        trace_cutS.stats = deepcopy(trace.stats)
        # remove the mean...
        trace_cutS.detrend(type='constant')
        # ...and the linear trend...
        trace_cutS.detrend(type='linear')
        trace_cutS.trim(starttime=s_start_t, endtime=s_end_t,
                        pad=True, fill_value=0)
        #trace_cutS.plot()

        rmsnoise2 = np.power(trace_noise.data, 2).sum()
        rmsnoise = np.sqrt(rmsnoise2)
        rmsS2 = np.power(trace_cutS.data, 2).sum()
        rmsS = np.sqrt(rmsS2)

        sn_ratio = rmsS/rmsnoise
        logging.info('%s %s: S/N: %.1f' % (trace.id, instrtype, sn_ratio))

        snratio_min = config.sn_min
        if sn_ratio < snratio_min:
            logging.warning('%s %s: S/N smaller than %g: '
                            'skipping trace' %
                            (trace.id, instrtype, snratio_min))
            continue
        # ***end signal/noise ratio

        orig_trace.stats.hypo_dist = hd
        out_st.append(trace)

    if len(out_st) == 0:
        logging.error('No traces left! Exiting.')
        ssp_exit()

    return out_st
