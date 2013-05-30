# -*- coding: utf-8 -*-
# ssp_process_traces.py
#
# Trace processing for source_spec 
# (c) 2012 Claudio Satriano <satriano@ipgp.fr>
# (c) 2013 Claudio Satriano <satriano@ipgp.fr>,
#          Emanuela Matrullo <matrullo@geologie.ens.fr>
import logging
import numpy as np
from copy import deepcopy, copy
from obspy.core import Stream
from ssp_setup import dprint
from ssp_util import remove_instr_response, hypo_dist, wave_arrival

def process_traces(config, st):
    ''' Removes mean, deconvolves, and ignores unwanted components '''
    out_st = Stream()
    out_st_noise = Stream()
    for orig_trace in st.traces:
        # copy trace for manipulation
        trace = copy(orig_trace)
        trace.stats = deepcopy(orig_trace.stats)

        traceId = trace.id
        stats = trace.stats
        comp  = stats.channel
        instrtype = stats.instrtype
        if config.ignore_vertical and comp[-1] == 'Z':
            continue
        station = stats.station
        dprint('%s %s' % (station, comp))

        # compute hypocentral distance
        hd = hypo_dist(trace)
        if hd == None:
            logging.warning('%s %s: Unable to compute hypocentral distance: skipping trace' % (traceId, instrtype))
            continue
        trace.stats.hypo_dist = hd

        # check if the trace has (significant) signal
        # since the count value is generally huge, we need to demean twice
        # to take into account for the rounding error
        trace.detrend(type='constant')
        trace.detrend(type='constant')
        rms2 = np.power(trace.data, 2).sum()
        rms = np.sqrt(rms2)
        rms_min = config.rmsmin
        if rms <= rms_min:
            logging.warning('%s %s: Trace RMS smaller than %g: skipping trace' % (traceId, instrtype, rms_min))
            continue

        # Remove instrument response
        if remove_instr_response(trace, config.correct_sensitivity_only,
                                 config.pre_filt) == None:
            logging.warning('%s %s: Unable to remove instrument response: skipping trace' % (traceId, instrtype))
            continue

        # Check if the trace has (significant) signal to noise ratio
        #### start signal/noise ratio
        p_arrival_time = wave_arrival(trace, config.vp, 'P')
        s_arrival_time = wave_arrival(trace, config.vs, 'S')
        if (p_arrival_time == None or s_arrival_time == None):
            logging.warning('%s %s: Unable to get arrival times: skipping trace' % (traceId, instrtype))
            continue

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
        trace_noise.trim(starttime=noise_start_t, endtime=noise_end_t, pad=True, fill_value=0)
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
        trace_cutS.trim(starttime=s_start_t, endtime=s_end_t, pad=True, fill_value=0)
        #trace_cutS.plot()
        
        rmsnoise2 = np.power(trace_noise.data, 2).sum()
        rmsnoise = np.sqrt(rmsnoise2)
        rmsS2 = np.power(trace_cutS.data, 2).sum()
        rmsS = np.sqrt(rmsS2)

        sn_ratio = rmsS/rmsnoise
        logging.info('%s %s: S/N: %.1f' % (traceId, instrtype, sn_ratio))

        snratio_min = config.sn_min
        if sn_ratio < snratio_min:
            logging.warning('%s %s: S/N smaller than %g: skipping trace' % (traceId, instrtype, snratio_min))
            continue
        #### end signal/noise ratio

        # Noise time window for weighting function:
        noise_start_t = p_arrival_time - config.pre_p_time
        noise_end_t = noise_start_t + config.s_win_length
        trace_noise = copy(trace)
        trace_noise.stats = deepcopy(trace.stats)
        # remove the mean...
        trace_noise.detrend(type='constant')
        # ...and the linear trend...
        trace_noise.detrend(type='linear')
        trace_noise.trim(starttime=noise_start_t, endtime=noise_end_t, pad=True, fill_value=0)

        orig_trace.stats.hypo_dist = hd
        out_st.append(trace)
        out_st_noise.append(trace_noise)

    return out_st, out_st_noise
