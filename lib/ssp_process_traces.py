# -*- coding: utf-8 -*-
# ssp_process_traces.py
#
# Trace processing for source_spec 
# (c) 2012 Claudio Satriano <satriano@ipgp.fr>
# (c) 2013 Claudio Satriano <satriano@ipgp.fr>,
#          Emanuela Matrullo <matrullo@geologie.ens.fr>
import logging
import numpy as np
from obspy.core import Stream
from ssp_setup import dprint
from ssp_util import remove_instr_response, hypo_dist

def process_traces(config, st, skip_vertical=True):
    ''' Deconvolves and cuts traces, and removes unwanted components '''
    out_st = Stream()
    for trace in st.traces:
        traceId = trace.getId()
        stats = trace.stats
        comp  = stats.channel
        if skip_vertical and comp[-1] == 'Z':
            continue
        station = stats.station
        dprint('%s %s' % (station, comp))

        # check if the trace has (significant) signal
        # since the count value is generally huge, we need to demean twice
        # to take into account for the rounding error
        trace.detrend(type='constant')
        trace.detrend(type='constant')
        rms2 = np.power(trace.data, 2).sum()
        rms_min = 1e-10 #TODO: parametrize?
        if rms2 <= rms_min:
            logging.warning('%s: Trace RMS smaller than %g: skipping trace' % (traceId, rms_min))
            continue

        # Remove instrument response
        if remove_instr_response(trace, config.correct_sensitivity_only,
                                 config.pre_filt) == None:
            logging.warning('%s: Unable to remove instrument response: skipping trace' % traceId)
            continue

        # compute hypocentral distance
        hd = hypo_dist(trace)
        if hd == None:
            logging.warning('%s: Unable to compute hypocentral distance: skipping trace' % traceId)
            continue
        trace.stats.hypo_dist = hd

        out_st.append(trace)

    return out_st
