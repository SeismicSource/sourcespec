# -*- coding: utf-8 -*-
# ssp_local_magnitude.py
#
# Local magnitude calculation for source_spec 
# (c) 2012 Claudio Satriano <satriano@ipgp.fr>
# (c) 2013 Claudio Satriano <satriano@ipgp.fr>,
#          Emanuela Matrullo <matrullo@geologie.ens.fr>
import logging
import numpy as np
from scipy.integrate import cumtrapz
from copy import deepcopy, copy
from obspy.signal import estimateMagnitude
from ssp_util import swave_arrival, cosine_taper

def local_magnitude(config, st, deconvolve=False):
    ''' Uses min/max amplitude for local magnitude estimation'''
    magnitudes = []
    for trace in st.traces:
        traceId = trace.getId()
        comp  = trace.stats.channel
        if comp[-1] == 'Z':
            continue

        # Skip traces which do not have hypo_dist defined
        try:
            trace.stats.hypo_dist
        except KeyError:
            continue

        # S time window
        s_arrival_time = swave_arrival(trace, config.vs)
        t1 = s_arrival_time - config.pre_s_time
        t2 = t1 + config.s_win_length

        trace_cut = copy(trace)
        trace_cut.stats = deepcopy(trace.stats)

        # Do a preliminary trim, in order to check if there is enough
        # data within the selected time window
        trace_cut.trim(starttime=t1, endtime=t2, pad=True, fill_value=0)
        npts = len(trace_cut.data)
        if npts == 0:
            logging.warning('%s: No data for the selected cut interval: skipping trace' % traceId)
            continue
        nzeros = len(np.where(trace_cut.data==0)[0])
        if nzeros > npts/4:
            logging.warning('%s: Too many gaps for the selected cut interval: skipping trace' % traceId)
            continue

        # If the check is ok, recover the full trace
        # (it will be cutted later on)
        trace_cut = copy(trace)
        trace_cut.stats = deepcopy(trace.stats)

        # remove the mean...
        trace_cut.detrend(type='constant')
        # ...and the linear trend...
        trace_cut.detrend(type='linear')
        # ...filter
        # TODO: parametrize?
        trace_cut.filter(type='bandpass',  freqmin=0.1, freqmax=20)

        instrtype = trace_cut.stats.instrtype
        if instrtype == 'acc':
            # integrate to velocity
            trace_cut.data = cumtrapz(trace_cut.data) * trace_cut.stats.delta
            trace_cut.stats.npts -= 1

        # trim...
        trace_cut.trim(starttime=t1, endtime=t2, pad=True, fill_value=0)
        # ...and taper
        cosine_taper(trace_cut.data, width=config.taper_halfwidth)

        delta_amp = trace_cut.data.max() - trace_cut.data.min() 
        delta_t = trace_cut.data.argmax() - trace_cut.data.argmin()    
        delta_t = delta_t / trace_cut.stats.sampling_rate

        # estimate local magnitude 
        if deconvolve:
            paz = trace_cut.stats.paz
        else:
            paz = {'poles': [],
                   'zeros': [],
                   'gain': 1.0, 'sensitivity': 1.0}

        ml = estimateMagnitude(paz, delta_amp, delta_t, trace_cut.stats.hypo_dist)

        magnitudes.append(ml)
        logging.info('%s %s: %s %.1f' % (traceId, trace.stats.instrtype, "Ml", ml))

    # average local magnitude
    Ml = np.mean(magnitudes)
    logging.info('Network Local Magnitude: %.2f' % Ml)
    return Ml
