# -*- coding: utf-8 -*-
# ssp_local_magnitude.py
#
# Local magnitude calculation for source_spec 
# (c) 2012 Claudio Satriano <satriano@ipgp.fr>
# (c) 2013 Claudio Satriano <satriano@ipgp.fr>,
#          Emanuela Matrullo <matrullo@geologie.ens.fr>
import numpy as np
from copy import deepcopy, copy
from obspy.signal import estimateMagnitude
from ssp_util import cosine_taper

def local_magnitude(st, deconvolve=False):
    magnitudes = []
    for trace in st.traces:
        traceId = trace.getId()
        # use min/max amplitude for local magnitude estimation
        comp  = trace.stats.channel
        if comp[-1] == 'Z':
            continue

        # copy trace for local magnitude estimation
        trace_cp = copy(trace)
        trace_cp.stats = deepcopy(trace.stats)

        # remove the mean...
        trace_cp.detrend(type='constant')
        # ...and the linear trend...
        trace_cp.detrend(type='linear')
        # ...filter
        trace_cp.filter(type='bandpass',  freqmin=0.1, freqmax=20)
        #trace_cp.filter(type='bandpass',  freqmin=0.5, freqmax=20)
        # ...and taper
        cosine_taper(trace_cp.data, width=0.5)

        delta_amp = trace_cp.data.max() - trace_cp.data.min() 
        #delta_amp = max(abs(trace_cp.data))
        delta_t = trace_cp.data.argmax() - trace_cp.data.argmin()    
        delta_t = delta_t / trace_cp.stats.sampling_rate

        #estimate Magnitude 
        instrtype = trace_cp.stats.instrtype
        if instrtype == 'acc':
            poles = [0.0j] #integration to velocity 
        else:
            poles = []
        paz = {'poles': poles,
               'zeros': [],
               'gain': 1.0, 'sensitivity': 1.0}
        if deconvolve:
            paz = trace_cp.stats.paz
            try:
                trace_cp.stats.hypo_dist
            except KeyError:
                continue
        ml = estimateMagnitude(paz, delta_amp, delta_t, trace_cp.stats.hypo_dist)

        magnitudes.append(ml)
        print '%s %s %s %.1f' % (traceId, trace.stats.instrtype, "Ml", ml)
    # end of loop on stations for building spectra 

    # average local magnitude
    Ml = np.mean(magnitudes)
    print '\nNetwork Local Magnitude: %.2f' % Ml
    return Ml
