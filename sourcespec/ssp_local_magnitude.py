# -*- coding: utf-8 -*-
"""
Local magnitude calculation for sourcespec.

:copyright:
    2012 Claudio Satriano <satriano@ipgp.fr>

    2013-2014 Claudio Satriano <satriano@ipgp.fr>,
              Emanuela Matrullo <matrullo@geologie.ens.fr>

    2015-2017 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement, Version 2.1
    (http://www.cecill.info/index.en.html)
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import logging
import numpy as np
from scipy.integrate import cumtrapz
from copy import deepcopy, copy
from obspy.signal.filter import envelope
from obspy.signal.invsim import estimate_magnitude
from obspy.signal.util import smooth
from obspy.signal.trigger import trigger_onset
from sourcespec.ssp_util import cosine_taper


def local_magnitude(config, st, deconvolve=False):
    """Compute local magnitude from min/max amplitude."""
    magnitudes = []
    for trace in st.traces:
        traceId = trace.get_id()
        comp = trace.stats.channel
        if comp[-1] in ['Z', 'H']:
            continue

        # Skip traces which do not have hypo_dist defined
        try:
            trace.stats.hypo_dist
        except (KeyError, AttributeError):
            continue

        trace_env = copy(trace)
        trace_env.stats = deepcopy(trace.stats)

        # remove the mean...
        trace_env.detrend(type='constant')
        # ...and the linear trend...
        trace_env.detrend(type='linear')
        # ...filter
        # TODO: parametrize?
        freqmin = 0.1
        freqmax = 20.
        nyquist = 1./(2. * trace.stats.delta)
        if freqmax >= nyquist:
            freqmax = nyquist * 0.999
            msg = '%s: maximum frequency for bandpass filtering ' % trace.id
            msg += 'in local magnitude computation is larger than or equal '
            msg += 'to Nyquist. Setting it to %s Hz' % freqmax
            logging.warning(msg)
        trace.filter(type='bandpass', freqmin=freqmin, freqmax=freqmax)
        trace_env.data = envelope(trace_env.data)
        trace_env.data = smooth(trace_env.data, 100)

        # Uncomment these lines to see the effect of envelope and smoothing
        #plt.figure()
        #plt.plot(trace_env, color='gray')
        #plt.grid(True)
        #plt.show()

        trace_noise = copy(trace_env)
        trace_noise.stats = deepcopy(trace_env.stats)
        trace_signal = copy(trace_env)
        trace_signal.stats = deepcopy(trace_env.stats)

        # Skip traces which do not have arrivals
        try:
            p_arrival_time = trace.stats.arrivals['P'][1]
        except:
            continue
        t1 = p_arrival_time - config.pre_mag_time
        t2 = p_arrival_time + config.pre_mag_time

        trace_noise.trim(starttime=t1, endtime=p_arrival_time,
                         pad=True, fill_value=0)
        trace_signal.trim(starttime=p_arrival_time, endtime=t2,
                          pad=True, fill_value=0)

        ampmin = trace_noise.data.mean()
        ampmax = trace_signal.data.mean()

        if ampmax <= ampmin:
            continue

        trigger = trigger_onset(trace_env.data, ampmax, ampmin,
                                max_len=9e99, max_len_delete=False)[0]

        df = trace.stats.sampling_rate
        triggeron = trigger[0]/df
        #triggeroff = trigger[1]/df  #Not used --C.S.
        #start_trace = trace.stats['starttime']
        #t_end = start_trace + triggeroff
        t1_end = p_arrival_time + triggeron
        if t1_end <= p_arrival_time:
            t_end = p_arrival_time + config.mag_win_length
        else:
            t_end = p_arrival_time + triggeron

        trace_cut = copy(trace)
        trace_cut.stats = deepcopy(trace.stats)

        # Do a preliminary trim, in order to check if there is enough
        # data within the selected time window
        #trace_cut.trim(starttime=t1, endtime=t2, pad=True, fill_value=0)
        trace_cut.trim(starttime=p_arrival_time, endtime=t_end,
                       pad=True, fill_value=0)
        #plt.figure()
        #plt.plot(trace_cut, color='gray')
        #plt.grid(True)
        #plt.show()

        npts = len(trace_cut.data)
        if npts == 0:
            logging.warning('%s: No data for the selected cut interval: '
                            'skipping trace' % traceId)
            continue
        nzeros = len(np.where(trace_cut.data == 0)[0])
        if nzeros > npts/4:
            logging.warning('%s: Too many gaps for the selected cut '
                            'interval: skipping trace' % traceId)
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
        #trace_cut.trim(starttime=t1, endtime=t2, pad=True, fill_value=0)
        trace_cut.trim(starttime=p_arrival_time, endtime=t_end,
                       pad=True, fill_value=0)

        # ...and taper
        cosine_taper(trace_cut.data, width=config.taper_halfwidth)

        delta_amp = trace_cut.data.max() - trace_cut.data.min()
        delta_t = trace_cut.data.argmax() - trace_cut.data.argmin()
        delta_t = delta_t / trace_cut.stats.sampling_rate

        # estimate local magnitude
        if deconvolve:
            paz = trace_cut.stats.paz
            if paz is None:
                logging.warning('%s: no poles and zeros for trace: '
                                'skipping trace' % traceId)
                continue
        else:
            paz = {'poles': [],
                   'zeros': [],
                   'gain': 1.0, 'sensitivity': 1.0}

        ml = estimate_magnitude(paz, delta_amp, delta_t,
                                trace_cut.stats.hypo_dist)
        magnitudes.append(ml)
        logging.info('%s %s: %s %.1f' %
                     (traceId, trace.stats.instrtype, "Ml", ml))

    # average local magnitude
    Ml = np.mean(magnitudes)
    logging.info('Network Local Magnitude: %.2f' % Ml)
    return Ml
