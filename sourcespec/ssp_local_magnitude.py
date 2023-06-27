# -*- coding: utf-8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Local magnitude calculation for sourcespec.

:copyright:
    2012 Claudio Satriano <satriano@ipgp.fr>

    2013-2014 Claudio Satriano <satriano@ipgp.fr>,
              Emanuela Matrullo <matrullo@geologie.ens.fr>

    2015-2023 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
import contextlib
import logging
import numpy as np
from copy import deepcopy
from obspy.signal.filter import envelope
from obspy.signal.invsim import WOODANDERSON
from obspy.signal.util import smooth
from obspy.signal.trigger import trigger_onset
from sourcespec.ssp_data_types import SpectralParameter
from sourcespec.ssp_util import cosine_taper
from sourcespec.ssp_util import remove_instr_response
logger = logging.getLogger(__name__.split('.')[-1])


def _check_nyquist(freqmax, trace):
    """Check if freqmax is smaller than Nyquist frequency."""
    nyquist = 1. / (2. * trace.stats.delta)
    if freqmax >= nyquist:
        freqmax = nyquist * 0.999
        msg = (
            f'{trace.id}: maximum frequency for bandpass filtering in '
            'local magnitude computation is larger than or equal '
            f'to Nyquist. Setting it to {freqmax} Hz'
        )
        if msg not in _check_nyquist.messages:
            _check_nyquist.messages.append(msg)
            logger.warning(msg)
    return freqmax
_check_nyquist.messages = []  # noqa


def _get_cut_times(config, tr):
    """Get trace cut times between P arrival and end of envelope coda."""
    tr_env = tr.copy()
    # remove the mean...
    tr_env.detrend(type='constant')
    # ...and the linear trend...
    tr_env.detrend(type='linear')
    # ...filter
    freqmin = 1.
    freqmax = _check_nyquist(freqmax=20, trace=tr)
    cosine_taper(tr_env.data, width=config.taper_halfwidth)
    tr_env.filter(type='bandpass', freqmin=freqmin, freqmax=freqmax)
    tr_env.data = envelope(tr_env.data)
    tr_env.data = smooth(tr_env.data, 100)

    # Skip traces which do not have arrivals
    try:
        p_arrival_time = tr.stats.arrivals['P'][1]
    except Exception as e:
        raise RuntimeError(
            f'{tr.id}: Trace has no P arrival: skipping trace'
        )from e
    t1 = p_arrival_time - config.win_length
    t2 = p_arrival_time + config.win_length

    tr_noise = tr_env.copy()
    tr_signal = tr_env.copy()
    tr_noise.trim(starttime=t1, endtime=p_arrival_time,
                  pad=True, fill_value=0)
    tr_signal.trim(starttime=p_arrival_time, endtime=t2,
                   pad=True, fill_value=0)
    ampmin = tr_noise.data.mean()
    ampmax = tr_signal.data.mean()
    if ampmax <= ampmin:
        raise RuntimeError(
            f'{tr.id}: Trace has too high noise before P arrival: '
            'skipping trace'
        )

    trigger = trigger_onset(
        tr_env.data, ampmax, ampmin, max_len=9e99, max_len_delete=False)[0]
    t0 = p_arrival_time
    t1 = t0 + trigger[-1] * tr.stats.delta
    t1 = min(t1, tr.stats.endtime)
    return t0, t1


def _process_trace(config, tr, t0, t1):
    """Convert to Wood-Anderson, filter, trim."""
    tr_process = tr.copy()
    # Do a preliminary trim, in order to check if there is enough
    # data within the selected time window
    tr_process.trim(starttime=t0, endtime=t1, pad=True, fill_value=0)
    npts = len(tr_process.data)
    if npts == 0:
        raise RuntimeError(
            f'{tr.id}: No data for the selected cut interval: skipping trace'
        )
    nzeros = len(np.where(tr_process.data == 0)[0])
    if nzeros > npts / 4:
        raise RuntimeError(
            f'{tr.id}: Too many gaps for the selected cut interval: '
            'skipping trace'
        )

    # If the check is ok, recover the full trace
    # (it will be cut later on)
    tr_process = tr.copy()
    # remove the mean...
    tr_process.detrend(type='constant')
    # ...and the linear trend...
    tr_process.detrend(type='linear')
    freqmin = config.ml_bp_freqmin
    freqmax = _check_nyquist(freqmax=config.ml_bp_freqmax, trace=tr)
    # ...remove response...
    pre_filt = (freqmin, freqmin * 1.1, freqmax * 0.9, freqmax)
    if config.correct_instrumental_response:
        remove_instr_response(tr_process, pre_filt)
    # ...filter
    tr_process.filter(type='bandpass', freqmin=freqmin, freqmax=freqmax)
    # Convert to Wood-Anderson
    # note: conversion to Wood-Anderson integrates the signal once
    if tr.stats.instrtype == 'acc':
        WA_double_int = deepcopy(WOODANDERSON)
        # we remove a zero to add an integration step
        WA_double_int['zeros'].pop()
        tr_process.simulate(paz_simulate=WA_double_int)
    else:
        tr_process.simulate(paz_simulate=WOODANDERSON)
    # trim...
    tr_process.trim(starttime=t0, endtime=t1, pad=True, fill_value=0)
    # ...and taper
    cosine_taper(tr_process.data, width=config.taper_halfwidth)
    return tr_process


def _compute_local_magnitude(config, amp, h_dist):
    """Compute local magnitude using the Richter formula."""
    a = config.a
    b = config.b
    c = config.c
    return (
        np.log10(amp) + a * np.log10(h_dist / 100.0) + b * (h_dist - 100.0) + c
    )


def local_magnitude(config, st, proc_st, sspec_output):
    """Compute local magnitude from max absolute W-A amplitude."""
    logger.info('Computing local magnitude...')
    # We only use traces selected for proc_st
    trace_ids = {tr.id for tr in proc_st}
    station_parameters = sspec_output.station_parameters
    for tr_id in sorted(trace_ids):
        tr = st.select(id=tr_id)[0]

        # only analyze traces for which we have the other
        # source parameters defined
        key = f'{tr_id[:-1]}H'
        try:
            station_pars = station_parameters[key]
        except KeyError:
            continue

        # Make sure that the param_Ml object is always defined, even when Ml
        # is not computed (i.e., "continue" below)
        try:
            param_Ml = station_pars.Ml
        except KeyError:
            param_Ml = SpectralParameter(
                id='Ml', value=np.nan, format='{:.2f}')
            station_pars.Ml = param_Ml

        # only compute Ml for horizontal components
        comp = tr.stats.channel
        if comp[-1] in ['Z', 'H']:
            continue

        try:
            t0, t1 = _get_cut_times(config, tr)
            tr_process = _process_trace(config, tr, t0, t1)
        except RuntimeError as msg:
            logger.warning(msg)
            continue

        # Local magnitude
        # amp must be in millimeters for local magnitude computation
        amp = np.abs(tr_process.max()) * 1e3
        h_dist = tr_process.stats.hypo_dist
        ml = _compute_local_magnitude(config, amp, h_dist)
        statId = f'{tr_id} {tr.stats.instrtype}'
        logger.info(f'{statId}: Ml {ml:.1f}')
        # compute average with the other component, if available
        with contextlib.suppress(KeyError):
            old_ml = param_Ml.value
            if not np.isnan(old_ml):
                ml = 0.5 * (ml + old_ml)
        param_Ml.value = ml
