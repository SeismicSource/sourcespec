# -*- coding: utf-8 -*-
"""
Arrival time calculation for sourcespec.

:copyright:
    2012-2021 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement, Version 2.1
    (http://www.cecill.info/index.en.html)
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import os
from glob import glob
import logging
import warnings
from sourcespec.ssp_setup import ssp_exit
from obspy.taup import TauPyModel
model = TauPyModel(model='iasp91')
logger = logging.getLogger(__name__.split('.')[-1])


def _wave_arrival_nll(trace, phase, NLL_time_dir):
    """Arrival time using a NLL grid."""
    if trace.stats.hypo.origin_time is None:
        return
    if NLL_time_dir is None:
        return
    try:
        from nllgrid import NLLGrid
    except ImportError:
        logger.error('Error: the "nllgrid" python module is required '
                     'for "NLL_time_dir".')
        ssp_exit()
    grdfile = '*.{}.{}.time.hdr'.format(phase, trace.stats.station)
    grdfile = os.path.join(NLL_time_dir, grdfile)
    try:
        grdfile = glob(grdfile)[0]
    except IndexError:
        return
    grd = NLLGrid(grdfile)
    if grd.type != 'TIME':
        return
    hypo_x, hypo_y = grd.project(trace.stats.hypo.longitude,
                                 trace.stats.hypo.latitude)
    hypo_z = trace.stats.hypo.depth
    tt = grd.get_value(hypo_x, hypo_y, hypo_z)
    return trace.stats.hypo.origin_time + tt


def _wave_arrival_vel(trace, phase, vel):
    """Arrival time using a constant velocity (in km/s)."""
    if trace.stats.hypo.origin_time is None:
        return
    if vel is None:
        return
    tt = trace.stats.hypo_dist / vel
    return trace.stats.hypo.origin_time + tt


def _wave_arrival_taup(trace, phase):
    """Arrival time using taup."""
    if trace.stats.hypo.origin_time is None:
        return
    phase_list = [phase.lower(), phase]
    with warnings.catch_warnings(record=True) as warns:
        try:
            arrivals = model.get_travel_times(
                        source_depth_in_km=trace.stats.hypo.depth,
                        distance_in_degree=trace.stats.gcarc,
                        phase_list=phase_list)
        except Exception:
            trace.stats.hypo.depth = 0.
            arrivals = model.get_travel_times(
                        source_depth_in_km=trace.stats.hypo.depth,
                        distance_in_degree=trace.stats.gcarc,
                        phase_list=phase_list)
        for w in warns:
            message = str(w.message)
            # Ignore a specific obspy.taup warning we do not care about
            if '#2280' in message:
                continue
            logger.warning(message)
    tt = min(a.time for a in arrivals)
    return trace.stats.hypo.origin_time + tt


def _validate_pick(pick, theo_pick_time, tolerance, trace_id):
    if theo_pick_time is None:
        return True
    delta_t = pick.time - theo_pick_time
    if abs(delta_t) > tolerance:  # seconds
        msg = '%s: measured %s pick time - theoretical time = %.1f s.' %\
              (trace_id, pick.phase, delta_t)
        logger.warning(msg)
        return False
    return True


def wave_arrival(trace, phase, tolerance=4., vel=None, NLL_time_dir=None):
    """
    Obtain arrival time for a given phase and a given trace.

    Returns the theoretical arrival time if no pick is available
    or if the pick is too different from the theoretical arrival.
    """
    try:
        arrival = trace.stats.arrivals[phase][1]
        return arrival
    except KeyError:
        pass
    theo_pick_time = _wave_arrival_nll(trace, phase, NLL_time_dir)
    if theo_pick_time is None:
        theo_pick_time = _wave_arrival_vel(trace, phase, vel)
    if theo_pick_time is None:
        theo_pick_time = _wave_arrival_taup(trace, phase)
    for pick in (p for p in trace.stats.picks if p.phase == phase):
        if _validate_pick(pick, theo_pick_time, tolerance, trace.id):
            logger.info('%s: found %s pick' % (trace.id, phase))
            trace.stats.arrivals[phase] = (phase, pick.time)
            return pick.time
    logger.info('%s: using theoretical %s pick' % (trace.id, phase))
    trace.stats.arrivals[phase] = (phase + 'theo', theo_pick_time)
    return theo_pick_time
