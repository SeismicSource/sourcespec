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
from math import asin, degrees
from obspy.taup import TauPyModel
model = TauPyModel(model='iasp91')
logger = logging.getLogger(__name__.split('.')[-1])


def _get_nll_grd(phase, station, type, NLL_time_dir):
    # Lazy-import here, since nllgrid is not an installation requirement
    from nllgrid import NLLGrid
    for _station in station, 'DEFAULT':
        key = '{}_{}_{}'.format(phase, _station, type)
        try:
            # first try to lookup in cache
            grd = _get_nll_grd.grds[key]
            return grd
        except KeyError:
            pass
        try:
            grdfile = '*.{}.{}.{}.hdr'.format(phase, _station, type)
            grdfile = os.path.join(NLL_time_dir, grdfile)
            grdfile = glob(grdfile)[0]
            grd = NLLGrid(grdfile)
            # cache NLL grid
            _get_nll_grd.grds[key] = grd
            return grd
        except IndexError:
            # IndexError from glob()[0]
            pass
    raise RuntimeError


# dictionary to cache NLL grids
_get_nll_grd.grds = dict()


def _wave_arrival_nll(trace, phase, NLL_time_dir, focmec):
    """Travel time and takeoff angle using a NLL grid."""
    if NLL_time_dir is None:
        raise RuntimeError
    station = trace.stats.station
    travel_time = takeoff_angle = None
    grdtypes = ['time']
    if focmec:
        grdtypes.append('angle')
    for type in grdtypes:
        try:
            grd = _get_nll_grd(phase, station, type, NLL_time_dir)
        except RuntimeError:
            logger.warning(
                '{}: Cannot find NLL {} grid. '
                'Falling back to another method'.format(trace.id, type))
            raise RuntimeError
        if grd.station == 'DEFAULT':
            sta_x, sta_y = grd.project(
                trace.stats.coords.longitude, trace.stats.coords.latitude)
            grd.sta_x, grd.sta_y = sta_x, sta_y
        hypo_x, hypo_y = grd.project(
            trace.stats.hypo.longitude, trace.stats.hypo.latitude)
        hypo_z = trace.stats.hypo.depth
        if type == 'time':
            travel_time = grd.get_value(hypo_x, hypo_y, hypo_z)
        elif type == 'angle':
            azimuth, takeoff_angle, quality = grd.get_value(
                hypo_x, hypo_y, hypo_z)
    return travel_time, takeoff_angle


def _wave_arrival_vel(trace, phase, vel):
    """Travel time and takeoff angle using a constant velocity (in km/s)."""
    if vel is None:
        raise RuntimeError
    travel_time = trace.stats.hypo_dist / vel
    takeoff_angle = degrees(asin(trace.stats.epi_dist/trace.stats.hypo_dist))
    # takeoff angle is 180° upwards and 0° downwards
    takeoff_angle = 180. - takeoff_angle
    return travel_time, takeoff_angle


def _wave_arrival_taup(trace, phase):
    """Travel time and takeoff angle using taup."""
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
    # get first arrival
    first_arrival = sorted(arrivals, key=lambda a: a.time)[0]
    travel_time = first_arrival.time
    takeoff_angle = first_arrival.takeoff_angle
    return travel_time, takeoff_angle


def _validate_pick(pick, theo_pick_time, tolerance, trace_id):
    if theo_pick_time is None:
        return True
    delta_t = pick.time - theo_pick_time
    if abs(delta_t) > tolerance:  # seconds
        logger.warning(
            '{}: measured {} pick time - theoretical time = {:.1f} s'.format(
                trace_id, pick.phase, delta_t
            ))
        return False
    return True


def add_arrivals_to_trace(trace, config):
    """
    Add P and S arrival times and takeoff angles to trace.

    Uses the theoretical arrival time if no pick is available
    or if the pick is too different from the theoretical arrival.
    """
    tolerance = config.p_arrival_tolerance
    NLL_time_dir = config.NLL_time_dir
    vel = {'P': config.vp_tt, 'S': config.vs_tt}
    focmec = config.rps_from_focal_mechanism
    for phase in 'P', 'S':
        try:
            travel_time, takeoff_angle = _wave_arrival_nll(
                trace, phase, NLL_time_dir, focmec)
            method = 'NonLinLoc grid'
        except RuntimeError:
            try:
                travel_time, takeoff_angle = _wave_arrival_vel(
                    trace, phase, vel[phase])
                method = 'constant V{}: {:.1f} km/s'.format(
                    phase.lower(), vel[phase])
            except RuntimeError:
                try:
                    travel_time, takeoff_angle = _wave_arrival_taup(
                        trace, phase)
                    method = 'global velocity model (iasp91)'
                except RuntimeError:
                    return
        if trace.stats.hypo.origin_time is None:
            theo_pick_time = None
        else:
            theo_pick_time = trace.stats.hypo.origin_time + travel_time
        found_pick = False
        for pick in (p for p in trace.stats.picks if p.phase == phase):
            if _validate_pick(pick, theo_pick_time, tolerance, trace.id):
                logger.info('{}: found {} pick'.format(trace.id, phase))
                trace.stats.arrivals[phase] = (phase, pick.time)
                found_pick = True
                break
        if not found_pick and theo_pick_time is not None:
            logger.info('{}: using theoretical {} pick from {}'.format(
                trace.id, phase, method))
            trace.stats.arrivals[phase] = (phase + 'theo', theo_pick_time)
        trace.stats.takeoff_angles[phase] = takeoff_angle
        if focmec:
            logger.info(
                '{}: {} takeoff angle: {:.1f} computed from {}'.format(
                    trace.id, phase, takeoff_angle, method
                ))
