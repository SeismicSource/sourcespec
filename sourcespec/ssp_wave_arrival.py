# -*- coding: utf-8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Arrival time calculation for sourcespec.

:copyright:
    2012-2022 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
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


def _wave_arrival_vel(trace, vel):
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
    kwargs = dict(
        source_depth_in_km=trace.stats.hypo.depth,
        distance_in_degree=trace.stats.gcarc,
        phase_list=phase_list)
    with warnings.catch_warnings(record=True) as warns:
        try:
            arrivals = model.get_travel_times(**kwargs)
        except Exception:
            trace.stats.hypo.depth = 0.
            kwargs['source_depth_in_km'] = 0.
            arrivals = model.get_travel_times(**kwargs)
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


def _wave_arrival(trace, phase, config):
    """Get travel time and takeoff angle."""
    NLL_time_dir = config.NLL_time_dir
    focmec = config.rp_from_focal_mechanism
    vel = {'P': config.vp_tt, 'S': config.vs_tt}
    try:
        travel_time, takeoff_angle =\
            _wave_arrival_nll(trace, phase, NLL_time_dir, focmec)
        method = 'NonLinLoc grid'
        return travel_time, takeoff_angle, method
    except RuntimeError:
        pass
    try:
        travel_time, takeoff_angle =\
            _wave_arrival_vel(trace, vel[phase])
        method = 'constant V{}: {:.1f} km/s'.format(phase.lower(), vel[phase])
        return travel_time, takeoff_angle, method
    except RuntimeError:
        pass
    try:
        travel_time, takeoff_angle = _wave_arrival_taup(trace, phase)
        method = 'global velocity model (iasp91)'
        return travel_time, takeoff_angle, method
    except RuntimeError:
        raise


def _validate_pick(pick, theo_pick_time, tolerance, trace_id):
    """Check if a pick is valid, i.e., close enough to theoretical one."""
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


def _get_theo_pick_time(trace, travel_time):
    try:
        theo_pick_time = trace.stats.hypo.origin_time + travel_time
    except TypeError:
        theo_pick_time = None
    return theo_pick_time


def _travel_time_from_pick(trace, pick_time):
    try:
        travel_time = pick_time - trace.stats.hypo.origin_time
    except TypeError:
        travel_time = None
    return travel_time


def _find_picks(trace, phase, theo_pick_time, tolerance):
    """Search for valid picks in trace stats. Return pick time if found."""
    for pick in (p for p in trace.stats.picks if p.phase == phase):
        if _validate_pick(pick, theo_pick_time, tolerance, trace.id):
            trace.stats.arrivals[phase] = (phase, pick.time)
            return pick.time
    return None


def add_arrivals_to_trace(trace, config):
    """
    Add P and S arrival times and takeoff angles to trace.

    Uses the theoretical arrival time if no pick is available
    or if the pick is too different from the theoretical arrival.
    """
    tolerance = config.p_arrival_tolerance
    for phase in 'P', 'S':
        key = '{}_{}'.format(trace.id, phase)
        # First, see if there are cached values
        try:
            trace.stats.arrivals[phase] =\
                add_arrivals_to_trace.pick_cache[key]
            trace.stats.travel_times[phase] =\
                add_arrivals_to_trace.travel_time_cache[key]
            trace.stats.takeoff_angles[phase] =\
                add_arrivals_to_trace.angle_cache[key]
            continue
        except KeyError:
            pass
        # If no cache is available, compute travel_time and takeoff_angle
        try:
            travel_time, takeoff_angle, method =\
                _wave_arrival(trace, phase, config)
            theo_pick_time = _get_theo_pick_time(trace, travel_time)
            pick_time = _find_picks(trace, phase, theo_pick_time, tolerance)
        except RuntimeError:
            continue
        if pick_time is not None:
            logger.info('{}: found {} pick'.format(trace.id, phase))
            travel_time = \
                _travel_time_from_pick(trace, pick_time) or travel_time
            pick_phase = phase
        elif theo_pick_time is not None:
            logger.info('{}: using theoretical {} pick from {}'.format(
                trace.id, phase, method))
            pick_time = theo_pick_time
            pick_phase = phase + 'theo'
        else:
            continue
        if config.rp_from_focal_mechanism:
            logger.info(
                '{}: {} takeoff angle: {:.1f} computed from {}'.format(
                    trace.id, phase, takeoff_angle, method
                ))
        add_arrivals_to_trace.pick_cache[key] =\
            trace.stats.arrivals[phase] = (pick_phase, pick_time)
        add_arrivals_to_trace.travel_time_cache[key] =\
            trace.stats.travel_times[phase] = travel_time
        add_arrivals_to_trace.angle_cache[key] =\
            trace.stats.takeoff_angles[phase] = takeoff_angle


add_arrivals_to_trace.pick_cache = dict()
add_arrivals_to_trace.travel_time_cache = dict()
add_arrivals_to_trace.angle_cache = dict()
