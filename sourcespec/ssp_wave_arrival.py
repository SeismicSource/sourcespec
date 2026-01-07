# -*- coding: utf-8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Arrival time calculation for sourcespec.

:copyright:
    2012-2025 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
import os
import contextlib
from glob import glob
import logging
import warnings
from math import asin, degrees
from obspy.taup import TauPyModel
from .ssp_velocity_model import CrustalVelocityModel
model = TauPyModel(model='iasp91')
logger = logging.getLogger(__name__.rsplit('.', maxsplit=1)[-1])


def _get_nll_grd(phase, station, grd_type, NLL_time_dir):
    """
    Locate and cache an NLLGrid file for a given phase/station/grid type.

    :param phase: Phase name
    :type phase: str
    :param station: Station code or 'DEFAULT'
    :type station: str
    :param grd_type: Grid type (e.g. 'time' or 'angle')
    :type grd_type: str
    :param NLL_time_dir: Directory containing NLL grid files
    :type NLL_time_dir: str

    :return: NLLGrid instance for the requested grid
    :rtype: nllgrid.NLLGrid

    :raises RuntimeError: if no matching grid file can be found
    """
    # Lazy-import here, since nllgrid is not an installation requirement
    # pylint: disable=import-outside-toplevel
    from nllgrid import NLLGrid
    for _station in station, 'DEFAULT':
        key = f'{phase}_{_station}_{grd_type}'
        with contextlib.suppress(KeyError):
            return _get_nll_grd.grds[key]
        with contextlib.suppress(IndexError):
            grdfile = f'*.{phase}.{_station}.{grd_type}.hdr'
            grdfile = os.path.join(NLL_time_dir, grdfile)
            grdfile = glob(grdfile)[0]
            grd = NLLGrid(grdfile)
            # cache NLL grid
            _get_nll_grd.grds[key] = grd
            return grd
    raise RuntimeError
# dictionary to cache NLL grids
_get_nll_grd.grds = {}  # noqa


def _wave_arrival_nll(trace, phase, NLL_time_dir, focmec):
    """
    Travel time and takeoff angle using an NLL grid.

    :param trace: ObsPy Trace object
    :type trace: :class:`obspy.core.trace.Trace`
    :param phase: Phase name (e.g. 'P' or 'S')
    :type phase: str
    :param NLL_time_dir: Directory containing NLL grid files
    :type NLL_time_dir: str
    :param focmec: Whether to compute takeoff angle (focal mechanism usage)
    :type focmec: bool

    :return: Tuple (travel_time, takeoff_angle, method)
    :rtype: tuple(float or None, float or None, str)

    :raises RuntimeError: if NLL grid cannot be used/found
    """
    if NLL_time_dir is None:
        raise RuntimeError
    station = trace.stats.station
    travel_time = takeoff_angle = None
    grd_types = ['time']
    if focmec:
        grd_types.append('angle')
    for grd_type in grd_types:
        try:
            grd = _get_nll_grd(phase, station, grd_type, NLL_time_dir)
        except RuntimeError as e:
            logger.warning(
                f'{trace.id}: Cannot find NLL {grd_type} grid. '
                'Falling back to another method')
            raise RuntimeError from e
        if grd.station == 'DEFAULT':
            sta_x, sta_y = grd.project(
                trace.stats.coords.longitude, trace.stats.coords.latitude)
            grd.sta_x, grd.sta_y = sta_x, sta_y
        lon = trace.stats.event.hypocenter.longitude.value_in_deg
        lat = trace.stats.event.hypocenter.latitude.value_in_deg
        hypo_x, hypo_y = grd.project(lon, lat)
        hypo_z = trace.stats.event.hypocenter.depth.value_in_km
        if grd_type == 'time':
            travel_time = grd.get_value(hypo_x, hypo_y, hypo_z)
        elif grd_type == 'angle':
            _azimuth, takeoff_angle, _quality = grd.get_value(
                hypo_x, hypo_y, hypo_z)
    method = 'nlloc grid'
    return travel_time, takeoff_angle, method


def _wave_arrival_from_const_vel(trace, phase, vel):
    """
    Travel time and takeoff angle using a constant velocity (in km/s).

    :param trace: ObsPy Trace object
    :type trace: :class:`obspy.core.trace.Trace`
    :param phase: Phase for which to get arrival (e.g. 'P' or 'S')
    :type phase: str
    :param vel: Constant velocity in km/s
    :type vel: float

    :return: Tuple (travel_time, takeoff_angle, method)
    :rtype: tuple(float, float, str)

    :raises RuntimeError: if vel is None
    """
    if vel is None:
        raise RuntimeError
    travel_time = trace.stats.hypo_dist / vel
    takeoff_angle = degrees(asin(trace.stats.epi_dist / trace.stats.hypo_dist))
    # takeoff angle is 180° upwards and 0° downwards
    takeoff_angle = 180. - takeoff_angle
    method = f'constant V{phase.lower()}: {vel:.1f} km/s'
    return travel_time, takeoff_angle, method


def _wave_arrival_from_layer_model(trace, phase, vel):
    """
    Travel time and takeoff angle based on 1D velocity model

    :param trace: ObsPy Trace object
    :type trace: :class:`obspy.core.trace.Trace`
    :param phase: Phase for which to get arrival
    :type phase: str
    :param vel: Velocity model
    :type vel: dict with keys 'depths', 'P', 'S'

    :return: Travel time, takeoff angle and method used
    :rtype: tuple(float, float, str)
    """
    try:
        vmodel = CrustalVelocityModel(vel['depths'], vel['P'], vel['S'])
        hypo_z = trace.stats.event.hypocenter.depth.value_in_km
        Repi = trace.stats.epi_dist
        sta_z = -trace.stats.coords.elevation
        result = vmodel.calc_tt_and_angles(hypo_z, sta_z, Repi, wave=phase[:1])
        travel_time, takeoff_angle, _, _ = result
    except Exception as e:
        raise RuntimeError from e
    method = f'1D {phase.upper()} velocity model'
    return travel_time, takeoff_angle, method


def _wave_arrival_from_velocity_model(config, trace, phase):
    """
    Compute travel time and takeoff angle using layered or constant velocity.

    Uses a 1D layered model if layer_top_depths is set in config (via
    _wave_arrival_layer_model). Otherwise uses a constant velocity taken from
    config.vp or config.vs for the requested phase.

    :param config: Configuration object with velocity settings
    :type config: object
    :param trace: ObsPy Trace object
    :type trace: :class:`obspy.core.trace.Trace`
    :param phase: Phase name ('P' or 'S')
    :type phase: str

    :return: Tuple (travel_time, takeoff_angle, method)
    :rtype: tuple(float or None, float or None, str)
    """
    vel = {
        'depths': config.layer_top_depths,
        'P': config.vp,
        'S': config.vs
    }
    if vel['depths'] is not None:
        return _wave_arrival_from_layer_model(trace, phase, vel)
    if vel[phase] is not None:
        v_phase = vel[phase][0]
        return _wave_arrival_from_const_vel(trace, phase, v_phase)
    # if neither layered nor constant velocity is available, raise error
    raise RuntimeError


def _wave_arrival_taup(trace, phase):
    """
    Travel time and takeoff angle using TauP (global velocity model).

    :param trace: ObsPy Trace object
    :type trace: :class:`obspy.core.trace.Trace`
    :param phase: Phase name (e.g. 'P' or 'S')
    :type phase: str

    :return: Tuple (travel_time, takeoff_angle, method)
    :rtype: tuple(float, float, str)
    """
    phase_list = [phase.lower(), phase]
    hypo_depth = trace.stats.event.hypocenter.depth.value_in_km
    kwargs = {
        'source_depth_in_km': hypo_depth,
        'distance_in_degree': trace.stats.gcarc,
        'phase_list': phase_list
    }
    with warnings.catch_warnings(record=True) as warns:
        try:
            arrivals = model.get_travel_times(**kwargs)
        except Exception:
            trace.stats.event.hypocenter.depth = 0.
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
    travel_time = float(first_arrival.time)
    takeoff_angle = float(first_arrival.takeoff_angle)
    method = 'global velocity model (iasp91)'
    return travel_time, takeoff_angle, method


def _wave_arrival(trace, phase, config):
    """
    Get travel time and takeoff angle using configured methods.

    Tries NLL grid, then layered/constant velocity models,
    then TauP as fallback.

    :param trace: ObsPy Trace object
    :type trace: :class:`obspy.core.trace.Trace`
    :param phase: Phase name ('P' or 'S')
    :type phase: str
    :param config: Configuration object with model options
    :type config: object

    :return: Tuple (travel_time, takeoff_angle, method)
    :rtype: tuple(float or None, float or None, str)
    """
    NLL_time_dir = config.NLL_time_dir
    focmec = config.rp_from_focal_mechanism
    with contextlib.suppress(RuntimeError):
        return _wave_arrival_nll(trace, phase, NLL_time_dir, focmec)
    with contextlib.suppress(RuntimeError):
        return _wave_arrival_from_velocity_model(config, trace, phase)
    # if _wave_arrival_taup() fails, it will raise a RuntimeError
    return _wave_arrival_taup(trace, phase)


def _validate_pick(pick, theo_pick_time, tolerance, trace_id):
    """
    Check if a pick is valid, i.e., close enough to theoretical one.

    :param pick: Pick object (with .time and .phase)
    :type pick: object
    :param theo_pick_time: Theoretical pick time (UTCDateTime) or None
    :type theo_pick_time: obspy.UTCDateTime or None
    :param tolerance: Maximum allowed difference in seconds
    :type tolerance: float
    :param trace_id: Trace identifier for logging
    :type trace_id: str

    :return: True if pick within tolerance or if no theoretical time provided
    :rtype: bool
    """
    if theo_pick_time is None:
        return True
    delta_t = pick.time - theo_pick_time
    if abs(delta_t) > tolerance:  # seconds
        logger.warning(
            f'{trace_id}: measured {pick.phase} pick time - theoretical time '
            f'= {delta_t:.1f} s, above tolerance of {tolerance:.1f} s')
        return False
    return True


def _get_theo_pick_time(trace, travel_time):
    """
    Compute the theoretical pick time from hypocenter origin and travel time.

    :param trace: ObsPy Trace object
    :type trace: :class:`obspy.core.trace.Trace`
    :param travel_time: Travel time in seconds
    :type travel_time: float

    :return: Theoretical pick time (origin_time + travel_time) or None
    :rtype: obspy.UTCDateTime or None
    """
    if trace.stats.event.hypocenter.origin_time is None:
        msg = (
            f'{trace.id}: hypocenter origin time not set: '
            'unable to compute theoretical pick time')
        if msg not in _get_theo_pick_time.msg_cache:
            _get_theo_pick_time.msg_cache.append(msg)
            logger.warning(msg)
        return None
    return trace.stats.event.hypocenter.origin_time + travel_time
_get_theo_pick_time.msg_cache = [] # noqa


def _travel_time_from_pick(trace, pick_time):
    """
    Compute travel time from a pick time and the event origin time.

    :param trace: ObsPy Trace object
    :type trace: :class:`obspy.core.trace.Trace`
    :param pick_time: Pick time (UTCDateTime)
    :type pick_time: obspy.UTCDateTime

    :return: Travel time in seconds or None if origin time missing
    :rtype: float or None
    """
    try:
        travel_time = pick_time - trace.stats.event.hypocenter.origin_time
    except TypeError:
        travel_time = None
    return travel_time


def _find_picks(trace, phase, theo_pick_time, tolerance):
    """
    Search for valid picks in trace stats. Return pick time if found.

    :param trace: ObsPy Trace object
    :type trace: :class:`obspy.core.trace.Trace`
    :param phase: Phase name to search for
    :type phase: str
    :param theo_pick_time: Theoretical pick time or None
    :type theo_pick_time: obspy.UTCDateTime or None
    :param tolerance: Tolerance in seconds for validating picks
    :type tolerance: float

    :return: Pick time if a valid pick is found, otherwise None
    :rtype: obspy.UTCDateTime or None
    """
    for pick in (p for p in trace.stats.picks if p.phase.upper() == phase):
        if _validate_pick(pick, theo_pick_time, tolerance, trace.id):
            trace.stats.arrivals[phase] = (phase, pick.time)
            return pick.time
    return None


def add_arrival_to_trace(trace, phase, config):
    """
    Add arrival time, travel time and takeoff angle to trace for the
    given phase.

    Uses the theoretical arrival time if no pick is available
    or if the pick is too different from the theoretical arrival.

    :param trace: ObsPy Trace object to augment
    :type trace: :class:`obspy.core.trace.Trace`
    :param phase: Phase name ('P' or 'S')
    :type phase: str
    :param config: Configuration object with tolerances and model options
    :type config: object

    :return: None
    :rtype: None
    """
    tolerance = (
        config.p_arrival_tolerance
        if phase == 'P' else
        config.s_arrival_tolerance
    )
    key = f'{trace.id}_{phase}'
    trst = trace.stats
    # First, see if there are cached values
    with contextlib.suppress(KeyError):
        trst.arrivals[phase] = add_arrival_to_trace.pick_cache[key]
        trst.travel_times[phase] = add_arrival_to_trace.travel_time_cache[key]
        trst.takeoff_angles[phase] = add_arrival_to_trace.angle_cache[key]
        return
    # If no cache is available, compute travel_time and takeoff_angle
    travel_time, takeoff_angle, method = _wave_arrival(trace, phase, config)
    logger.info(
        f'{trace.id}: {phase} travel time: {travel_time:.2f} s, '
        f'computed from {method}')
    theo_pick_time = _get_theo_pick_time(trace, travel_time)
    pick_time = _find_picks(trace, phase, theo_pick_time, tolerance)
    if pick_time is not None:
        logger.info(f'{trace.id}: found {phase} pick')
        travel_time = \
            _travel_time_from_pick(trace, pick_time) or travel_time
        pick_phase = phase
    else:
        logger.info(f'{trace.id}: no {phase} pick found')
        if theo_pick_time is None:
            raise ValueError(
                f'{trace.id}: no theoretical {phase} pick time available')
        logger.info(
            f'{trace.id}: using theoretical {phase} pick from {method}')
        pick_time = theo_pick_time
        pick_phase = f'{phase}theo'
    if config.rp_from_focal_mechanism:
        logger.info(
            f'{trace.id}: {phase} takeoff angle: {takeoff_angle:.1f} '
            f'computed from {method}')
    add_arrival_to_trace.pick_cache[key] =\
        trst.arrivals[phase] = (pick_phase, pick_time)
    add_arrival_to_trace.travel_time_cache[key] =\
        trst.travel_times[phase] = travel_time
    add_arrival_to_trace.angle_cache[key] =\
        trst.takeoff_angles[phase] = takeoff_angle
add_arrival_to_trace.pick_cache = {}  # noqa
add_arrival_to_trace.travel_time_cache = {}  # noqa
add_arrival_to_trace.angle_cache = {}  # noqa
