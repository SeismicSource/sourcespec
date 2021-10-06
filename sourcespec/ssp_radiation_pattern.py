
# -*- coding: utf-8 -*-
"""
Compute radiation pattern.

:copyright:
    2021 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement, Version 2.1
    (http://www.cecill.info/index.en.html)
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import logging
import warnings
from math import pi, sin, cos
from obspy.taup import TauPyModel
logger = logging.getLogger(__name__.split('.')[-1])
model = TauPyModel(model='iasp91')


def toRad(angle):
    return angle/180. * pi


def radiation_pattern(strike, dip, rake, takeoff_angle, azimuth, wave):
    """
    Body wave radiation pattern.

    From Lay-Wallace, page 340.
    """
    strike = toRad(strike)
    dip = toRad(dip)
    rake = toRad(rake)
    takeoff_angle = toRad(takeoff_angle)
    azimuth = toRad(azimuth)
    phi = azimuth - strike

    if wave == 'P':
        R = _rad_patt_P(phi, dip, rake, takeoff_angle)
    elif wave == 'S':
        R = _rad_patt_S(phi, dip, rake, takeoff_angle)
    elif wave == 'SV':
        R = _rad_patt_SV(phi, dip, rake, takeoff_angle)
    elif wave == 'SH':
        R = _rad_patt_SH(phi, dip, rake, takeoff_angle)
    else:
        raise ValueError('Unknown wave type: {}'.format(wave))
    return R


def _rad_patt_P(phi, dip, rake, takeoff):
    R = cos(rake) * sin(dip) * sin(takeoff)**2 * sin(2*phi) -\
        cos(rake) * cos(dip) * sin(2*takeoff) * cos(phi) +\
        sin(rake) * sin(2*dip) *\
        (cos(takeoff)**2 - sin(takeoff)**2 * sin(phi)**2) +\
        sin(rake) * cos(2*dip) * sin(2*takeoff) * sin(phi)
    return R


def _rad_patt_S(phi, dip, rake, takeoff):
    RSV = _rad_patt_SV(phi, dip, rake, takeoff)
    RSH = _rad_patt_SH(phi, dip, rake, takeoff)
    R = (RSV**2. + RSH**2.)**(1./2)
    return R


def _rad_patt_SV(phi, dip, rake, takeoff):
    R = sin(rake) * cos(2*dip) * cos(2*takeoff) * sin(phi) -\
        cos(rake) * cos(dip) * cos(2*takeoff) * cos(phi) +\
        0.5 * cos(rake) * sin(dip) * sin(2*takeoff) * sin(2*phi) -\
        0.5 * sin(rake) * sin(2*dip) * sin(2*takeoff) *\
        (1 + sin(phi)**2)
    return R


def _rad_patt_SH(phi, dip, rake, takeoff):
    R = cos(rake) * cos(dip) * cos(takeoff) * sin(phi) +\
        cos(rake) * sin(dip) * sin(takeoff) * cos(2*phi) +\
        sin(rake) * cos(2*dip) * cos(takeoff) * cos(phi) -\
        0.5 * sin(rake) * sin(2*dip) * sin(takeoff) * sin(2*phi)
    return R


def get_radiation_pattern_coefficient(stats, config):
    if not config.rps_from_focal_mechanism:
        return config.rps
    try:
        strike = stats.hypo.strike
        dip = stats.hypo.dip
        rake = stats.hypo.rake
    except Exception:
        logger.warning(
            'Cannot find focal mechanism. Using "rps" value from config file')
        return config.rps
    phase_list = ['s', 'S']
    with warnings.catch_warnings(record=True) as warns:
        try:
            arrivals = model.get_travel_times(
                        source_depth_in_km=stats.hypo.depth,
                        distance_in_degree=stats.gcarc,
                        phase_list=phase_list)
        except Exception:
            stats.hypo.depth = 0.
            arrivals = model.get_travel_times(
                        source_depth_in_km=stats.hypo.depth,
                        distance_in_degree=stats.gcarc,
                        phase_list=phase_list)
        for w in warns:
            message = str(w.message)
            # Ignore a specific obspy.taup warning we do not care about
            if '#2280' in message:
                continue
            logger.warning(message)
    # get first arrival
    arrivals = sorted(arrivals, key=lambda a: a.time)
    first_arrival = arrivals[0]
    takeoff_angle = first_arrival.takeoff_angle
    rps = radiation_pattern(
        strike, dip, rake, takeoff_angle, stats.azimuth, 'S')
    traceid = '.'.join(
        (stats.network, stats.station, stats.location, stats.channel))
    logger.info(
        '{}: Radiation pattern from focal mechanism: {:.2f}'.format(
            traceid, rps
        ))
    return rps
