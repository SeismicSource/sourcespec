# -*- coding: utf-8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Compute radiation pattern.

:copyright:
    2021-2025 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
import contextlib
import logging
from math import pi, sin, cos
from obspy.taup import TauPyModel
logger = logging.getLogger(__name__.rsplit('.', maxsplit=1)[-1])
model = TauPyModel(model='iasp91')


def toRad(angle):
    """Convert angle from degrees to radians."""
    return angle / 180. * pi


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
        raise ValueError(f'Unknown wave type: {wave}')
    return R


def _rad_patt_P(phi, dip, rake, takeoff):
    return (
        cos(rake) * sin(dip) * sin(takeoff)**2 * sin(2 * phi) -
        cos(rake) * cos(dip) * sin(2 * takeoff) * cos(phi) +
        sin(rake) * sin(2 * dip) *
        (cos(takeoff)**2 - sin(takeoff)**2 * sin(phi)**2) +
        sin(rake) * cos(2 * dip) * sin(2 * takeoff) * sin(phi)
    )


def _rad_patt_S(phi, dip, rake, takeoff):
    RSV = _rad_patt_SV(phi, dip, rake, takeoff)
    RSH = _rad_patt_SH(phi, dip, rake, takeoff)
    return (RSV**2. + RSH**2.)**(1. / 2)


def _rad_patt_SV(phi, dip, rake, takeoff):
    return (
        sin(rake) * cos(2 * dip) * cos(2 * takeoff) * sin(phi) -
        cos(rake) * cos(dip) * cos(2 * takeoff) * cos(phi) +
        0.5 * cos(rake) * sin(dip) * sin(2 * takeoff) * sin(2 * phi) -
        0.5 * sin(rake) * sin(2 * dip) * sin(2 * takeoff) *
        (1 + sin(phi)**2)
    )


def _rad_patt_SH(phi, dip, rake, takeoff):
    return (
        cos(rake) * cos(dip) * cos(takeoff) * sin(phi) +
        cos(rake) * sin(dip) * sin(takeoff) * cos(2 * phi) +
        sin(rake) * cos(2 * dip) * cos(takeoff) * cos(phi) -
        0.5 * sin(rake) * sin(2 * dip) * sin(takeoff) * sin(2 * phi)
    )


# Cache radiation patterns
RP_CACHE = {}
# Cache messages to avoid duplication
RP_MSG_CACHE = []


def get_radiation_pattern_coefficient(stats, config):
    """Get radiation pattern coefficient."""
    wave_type = config.wave_type  # P, S, SV, SH
    simple_wave_type = wave_type[0].lower()  # p or s
    if not config.rp_from_focal_mechanism:
        return config[f'rp{simple_wave_type}']
    try:
        fm = stats.event.focal_mechanism
        # will raise an Exception if any of these is None
        strike = float(fm.strike)
        dip = float(fm.dip)
        rake = float(fm.rake)
    except Exception:
        msg = (
            f'Cannot find focal mechanism. Using "rp{simple_wave_type}" '
            'value from config file'
        )
        if msg not in RP_MSG_CACHE:
            logger.warning(msg)
            RP_MSG_CACHE.append(msg)
        return config[f'rp{simple_wave_type}']
    traceid = '.'.join(
        (stats.network, stats.station, stats.location, stats.channel))
    key = f'{traceid}_{wave_type}'
    # try to get radiation pattern from cache
    with contextlib.suppress(KeyError):
        return RP_CACHE[key]
    try:
        takeoff_angle = stats.takeoff_angles[simple_wave_type.upper()]
    except KeyError:
        msg = (
            f'{traceid}: Cannot find takeoff angle. '
            f'Using "rp{simple_wave_type}" value from config file'
        )
        if msg not in RP_MSG_CACHE:
            logger.warning(msg)
            RP_MSG_CACHE.append(msg)
        return config[f'rp{simple_wave_type}']
    rp = radiation_pattern(
        strike, dip, rake, takeoff_angle, stats.azimuth, wave_type)
    # we are interested only in amplitude
    # (P, SV and SH radiation patterns have a sign)
    rp = abs(rp)
    # Apply water level
    rp = max(rp, config.rp_lower_bound)
    logger.info(
        f'{traceid}: {wave_type} radiation pattern from focal mechanism: '
        f'{rp:.2f}'
    )
    RP_CACHE[key] = rp
    return rp
