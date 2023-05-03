# -*- coding: utf-8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Compute radiation pattern.

:copyright:
    2021-2023 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
import contextlib
import logging
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
        raise ValueError(f'Unknown wave type: {wave}')
    return R


def _rad_patt_P(phi, dip, rake, takeoff):
    return (
        cos(rake) * sin(dip) * sin(takeoff)**2 * sin(2*phi) -
        cos(rake) * cos(dip) * sin(2*takeoff) * cos(phi) +
        sin(rake) * sin(2*dip) *
        (cos(takeoff)**2 - sin(takeoff)**2 * sin(phi)**2) +
        sin(rake) * cos(2*dip) * sin(2*takeoff) * sin(phi)
    )


def _rad_patt_S(phi, dip, rake, takeoff):
    RSV = _rad_patt_SV(phi, dip, rake, takeoff)
    RSH = _rad_patt_SH(phi, dip, rake, takeoff)
    return (RSV**2. + RSH**2.)**(1./2)


def _rad_patt_SV(phi, dip, rake, takeoff):
    return (
        sin(rake) * cos(2*dip) * cos(2*takeoff) * sin(phi) -
        cos(rake) * cos(dip) * cos(2*takeoff) * cos(phi) +
        0.5 * cos(rake) * sin(dip) * sin(2*takeoff) * sin(2*phi) -
        0.5 * sin(rake) * sin(2*dip) * sin(2*takeoff) *
        (1 + sin(phi)**2)
    )


def _rad_patt_SH(phi, dip, rake, takeoff):
    return (
        cos(rake) * cos(dip) * cos(takeoff) * sin(phi) +
        cos(rake) * sin(dip) * sin(takeoff) * cos(2*phi) +
        sin(rake) * cos(2*dip) * cos(takeoff) * cos(phi) -
        0.5 * sin(rake) * sin(2*dip) * sin(takeoff) * sin(2*phi)
    )


# Cache radiation patterns
rp_cache = {}
# Cache messages to avoid duplication
rp_msg_cache = []


def get_radiation_pattern_coefficient(stats, config):
    global rp_cache
    global rp_msg_cache
    wave_type = config.wave_type  # P, S, SV, SH
    simple_wave_type = wave_type[0].lower()  # p or s
    if not config.rp_from_focal_mechanism:
        return config[f'rp{simple_wave_type}']
    try:
        strike = stats.event.focal_mechanism.strike
        dip = stats.event.focal_mechanism.dip
        rake = stats.event.focal_mechanism.rake
    except Exception:
        msg = (
            f'Cannot find focal mechanism. Using "rp{simple_wave_type}" '
            'value from config file'
        )
        if msg not in rp_msg_cache:
            logger.warning(msg)
            rp_msg_cache.append(msg)
        return config[f'rp{simple_wave_type}']
    traceid = '.'.join(
        (stats.network, stats.station, stats.location, stats.channel))
    key = f'{traceid}_{wave_type}'
    # try to get radiation pattern from cache
    with contextlib.suppress(KeyError):
        return rp_cache[key]
    try:
        takeoff_angle = stats.takeoff_angles[simple_wave_type.upper()]
    except Exception:
        msg = (
            f'{traceid}: Cannot find takeoff angle. '
            f'Using "rp{simple_wave_type}" value from config file'
        )
        if msg not in rp_msg_cache:
            logger.warning(msg)
            rp_msg_cache.append(msg)
        return config[f'rp{simple_wave_type}']
    rp = radiation_pattern(
        strike, dip, rake, takeoff_angle, stats.azimuth, wave_type)
    # we are interested only in amplitude
    # (P, SV and SH radiation patterns have a sign)
    rp = abs(rp)
    logger.info(
        f'{traceid}: {wave_type} radiation pattern from focal mechanism: '
        f'{rp:.2f}'
    )
    rp_cache[key] = rp
    return rp
