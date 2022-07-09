# -*- coding: utf-8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Compute radiation pattern.

:copyright:
    2021-2022 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
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


# Cache radiation patterns
rps_cache = dict()


def get_radiation_pattern_coefficient(stats, config):
    if not config.rps_from_focal_mechanism:
        if config.wave_type[0] == 'S':
            return config.rps
        elif config.wave_type[0] == 'P':
            return config.rpp
    try:
        strike = stats.hypo.strike
        dip = stats.hypo.dip
        rake = stats.hypo.rake
    except Exception:
        logger.warning(
            'Cannot find focal mechanism. Using "rp[p/s]" value from config file')
        if config.wave_type[0] == 'S':
            return config.rps
        elif config.wave_type[0] == 'P':
            return config.rpp
    traceid = '.'.join(
        (stats.network, stats.station, stats.location, stats.channel))
    wave = config.wave_type
    key = '{}_{}'.format(traceid, wave)
    global rps_cache
    try:
        rps = rps_cache[key]
        return rps
    except KeyError:
        pass
    try:
        takeoff_angle = stats.takeoff_angles[config.wave_type[0]]
    except Exception:
        logger.warning(
            '{}: Cannot find takeoff angle. '
            'Using "rps" value from config file'.format(traceid))
        return config.rps
    rps = radiation_pattern(
        strike, dip, rake, takeoff_angle, stats.azimuth, wave)
    # we are interested only in amplitude
    # (SV and SH radiation patterns have a sign)
    rps = abs(rps)
    logger.info(
        '{}: {} radiation pattern from focal mechanism: {:.2f}'.format(
            traceid, wave, rps
        ))
    rps_cache[key] = rps
    return rps
