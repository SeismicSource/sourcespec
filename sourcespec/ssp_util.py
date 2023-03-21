# -*- coding: utf-8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Utility functions for sourcespec.

:copyright:
    2012-2023 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
import os
from glob import glob
import logging
import math
import numpy as np
from obspy.signal.invsim import cosine_taper as _cos_taper
from obspy.geodetics import gps2dist_azimuth, kilometers2degrees
from obspy.taup import TauPyModel
model = TauPyModel(model='iasp91')
v_model = model.model.s_mod.v_mod
logger = logging.getLogger(__name__.split('.')[-1])


# MISC ------------------------------------------------------------------------
def spec_minmax(amp, freq, amp_minmax=None, freq_minmax=None):
    amp_min = amp.min()
    amp_max = amp.max()
    if amp_minmax is None:
        amp_minmax = [amp_min, amp_max]
    else:
        if amp_min < amp_minmax[0]:
            amp_minmax[0] = amp_min
        if amp_max > amp_minmax[1]:
            amp_minmax[1] = amp_max
    freq_min = freq.min()
    freq_max = freq.max()
    if freq_minmax is None:
        freq_minmax = [freq_min, freq_max]
    else:
        if freq_min < freq_minmax[0]:
            freq_minmax[0] = freq_min
        if freq_max > freq_minmax[1]:
            freq_minmax[1] = freq_max
    return amp_minmax, freq_minmax


def select_trace(stream, traceid, instrtype):
    """Select trace from stream using traceid and instrument type."""
    return [tr for tr in stream.select(id=traceid)
            if tr.stats.instrtype == instrtype][0]


def _get_vel_from_config(wave, where, config):
    if wave not in ['P', 'S']:
        msg = f'Invalid wave type: {wave}'
        raise ValueError(msg)
    if where not in ['source', 'stations']:
        msg = f'Invalid location type: {wave}'
        raise ValueError(msg)
    if wave == 'P':
        if where == 'source':
            vel = config.vp_source
        elif where == 'stations':
            vel = (
                config.vp_stations if config.vp_stations is not None
                else config.vp_source)
    elif wave == 'S':
        if where == 'source':
            vel = config.vs_source
        elif where == 'stations':
            vel = (
                config.vs_stations if config.vs_stations is not None
                else config.vs_source)
    return vel


def _get_vel_from_NLL(lon, lat, depth_in_km, wave, config):
    # Lazy-import here, since nllgrid is not an installation requirement
    from nllgrid import NLLGrid
    grdfile = f'*.{wave}.mod.hdr'
    grdfile = os.path.join(config.NLL_model_dir, grdfile)
    try:
        grdfile = glob(grdfile)[0]
    except IndexError as e:
        raise FileNotFoundError(f'Unable to find model file {grdfile}') from e
    grd = NLLGrid(grdfile)
    x, y = grd.project(lon, lat)
    value = grd.get_value(x, y, depth_in_km)
    if grd.type == 'VELOCITY':
        vel = value
    elif grd.type == 'VELOCITY_METERS':
        vel = value/1e3
    elif grd.type == 'SLOWNESS':
        vel = 1./value
    elif grd.type == 'SLOW_LEN':
        vel = grd.dx / value
    elif grd.type == 'VEL2':
        vel = value**0.5
    elif grd.type == 'SLOW2':
        vel = 1./(value**0.5)
    elif grd.type == 'SLOW2_METERS':
        vel = (1./(value**0.5))/1e3
    else:
        raise ValueError(f'Unsupported grid type: {grd.type}')
    logger.info(f'Using {wave} velocity from NLL model')
    return vel


def _get_vel_from_taup(depth_in_km, wave):
    if depth_in_km < 0:
        depth_in_km = 1e-3
    return v_model.evaluate_above(depth_in_km, wave)[0]


def get_vel(lon, lat, depth_in_km, wave, config):
    """Get velocity at a given point from NLL grid, config or taup model."""
    # If depth is large, we assume that we are close to the source
    if depth_in_km >= 2:
        vel = _get_vel_from_config(wave, 'source', config)
    else:
        vel = _get_vel_from_config(wave, 'stations', config)
    if vel is None and config.NLL_model_dir is None:
        vel = _get_vel_from_taup(depth_in_km, wave)
        logger.info(
            f'Using {wave} velocity from global velocity model (iasp91)')
    if config.NLL_model_dir is not None:
        try:
            vel = _get_vel_from_NLL(lon, lat, depth_in_km, wave, config)
        except Exception as msg:
            logger.warning(msg)
            logger.warning(f'Using {wave} velocity from config')
    return vel
# -----------------------------------------------------------------------------


# SOURCE PARAMETERS -----------------------------------------------------------
def moment_to_mag(moment):
    """Convert moment to magnitude."""
    return (np.log10(moment) - 9.1) / 1.5


def mag_to_moment(magnitude):
    """Convert magnitude to moment."""
    return np.power(10, (1.5 * magnitude + 9.1))


def source_radius(fc_in_hz, vs_in_m_per_s):
    """
    Compute source radius in meters.

    Madariaga (2009), doi:10.1007/978-1-4419-7695-6_22, eq. 31
    """
    return 0.3724 * vs_in_m_per_s / fc_in_hz


def bsd(Mo_in_N_m, ra_in_m):
    """
    Compute Brune stress drop in MPa.

    Madariaga (2009), doi:10.1007/978-1-4419-7695-6_22, eq. 27
    """
    return 7./16 * Mo_in_N_m / ra_in_m**3 * 1e-6


def quality_factor(travel_time_in_s, t_star_in_s):
    """Compute quality factor from travel time and t_star."""
    return np.inf if t_star_in_s == 0 else travel_time_in_s/t_star_in_s
# -----------------------------------------------------------------------------


# SIGNAL ANALYSIS -------------------------------------------------------------
def cosine_taper(signal, width, left_taper=False):
    # TODO: this taper looks more like a hanning...
    npts = len(signal)
    p = 2 * width
    tap = _cos_taper(npts, p)
    if left_taper:
        tap[npts // 2:] = 1.
    signal *= tap


# modified from: http://stackoverflow.com/q/5515720
def smooth(x, window_len=11, window='hanning'):
    if x.ndim != 1:
        raise ValueError('smooth only accepts 1 dimension arrays.')
    if x.size < window_len:
        raise ValueError('Input vector needs to be bigger than window size.')
    if window_len < 3:
        return x
    if window not in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError("Window is one of 'flat', 'hanning', 'hamming', "
                         "'bartlett', 'blackman'")
    s = np.r_[2*x[0]-x[window_len-1::-1], x, 2*x[-1]-x[-1:-window_len:-1]]
    if window == 'flat':  # moving average
        w = np.ones(window_len, 'd')
    else:
        w = eval(f'np.{window}(window_len)')
    y = np.convolve(w/w.sum(), s, mode='same')
    yy = y[window_len:-window_len+1]
    # check if there are NaN values
    nanindexes = np.where(np.isnan(yy))
    yy[nanindexes] = x[nanindexes]
    return yy


def remove_instr_response(trace, pre_filt=(0.5, 0.6, 40., 45.)):
    """
    Remove instrument response from a trace.

    Trace is converted to the sensor units (m for a displacement sensor,
    m/s for a short period or broadband velocity sensor, m/s**2 for a strong
    motion sensor).

    :param trace: Trace to be corrected.
    :type trace: :class:`~obspy.core.trace.Trace`
    :param pre_filt: Pre-filter frequencies (``None`` means no pre-filtering).
    :type pre_filt: tuple of four floats
    """
    trace_info = trace.stats.info
    inventory = trace.stats.inventory
    if not inventory:
        # empty inventory
        raise RuntimeError(f'{trace_info}: no instrument response for trace')
    # remove the mean...
    trace.detrend(type='constant')
    # ...and the linear trend
    trace.detrend(type='linear')
    # Define output units based on nominal units in inventory
    # Note: ObsPy >= 1.3.0 supports the 'DEF' output unit, which will make
    #       this step unnecessary
    if trace.stats.units.lower() == 'm':
        output = 'DISP'
    if trace.stats.units.lower() == 'm/s':
        output = 'VEL'
    if trace.stats.units.lower() == 'm/s**2':
        output = 'ACC'
    # Finally remove instrument response,
    # trace is converted to the sensor units
    trace.remove_response(
        inventory=inventory, output=output, pre_filt=pre_filt)
    if any(np.isnan(trace.data)):
        raise RuntimeError(
            f'{trace_info}: NaN values in trace after '
            'instrument response removal')
# -----------------------------------------------------------------------------


# GEODETICS AND COORDINATES ---------------------------------------------------
def toRad(degrees):
    return math.pi * degrees / 180


def toDeg(radians):
    return 180 * radians / math.pi


def hypo_dist(trace):
    """Compute hypocentral and epicentral distance (in km) for a trace."""
    try:
        coords = trace.stats.coords
        hypo = trace.stats.hypo
    except (KeyError, AttributeError):
        return None
    if None in (coords, hypo):
        return None
    stla = coords.latitude
    stlo = coords.longitude
    stel = coords.elevation
    evla = hypo.latitude
    evlo = hypo.longitude
    evdp = hypo.depth
    if None in (stla, stlo, stel, evla, evlo, evdp):
        return None
    epi_dist, az, baz = gps2dist_azimuth(
        hypo.latitude, hypo.longitude,
        trace.stats.coords.latitude, trace.stats.coords.longitude)
    epi_dist /= 1e3   # in km
    gcarc = kilometers2degrees(epi_dist)
    hypo_dist = math.sqrt(epi_dist**2 + (stel+evdp)**2)
    trace.stats.azimuth = az
    trace.stats.back_azimuth = baz
    trace.stats.epi_dist = epi_dist
    trace.stats.hypo_dist = hypo_dist
    trace.stats.gcarc = gcarc
    return hypo_dist
# -----------------------------------------------------------------------------
