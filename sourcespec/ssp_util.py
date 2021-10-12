# -*- coding: utf-8 -*-
"""
Utility functions for sourcespec.

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
import math
import numpy as np
from obspy.signal.invsim import cosine_taper as _cos_taper
from obspy.geodetics import gps2dist_azimuth, kilometers2degrees
from sourcespec.ssp_setup import ssp_exit
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


def moment_to_mag(moment):
    """Convert moment to magnitude."""
    return (np.log10(moment) - 9.1) / 1.5


def mag_to_moment(magnitude):
    """Convert magnitude to moment."""
    return np.power(10, (1.5 * magnitude + 9.1))


def select_trace(stream, traceid, instrtype):
    """Select trace from stream using traceid and instrument type."""
    return [tr for tr in stream.select(id=traceid)
            if tr.stats.instrtype == instrtype][0]


def get_vel(lon, lat, depth, wave, config):
    """
    Get velocity at a given point from NonLinLoc grid.

    Fall back to config.vp or config.vs if no grid is defined.
    """
    if wave == 'P':
        vel = config.vp
    elif wave == 'S':
        vel = config.vs
    if config.NLL_model_dir is None:
        return vel
    # Lazy-import here, since nllgrid is not an installation requirement
    from nllgrid import NLLGrid
    grdfile = '*.{}.mod.hdr'.format(wave)
    grdfile = os.path.join(config.NLL_model_dir, grdfile)
    try:
        grdfile = glob(grdfile)[0]
    except IndexError:
        return vel
    grd = NLLGrid(grdfile)
    x, y = grd.project(lon, lat)
    if grd.type == 'SLOW_LEN':
        slow_len = grd.get_value(lon, lat, depth)
        vel = grd.dx / slow_len
    return vel
# -----------------------------------------------------------------------------


# SIGNAL ANALYSIS -------------------------------------------------------------
def cosine_taper(signal, width):
    # TODO: this taper looks more like a hanning...
    npts = len(signal)
    p = 2 * width
    tap = _cos_taper(npts, p)
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
        w = eval('np.'+window+'(window_len)')
    y = np.convolve(w/w.sum(), s, mode='same')

    yy = y[window_len:-window_len+1]
    # check if there are NaN values
    nanindexes = np.where(np.isnan(yy))
    yy[nanindexes] = x[nanindexes]
    return yy


def remove_instr_response(trace, correct='True',
                          pre_filt=(0.5, 0.6, 40., 45.)):
    if correct == 'False':
        return trace
    traceId = trace.get_id()
    paz = trace.stats.paz
    if paz is None:
        logger.warning('%s: no poles and zeros for trace' % traceId)
        return None

    # remove the mean...
    trace.detrend(type='constant')
    # ...and the linear trend
    trace.detrend(type='linear')

    # Finally remove instrument response
    # If we don't have any pole or zero defined,
    # then be smart and just use the sensitivity ;)
    if len(paz.poles) == 0 and len(paz.zeros) == 0:
        correct = 'sensitivity_only'
    if correct == 'sensitivity_only':
        trace.data /= paz.sensitivity
    # Otherwhise we need to call trace.simulate(), which is quite slow...
    else:
        with warnings.catch_warnings(record=True) as w:
            # N.B. using "sacsim=True" makes this two times slower!
            # (because of "c_sac_taper()")
            # TODO: fill up a bug on obspy.org
            trace.simulate(paz_remove=paz, paz_simulate=None,
                           remove_sensitivity=True, simulate_sensitivity=None,
                           pre_filt=pre_filt, sacsim=False)
            if len(w) > 0:
                logger.warning('%s: remove_instr_response: %s' %
                               (trace.stats.station, w[-1].message))
    return trace
# -----------------------------------------------------------------------------


# GEODETICS AND COORDINATES ---------------------------------------------------
def toRad(degrees):
    radians = math.pi * degrees / 180
    return radians


def toDeg(radians):
    degrees = 180 * radians / math.pi
    return degrees


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
