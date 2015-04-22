# -*- coding: utf-8 -*-
# ssp_util.py
#
# (c) 2012-2015 Claudio Satriano <satriano@ipgp.fr>
'''
Utility functions for source_spec.
'''
from __future__ import division
import logging
import warnings
import math
import numpy as np
from obspy.signal import cosTaper


# MISC ------------------------------------------------------------------------
def spec_minmax(amp, freq, amp_minmax=None, freq_minmax=None):
    amp_min = amp.min()
    amp_max = amp.max()
    if amp_minmax is None:
        amp_minmax = [amp_min, amp_max]
    else:
        if amp_min < amp_minmax[0]: amp_minmax[0] = amp_min
        if amp_max > amp_minmax[1]: amp_minmax[1] = amp_max
    freq_min = freq.min()
    freq_max = freq.max()
    if freq_minmax is None:
        freq_minmax = [freq_min, freq_max]
    else:
        if freq_min < freq_minmax[0]: freq_minmax[0] = freq_min
        if freq_max > freq_minmax[1]: freq_minmax[1] = freq_max
    return amp_minmax, freq_minmax


def wave_arrival(trace, vel, phase):
    '''
    Obtain arrival time for a given phase and a given trace.

    Returns the theoretical arrival time if no pick is available
    or if the pick is too different from the theoretical arrival.
    '''
    try:
        arrival_secs = trace.stats.arrivals[phase][1]
        return trace.stats.starttime + arrival_secs
    except KeyError:
        pass
    if trace.stats.hypo.origin_time is not None:
        theo_pick_time =\
            trace.stats.hypo.origin_time + trace.stats.hypo_dist / vel
        trace.stats.arrivals[phase] = (phase + 'theo', theo_pick_time - trace.stats.starttime)
    else:
        theo_pick_time = None
    for pick in trace.stats.picks:
        if pick.phase == phase:
            if theo_pick_time is not None:
                delta_t = pick.time - theo_pick_time
                if abs(delta_t) > 4.: # seconds #TODO parametrize?
                    logging.warning(
                        '%s: measured %s pick time - theoretical time = %.1f s. Using theoretical'
                        % (trace.id, phase, delta_t))
                    continue
            logging.info('%s: found %s pick' % (trace.id, phase))
            trace.stats.arrivals[phase] = (phase, pick.time - trace.stats.starttime)
            return pick.time
    logging.info('%s: using theoretical %s pick' % (trace.id, phase))
    return theo_pick_time


def moment_to_mag(data):
    ''' Convert moment to magnitude.'''
    return (np.log10(data) - 9.1 ) / 1.5


def mag_to_moment(data):
    ''' Convert magnitude to moment.'''
    return np.power(10, (1.5 * data + 9.1))


def select_trace(stream, traceid, instrtype):
    '''
    Select trace from stream using traceid
    and instrument type.
    '''
    return [tr for tr in stream.select(id=traceid)
            if tr.stats.instrtype == instrtype][0]
# -----------------------------------------------------------------------------


# SIGNAL ANALYSIS -------------------------------------------------------------
def cosine_taper(signal, width):
    #TODO: this taper looks more like a hanning...
    npts = len(signal)
    p = 2 * width
    tap = cosTaper(npts, p)
    signal *= tap


# modified from: http://stackoverflow.com/questions/5515720/python-smooth-time-series-data
def smooth(x, window_len=11, window='hanning'):
    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."
    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."
    if window_len<3:
        return x
    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"
    s = np.r_[2*x[0]-x[window_len-1::-1],x,2*x[-1]-x[-1:-window_len:-1]]
    if window == 'flat': #moving average
        w = np.ones(window_len,'d')
    else:
        w = eval('np.'+window+'(window_len)')
    y = np.convolve(w/w.sum(),s,mode='same')

    yy = y[window_len:-window_len+1]
    # check if there are NaN values
    nanindexes = np.where(np.isnan(yy))
    yy[nanindexes] = x[nanindexes]
    return yy


def remove_instr_response(trace, correct='True', pre_filt=(0.5, 0.6, 40., 45.)):
    if correct == 'False':
        return trace
    traceId = trace.getId()
    paz = trace.stats.paz
    if paz is None:
        logging.warning('%s: no poles and zeros for trace' % traceId)
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
            if len(w)>0:
                logging.warning('%s: remove_instr_response: %s'
                        % (trace.stats.station, w[-1].message))
    return trace
# -----------------------------------------------------------------------------


# GEODETICS AND COORDINATES ---------------------------------------------------
def toRad(degrees):
    radians = math.pi * degrees / 180
    return radians


def toDeg(radians):
    degrees = 180 * radians / math.pi
    return degrees


def calc_dist(lat1, lon1, lat2, lon2):
    '''
    Distance between two point on the earth.
    Haversine formula:
    http://www.movable-type.co.uk/scripts/latlong.html
    '''
    R = 6371 # km
    dLat = toRad(lat2-lat1)
    dLon = toRad(lon2-lon1)
    a = math.sin(dLat/2) * math.sin(dLat/2) + \
        math.cos(toRad(lat1)) * math.cos(toRad(lat2)) * \
        math.sin(dLon/2) * math.sin(dLon/2)
    gcarc = 2 * math.atan2(math.sqrt(a), math.sqrt(1-a))
    dist = R * gcarc
    return dist, toDeg(gcarc)


def hypo_dist(trace):
    try:
        coords = trace.stats.coords
        hypo = trace.stats.hypo
    except (KeyError, AttributeError):
        return None
    if coords is None or hypo is None:
        return None
    stla = coords.latitude
    stlo = coords.longitude
    stel = coords.elevation
    evla = hypo.latitude
    evlo = hypo.longitude
    evdp = hypo.depth
    if stla is None or stlo is None or stel is None\
        or evla is None or evlo is None or evdp is None:
            return None
    epi_dist, gcarc = calc_dist(stla, stlo, evla, evlo)
    hypo_dist = math.sqrt(epi_dist**2 + (stel+evdp)**2)
    trace.stats.hypo_dist = hypo_dist
    trace.stats.gcarc = gcarc
    return hypo_dist
# -----------------------------------------------------------------------------
