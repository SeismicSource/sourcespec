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
logger = logging.getLogger(__name__.rsplit('.', maxsplit=1)[-1])


# MISC ------------------------------------------------------------------------
def spec_minmax(amp, freq, amp_minmax=None, freq_minmax=None):
    """Get minimum and maximum values of spectral amplitude and frequency."""
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
# -----------------------------------------------------------------------------


# MEDIUM PROPERTIES -----------------------------------------------------------
def _get_medium_property_from_config(medium_property, where, config):
    """
    Get medium property from config.

    :param medium_property: Property to be retrieved
                     (``'vp'``, ``'vs'`` or ``'rho'``).
    :type medium_property: str
    :param where: Location type (``'source'`` or ``'stations'``).
    :type where: str
    :param config: Configuration object.
    :type config: :class:`~sourcespec.config.Config`
    :return: Property value.
    :rtype: float
    """
    if medium_property not in ['vp', 'vs', 'rho']:
        raise ValueError(f'Invalid property: {medium_property}')
    if where not in ['source', 'stations']:
        raise ValueError(f'Invalid location type: {where}')
    value = config[f'{medium_property}_{where}']
    # Use source property if stations property is not defined
    if where == 'stations' and value is None:
        value = config[f'{medium_property}_source']
    return value


def _get_vel_from_NLL(lon, lat, depth_in_km, wave, config):
    """
    Get velocity from NLL model.

    :param lon: Longitude (degrees).
    :type lon: float
    :param lat: Latitude (degrees).
    :type lat: float
    :param depth_in_km: Depth (km).
    :type depth_in_km: float
    :param wave: Wave type (``'P'`` or ``'S'``).
    :type wave: str
    :param config: Configuration object.
    :type config: :class:`~sourcespec.config.Config`
    :return: Velocity (km/s).
    :rtype: float
    """
    # Lazy-import here, since nllgrid is not an installation requirement
    # pylint: disable=import-outside-toplevel
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
        vel = value / 1e3
    elif grd.type == 'SLOWNESS':
        vel = 1. / value
    elif grd.type == 'SLOW_LEN':
        vel = grd.dx / value
    elif grd.type == 'VEL2':
        vel = value**0.5
    elif grd.type == 'SLOW2':
        vel = 1. / (value**0.5)
    elif grd.type == 'SLOW2_METERS':
        vel = (1. / (value**0.5)) / 1e3
    else:
        raise ValueError(f'Unsupported grid type: {grd.type}')
    logger.info(f'Using {wave} velocity from NLL model')
    return vel


def _get_medium_property_from_taup(depth_in_km, medium_property):
    """
    Get medium property (P- or S-wave velocity, density) at a given depth
    from taup model.

    :param depth_in_km: Depth (km).
    :type depth_in_km: float
    :param medium_property: Property to be retrieved
                     (``'vp'``, ``'vs'`` or ``'rho'``).
    :type medium_property: str
    :return: Property value.
    :rtype: float
    """
    # avoid negative depths
    depth_in_km = max(depth_in_km, 1e-2)
    try:
        prop = {'vp': 'p', 'vs': 's', 'rho': 'r'}[medium_property]
    except KeyError as e:
        raise ValueError(f'Invalid property: {medium_property}') from e
    value = v_model.evaluate_above(depth_in_km, prop)[0]
    if medium_property == 'rho':
        value *= 1e3  # convert g/cm**3 to kg/m**3
    return value


def get_medium_property(lon, lat, depth_in_km, medium_property, config):
    """
    Get medium property (P- or S-wave velocity, density) at a given point
    from NLL grid, config or taup model.

    :param lon: Longitude (degrees).
    :type lon: float
    :param lat: Latitude (degrees).
    :type lat: float
    :param depth_in_km: Depth (km).
    :type depth_in_km: float
    :param medium_property: Property to be retrieved
                     (``'vp'``, ``'vs'`` or ``'rho'``).
    :type medium_property: str
    :param config: Configuration object.
    :type config: :class:`~sourcespec.config.Config`
    :return: Property value.
    :rtype: float
    """
    if medium_property not in ['vp', 'vs', 'rho']:
        raise ValueError(f'Invalid property: {medium_property}')
    try:
        depth_in_km = float(depth_in_km)
    except Exception as e:
        raise ValueError(f'Invalid depth: {depth_in_km}') from e
    # If depth is large, we assume that we are close to the source
    if depth_in_km >= 2:
        value = _get_medium_property_from_config(
            medium_property, 'source', config)
    else:
        value = _get_medium_property_from_config(
            medium_property, 'stations', config)
    if config.NLL_model_dir is not None and medium_property in ['vp', 'vs']:
        wave = 'P' if medium_property == 'vp' else 'S'
        try:
            value = _get_vel_from_NLL(lon, lat, depth_in_km, wave, config)
        except Exception as msg:
            logger.warning(msg)
            logger.warning(f'Using {wave} velocity from config')
    if value is None:
        value = _get_medium_property_from_taup(depth_in_km, medium_property)
        logger.info(
            f'Using {medium_property} from global model (iasp91)')
    return value


def medium_property_string(medium_property, value):
    """
    Return a string with the property name and value.

    :param property: Property name. Must contain one of the following:
                     ``'vp'``, ``'vs'``, ``'rho'``, ``'depth'``
    :type property: str
    :param value: Property value.
    :type value: float
    :return: Property string.
    :rtype: str
    """
    if 'vp' in medium_property or 'vs' in medium_property:
        value = round(value, 2)
        unit = 'km/s'
    elif 'rho' in medium_property:
        value = round(value, 1)
        unit = 'kg/m3'
    elif 'depth' in medium_property:
        value = round(value, 1)
        unit = 'km'
    else:
        raise ValueError(f'Invalid property: {medium_property}')
    if value.is_integer():
        value = int(value)
    return f'{medium_property}: {value} {unit}'
# -----------------------------------------------------------------------------


# GEOMETRICAL SPREADING -------------------------------------------------------
def geom_spread_r_power_n(hypo_dist_in_km, exponent):
    """
    rⁿ geometrical spreading coefficient.

    :param hypo_dist_in_km: Hypocentral distance (km).
    :type hypo_dist_in_km: float
    :param exponent: Exponent.
    :type exponent: float
    :return: Geometrical spreading correction (in m)
    :rtype: float
    """
    dist = hypo_dist_in_km * 1e3
    return dist**exponent


def _boatwright_above_cutoff_dist(freqs, cutoff_dist, dist):
    """
    Geometrical spreading coefficient from Boatwright et al. (2002), eq. 8.

    Except that we take the square root of eq. 8, since we correct amplitude
    and not energy.

    This is the part of the equation that is valid for distances above the
    cutoff distance.

    :param freqs: Frequencies (Hz).
    :type freqs: numpy.ndarray
    :param cutoff_dist: Cutoff distance (m).
    :type cutoff_dist: float
    :param dist: Distance (m).
    :type dist: float
    :return: Geometrical spreading correction (in m)
    :rtype: numpy.ndarray
    """
    exponent = np.ones_like(freqs)
    low_freq = freqs <= 0.2
    mid_freq = np.logical_and(freqs > 0.2, freqs <= 0.25)
    high_freq = freqs >= 0.25
    exponent[low_freq] = 0.5
    exponent[mid_freq] = 0.5 + 2 * np.log10(5 * freqs[mid_freq])
    exponent[high_freq] = 0.7
    return cutoff_dist * (dist / cutoff_dist)**exponent


def geom_spread_boatwright(hypo_dist_in_km, cutoff_dist_in_km, freqs):
    """"
    Geometrical spreading coefficient from Boatwright et al. (2002), eq. 8.

    Except that we take the square root of eq. 8, since we correct amplitude
    and not energy.

    :param hypo_dist_in_km: Hypocentral distance (km).
    :type hypo_dist_in_km: float
    :param cutoff_dist_in_km: Cutoff distance (km).
    :type cutoff_dist_in_km: float
    :param freqs: Frequencies (Hz).
    :type freqs: numpy.ndarray
    :return: Geometrical spreading correction (in m)
    :rtype: numpy.ndarray
    """
    dist = hypo_dist_in_km * 1e3
    cutoff_dist = cutoff_dist_in_km * 1e3
    if dist <= cutoff_dist:
        return dist
    return _boatwright_above_cutoff_dist(freqs, cutoff_dist, dist)


def _compute_dtdd(angular_distance, aperture, source_depth_in_km, phase_list):
    """
    Compute the local derivative of takeoff angle with respect to angular
    distance.

    This is used to compute the geometrical spreading coefficient for
    teleseismic body waves, following Okal (1992), eq. 4.

    :param angular_distance: Angular distance (degrees).
    :type angular_distance: float
    :param aperture: Aperture of the ray tube (degrees).
    :type aperture: float
    :param source_depth_in_km: Source depth (km).
    :type source_depth_in_km: float
    :param phase_list: List of phases.
    :type phase_list: list of str
    :return: Local derivative of takeoff angle with respect to angular
             distance.
    :rtype: float
    """
    distances = np.linspace(
        angular_distance - aperture, angular_distance + aperture, 3)
    # pylint: disable=no-member
    takeoff_angles = np.array([
        model.get_travel_times(
            source_depth_in_km, d, phase_list)[0].takeoff_angle
        for d in distances])
    return np.gradient(takeoff_angles, distances)[1]


def geom_spread_teleseismic(
        angular_distance, source_depth_in_km, station_depth_in_km, phase):
    """
    Calculate geometrical spreading coefficient for teleseismic body waves.

    Implements eq (4) in Okal (1992) for a spherically symmetric Earth.
    This equations is derived from the conservation of the kinetic energy flux
    along a ray tube between the source and the receiver.

    :param angular_distance: Angular distance (degrees).
    :type angular_distance: float
    :param source_depth_in_km: Source depth (km).
    :type source_depth_in_km: float
    :param station_depth_in_km: Station depth (km).
    :type station_depth_in_km: float
    :param phase: Phase type (``'P'`` or ``'S'``).
    :type phase: str
    :return: Geometrical spreading correction (in m)
    :rtype: float
    """
    if phase == 'P':
        v_source = _get_medium_property_from_taup(source_depth_in_km, 'vp')
        v_station = _get_medium_property_from_taup(station_depth_in_km, 'vp')
        phase_list = ['p', 'P', 'pP', 'sP']
    elif phase == 'S':
        v_source = _get_medium_property_from_taup(source_depth_in_km, 'vs')
        v_station = _get_medium_property_from_taup(station_depth_in_km, 'vs')
        phase_list = ['s', 'S', 'sS', 'pS']
    else:
        raise ValueError(f'Invalid phase: {phase}')
    rho_source = _get_medium_property_from_taup(source_depth_in_km, 'rho')
    rho_station = _get_medium_property_from_taup(station_depth_in_km, 'rho')
    delta = np.deg2rad(angular_distance)
    arrival = model.get_travel_times(
        source_depth_in_km, angular_distance, phase_list)[0]
    # pylint: disable=no-member
    takeoff_angle = np.deg2rad(arrival.takeoff_angle)
    incident_angle = np.deg2rad(arrival.incident_angle)
    # Cmpute the local derivative of the takeoff angle (dt) with respect to the
    # angular distance (dd, also called the aperture of the ray tube).
    # We use a finite difference approximation, and we increase the aperture
    # until we get a non-zero value.
    aperture = 0.1  # degrees
    for _ in range(10):
        dtdd = _compute_dtdd(
            angular_distance, aperture, source_depth_in_km, phase_list)
        if dtdd != 0:
            break
        aperture *= 2
    if dtdd == 0:
        raise ValueError(
            f'Unable to compute geometrical spreading coefficient for '
            f'{phase} wave at {angular_distance:.2f} degrees: dt/dd = 0')
    # we have now all the ingredients to calculate the spreading coefficient,
    # eq (4) in Okal, 1992
    spreading_coeff = (
        (rho_source * v_source) / (rho_station * v_station) *
        np.sin(takeoff_angle) / np.sin(delta) *
        1 / np.cos(incident_angle) *
        np.abs(dtdd)
    )**0.5
    earth_radius = 6371e3  # m
    # We return the inverse of Okal's coefficient, since we use it to correct
    # the amplitude
    return earth_radius / spreading_coeff
# -----------------------------------------------------------------------------


# SOURCE PARAMETERS -----------------------------------------------------------
def moment_to_mag(moment):
    """Convert moment to magnitude."""
    return (np.log10(moment) - 9.1) / 1.5


def mag_to_moment(magnitude):
    """Convert magnitude to moment."""
    return np.power(10, (1.5 * magnitude + 9.1))


def source_radius(fc_in_hz, vs_in_m_per_s, k_coeff=0.3724):
    """
    Compute source radius in meters.

    Madariaga (2009), doi:10.1007/978-1-4419-7695-6_22, eq. 31
    Kaneko and Shearer (2014), doi:10.1093/gji/ggu030, eq. 2, 15, 16
    """
    return k_coeff * vs_in_m_per_s / fc_in_hz


def static_stress_drop(Mo_in_N_m, ra_in_m):
    """
    Compute static stress drop in MPa.

    Madariaga (2009), doi:10.1007/978-1-4419-7695-6_22, eq. 27
    """
    return 7. / 16 * Mo_in_N_m / ra_in_m**3 * 1e-6


def quality_factor(travel_time_in_s, t_star_in_s):
    """Compute quality factor from travel time and t_star."""
    return np.inf if t_star_in_s == 0 else travel_time_in_s / t_star_in_s
# -----------------------------------------------------------------------------


# SIGNAL ANALYSIS -------------------------------------------------------------
def cosine_taper(signal, width, left_taper=False):
    """Apply a cosine taper to the signal."""
    # TODO: this taper looks more like a hanning...
    npts = len(signal)
    p = 2 * width
    tap = _cos_taper(npts, p)
    if left_taper:
        tap[npts // 2:] = 1.
    signal *= tap


# modified from: http://stackoverflow.com/q/5515720
def smooth(signal, window_len=11, window='hanning'):
    """Smooth the signal using a window with requested size."""
    if signal.ndim != 1:
        raise ValueError('smooth only accepts 1 dimension arrays.')
    if signal.size < window_len:
        raise ValueError('Input vector needs to be bigger than window size.')
    if window_len < 3:
        return signal
    if window not in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError("Window is one of 'flat', 'hanning', 'hamming', "
                         "'bartlett', 'blackman'")
    s = np.r_[
        2 * signal[0] - signal[window_len - 1::-1],
        signal,
        2 * signal[-1] - signal[-1:-window_len:-1]
    ]
    if window == 'flat':  # moving average
        w = np.ones(window_len, 'd')
    else:
        w = eval(f'np.{window}(window_len)')  # pylint: disable=eval-used
    y = np.convolve(w / w.sum(), s, mode='same')
    yy = y[window_len:-window_len + 1]
    # check if there are NaN values
    nanindexes = np.where(np.isnan(yy))
    yy[nanindexes] = signal[nanindexes]
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
    """Convert degrees to radians."""
    return math.pi * degrees / 180


def toDeg(radians):
    """Convert radians to degrees."""
    return 180 * radians / math.pi


def _validate_coord(coordinate, coordinate_name):
    try:
        coordinate = float(coordinate)
    except Exception as e:
        raise ValueError(f'Invalid {coordinate_name}: {coordinate}') from e
    return coordinate


def station_to_event_position(trace):
    """
    Compute station position with respect to the event, in terms of hypocentral
    distance (km), epicentral distance (km), great-circle distance (degrees),
    azimuth and back-azimuth.

    Values are stored in the trace stats dictionary.
    """
    coords = trace.stats.coords
    stla = _validate_coord(coords.latitude, 'station latitude')
    stlo = _validate_coord(coords.longitude, 'station longitude')
    stel = _validate_coord(coords.elevation, 'station elevation')
    hypo = trace.stats.event.hypocenter
    evla = _validate_coord(hypo.latitude.value_in_deg, 'event latitude')
    evlo = _validate_coord(hypo.longitude.value_in_deg, 'event longitude')
    evdp = _validate_coord(hypo.depth.value_in_km, 'event depth')
    epi_dist, az, baz = gps2dist_azimuth(evla, evlo, stla, stlo)
    epi_dist /= 1e3  # in km
    gcarc = kilometers2degrees(epi_dist)
    hypo_dist = math.sqrt(epi_dist**2 + (stel + evdp)**2)
    trace.stats.azimuth = az
    trace.stats.back_azimuth = baz
    trace.stats.epi_dist = epi_dist
    trace.stats.hypo_dist = hypo_dist
    trace.stats.gcarc = gcarc
# -----------------------------------------------------------------------------
