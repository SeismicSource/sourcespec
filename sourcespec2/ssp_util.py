# -*- coding: utf-8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Utility functions for sourcespec.

:copyright:
    2012-2025 Claudio Satriano <satriano@ipgp.fr>
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
from .setup import config
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
class MediumProperties():
    """
    Class to retrieve medium properties from config.

    :param lon: Longitude (degrees).
    :type lon: float
    :param lat: Latitude (degrees).
    :type lat: float
    :param depth_in_km: Depth (km).
    :type depth_in_km: float
    """

    def __init__(self, lon, lat, depth_in_km):
        self.lon = lon
        self.lat = lat
        self.depth_in_km = depth_in_km

    def get_from_config_param_source(self, mproperty):
        """
        Get medium property at the source from config parameter.

        :param mproperty: Property to be retrieved
            (``'vp'``, ``'vs'`` or ``'rho'``).
        :type mproperty: str
        :return: Property value.
        :rtype: float
        """
        if mproperty not in ['vp', 'vs', 'rho']:
            raise ValueError(f'Invalid property: {mproperty}')
        values = config[f'{mproperty}_source']
        if values is None:
            return None
        if config.layer_top_depths is None:
            return values[0]
        values = np.array(values)
        depths = np.array(config.layer_top_depths)
        try:
            # find the last value that is smaller than the source depth
            value = values[depths <= self.depth_in_km][-1]
        except IndexError:
            value = values[0]
        return value

    def get_from_config_param_station(self, mproperty):
        """
        Get medium property at the station from config parameter.

        :param mproperty: Property to be retrieved
            (``'vp'``, ``'vs'`` or ``'rho'``).
        :type mproperty: str
        :return: Property value.
        :rtype: float
        """
        if mproperty not in ['vp', 'vs', 'rho']:
            raise ValueError(f'Invalid property: {mproperty}')
        value = config[f'{mproperty}_stations']
        if value is None:
            value = self.get_from_config_param_source(mproperty)
        return value

    def get_vel_from_NLL(self, wave):
        """
        Get velocity from NLL model.

        :param wave: Wave type (``'P'`` or ``'S'``).
        :type wave: str
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
            raise FileNotFoundError(
                f'Unable to find model file {grdfile}') from e
        grd = NLLGrid(grdfile)
        x, y = grd.project(self.lon, self.lat)
        value = grd.get_value(x, y, self.depth_in_km)
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

    def get_from_taup(self, mproperty):
        """
        Get medium property (P- or S-wave velocity, density) at a given depth
        from taup model.

        :param mproperty: Property to be retrieved
            (``'vp'``, ``'vs'`` or ``'rho'``).
        :type mproperty: str
        :return: Property value.
        :rtype: float
        """
        # avoid negative depths
        depth_in_km = max(self.depth_in_km, 1e-2)
        try:
            prop = {'vp': 'p', 'vs': 's', 'rho': 'r'}[mproperty]
        except KeyError as e:
            raise ValueError(f'Invalid property: {mproperty}') from e
        value = v_model.evaluate_above(depth_in_km, prop)[0]
        if mproperty == 'rho':
            value *= 1e3  # convert g/cm**3 to kg/m**3
        return value

    def get(self, mproperty, where):
        """
        Get medium property (P- or S-wave velocity, density) at a given point
        from NLL grid, config or taup model.

        :param mproperty: Property to be retrieved
            (``'vp'``, ``'vs'`` or ``'rho'``).
        :type mproperty: str
        :param where: Where to retrieve medium property
            (``'source'`` or ``'stations'``).
        :type where: str
        :return: Property value.
        :rtype: float
        """
        if mproperty not in ['vp', 'vs', 'rho']:
            raise ValueError(f'Invalid property: {mproperty}')
        if where == 'source':
            value = self.get_from_config_param_source(
                mproperty)
        elif where == 'stations':
            value = self.get_from_config_param_station(
                mproperty)
        else:
            raise ValueError(f'Invalid location: {where}')
        if (
            config.NLL_model_dir is not None and
            mproperty in ['vp', 'vs']
        ):
            wave = 'P' if mproperty == 'vp' else 'S'
            try:
                value = self.get_vel_from_NLL(wave)
            except Exception as msg:
                logger.warning(msg)
                logger.warning(f'Using {wave} velocity from config')
        if value is None:
            value = self.get_from_taup(mproperty)
            logger.info(
                f'Using {mproperty} from global model (iasp91)')
        return float(value)

    def to_string(self, mproperty, value):
        """
        Return a string with the property name and value.

        :param mproperty: Property name. Must contain one of the following:
                        ``'vp'``, ``'vs'``, ``'rho'``, ``'depth'``
        :type mproperty: str
        :param value: Property value.
        :type value: float
        :return: Property string.
        :rtype: str
        """
        if 'vp' in mproperty or 'vs' in mproperty:
            value = round(value, 2)
            unit = 'km/s'
        elif 'rho' in mproperty:
            value = round(value, 1)
            unit = 'kg/m3'
        elif 'depth' in mproperty:
            value = round(value, 1)
            unit = 'km'
        else:
            raise ValueError(f'Invalid property: {mproperty}')
        # we need to use float(value) to be sure that is_integer() is defined
        if float(value).is_integer():
            value = int(value)
        return f'{mproperty}: {value} {unit}'
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
