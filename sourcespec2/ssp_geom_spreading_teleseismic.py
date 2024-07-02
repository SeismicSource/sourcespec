# -*- coding: utf-8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Teleseismic geometrical spreading coefficient.

Implements the geometrical spreading coefficient for teleseismic body waves
in a spherically symmetric Earth, following Okal (1992), eq. 4:

    G(Δ) = a / g(Δ)

where ``a`` is the Earth radius (m) and ``g(Δ)`` is defined as:

    g(Δ) = sqrt(
        (ρ_h c_h) / (ρ_r c_r)   # density × velocity at the source
        × sin(i_h) / sin(Δ)     # takeoff angle and angular distance
        × 1 / cos(i_r)          # incidence angle at the receiver
        × |d i_h / d Δ|         # aperture of the ray tube
    )

with:
    - ``Δ``: great-circle distance between source and receiver (radians)
    - ``ρ_h``, ``c_h``: density and P- or S-wave velocity at the
      hypocenter
    - ``ρ_r``, ``c_r``: density and P- or S-wave velocity at the receiver
    - ``i_h``: takeoff angle at the hypocenter, measured from the downward
      vertical (0° is straight down, 180° is straight up)
    - ``i_r``: incidence angle at the receiver, measured from the downward
      vertical
    - ``d i_h / d Δ``: variation of the takeoff angle within a ray tube of
      angular width ``Δ``

``G`` is expressed in meters and is the amplitude correction to apply: we
correct amplitude and not energy, hence the square root in ``g``.

The spreading curve ``G(Δ)`` is precomputed once for a given
(source depth, wave phase) pair on a fine angular-distance grid and then
interpolated at the angular distances requested by the caller.  The
precomputation is cached, so that it is performed only once per run,
regardless of the number of stations.

:copyright:
    2012-2026 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
import logging

import numpy as np
from obspy.taup import TauPyModel
from scipy.ndimage import median_filter
from scipy.signal import savgol_filter

logger = logging.getLogger(__name__.rsplit('.', maxsplit=1)[-1])

MODEL = TauPyModel(model='iasp91')
EARTH_RADIUS_M = 6371e3  # m

# Phase names passed to TauP for each wave type.  The lowercase phase is the
# upgoing ray (the one leaving the source upward), the uppercase phase is the
# downgoing ray.
_PHASE_LISTS = {'P': ['p', 'P'], 'S': ['s', 'S']}

# Module-level cache of the precomputed spreading curves, keyed by the
# parameters of the precomputation.
_spreading_cache = {}


def _medium_properties(depth_in_km, phase):
    """
    Return density (kg/m³) and wave velocity (m/s) at ``depth_in_km``.

    Values are taken from the iasp91 model used by TauP.  The density is
    converted from g/cm³ to kg/m³ and the velocity from km/s to m/s.
    Note that the Okal formula only uses the ratio
    ``(ρ_h c_h) / (ρ_r c_r)``, so the absolute units do not matter, as
    long as they are consistent between hypocenter and receiver.  SI units
    are used here for clarity.

    :param depth_in_km: Depth (km).
    :type depth_in_km: float
    :param phase: Wave type (``'P'`` or ``'S'``).
    :type phase: str
    :return: (density in kg/m³, velocity in m/s).
    :rtype: tuple(float, float)
    """
    # The TauP velocity model (s_mod.v_mod) contains all three
    # properties: 'p' velocity, 's' velocity and density ('r'), so it is
    # used for both wave types.
    v_mod = MODEL.model.s_mod.v_mod
    prop = 's' if phase == 'S' else 'p'
    # Avoid negative depths (e.g. stations above sea level).
    depth_in_km = max(depth_in_km, 1e-2)
    rho = v_mod.evaluate_above(depth_in_km, 'r')[0] * 1e3
    vel = v_mod.evaluate_above(depth_in_km, prop)[0] * 1e3
    return rho, vel


def _select_direct_branch(distances, depth_in_km, phase_list):
    """
    Select the direct ray branch at each distance.

    The direct branch is the arrival with the largest takeoff angle, i.e.
    the shallowest ray connecting source and receiver (the one leaving the
    source closest to the horizontal/upward direction).

    :param distances: Angular distances (degrees).
    :type distances: numpy.ndarray
    :param depth_in_km: Source depth (km).
    :type depth_in_km: float
    :param phase_list: Phase names passed to TauP.
    :type phase_list: list[str]
    :return: Dictionary with ``distance``, ``takeoff_angle`` and
        ``incident_angle`` arrays.
    :rtype: dict
    """
    dists = []
    takeoffs = []
    incidents = []
    for dist in distances:
        arrivals = list(MODEL.get_travel_times(
            depth_in_km, dist, phase_list))
        if not arrivals:
            continue
        arrival = max(arrivals, key=lambda a: a.takeoff_angle)
        dists.append(dist)
        takeoffs.append(arrival.takeoff_angle)
        incidents.append(getattr(arrival, 'incident_angle', np.nan))
    return {
        'distance': np.array(dists),
        'takeoff_angle': np.array(takeoffs),
        'incident_angle': np.array(incidents),
    }


def _takeoff_derivative(takeoff_angle_deg, delta_deg, window, polyorder):
    """
    Compute ``d i_h / d Δ``, the derivative of the takeoff angle with
    respect to the angular distance.

    This derivative is dimensionless: both the takeoff angle and the angular
    distance are angles, so deg/deg == rad/rad and no unit conversion is
    needed.  It is computed with a Savitzky-Golay first derivative, which is
    robust to the numerical noise appearing at branch switches and caustics.
    When the window is too large for the available data, a simple gradient
    is used as a fallback.

    :param takeoff_angle_deg: Takeoff angles (degrees).
    :type takeoff_angle_deg: numpy.ndarray
    :param delta_deg: Angular distances (degrees), uniform grid.
    :type delta_deg: numpy.ndarray
    :param window: Savitzky-Golay window length.
    :type window: int
    :param polyorder: Savitzky-Golay polynomial order.
    :type polyorder: int
    :return: ``d i_h / d Δ`` (dimensionless).
    :rtype: numpy.ndarray
    """
    # Savitzky-Golay requires an odd window.
    window = window if window % 2 else window + 1
    if takeoff_angle_deg.size >= window > polyorder:
        step = np.median(np.diff(delta_deg))
        return savgol_filter(
            takeoff_angle_deg, window_length=window, polyorder=polyorder,
            deriv=1, delta=step)
    return np.gradient(takeoff_angle_deg, delta_deg)


def _smooth_log(x, window_length, polyorder, despike_size):
    """
    Smooth a positive series in log space.

    Optionally median-despikes first, then applies a Savitzky-Golay filter.
    NaN values (and non-positive ones) are interpolated in log space before
    smoothing and restored afterwards.

    :param x: Input series.
    :type x: numpy.ndarray
    :param window_length: Savitzky-Golay window length.
    :type window_length: int
    :param polyorder: Savitzky-Golay polynomial order.
    :type polyorder: int
    :param despike_size: Median filter size for despiking; 1 means no
        despiking.
    :type despike_size: int
    :return: Smoothed series, same shape as ``x``.
    :rtype: numpy.ndarray
    """
    x = np.asarray(x, dtype=float)
    valid = np.isfinite(x) & (x > 0)
    # Savitzky-Golay requires an odd window larger than the polynomial order.
    window = window_length if window_length % 2 else window_length + 1
    window = max(window, polyorder + 1)
    if np.sum(valid) < window:
        return x
    logx = np.log(np.where(valid, x, np.nan))
    finite = np.isfinite(logx)
    idx = np.arange(logx.size)
    # Fill NaN values (e.g. distances where the phase does not exist) by
    # linear interpolation in log space.
    logx_filled = np.interp(idx, idx[finite], logx[finite])
    if despike_size > 1:
        logx_filled = median_filter(logx_filled, size=despike_size)
    logx_smooth = savgol_filter(
        logx_filled, window_length=window, polyorder=polyorder)
    out = np.exp(logx_smooth)
    out[~valid] = np.nan
    return out


class TeleseismicSpreading:
    """
    Precomputed teleseismic geometrical spreading curve.

    On instantiation, the spreading curve ``G(Δ)`` is computed on a fine
    angular-distance grid for a given source depth and wave phase, using the
    direct (shallowest) ray branch at each distance.  The curve is despiked
    and smoothed in log space before being stored.

    The instance is callable: given an angular distance (or an array of
    distances) and a station depth, it returns the interpolated spreading
    coefficient, corrected for the density and velocity at the station
    depth (see :meth:`__call__`).
    """

    def __init__(self, source_depth_in_km, phase, min_delta_deg=0.5,
                 max_delta_deg=100.0, step_deg=0.1, smooth_window=5,
                 smooth_poly=2, despike_size=10):
        """
        Precompute the spreading curve.

        :param source_depth_in_km: Source depth (km).
        :type source_depth_in_km: float
        :param phase: Wave type (``'P'`` or ``'S'``).
        :type phase: str
        :param min_delta_deg: Minimum angular distance of the grid (degrees).
        :type min_delta_deg: float
        :param max_delta_deg: Maximum angular distance of the grid (degrees).
        :type max_delta_deg: float
        :param step_deg: Grid step (degrees).
        :type step_deg: float
        :param smooth_window: Savitzky-Golay window for the derivative and
            the smoothing.
        :type smooth_window: int
        :param smooth_poly: Savitzky-Golay polynomial order.
        :type smooth_poly: int
        :param despike_size: Median filter size for despiking; 1 disables.
        :type despike_size: int
        """
        self.source_depth_in_km = source_depth_in_km
        self.phase = phase
        self.smooth_window = smooth_window
        self.smooth_poly = smooth_poly
        self.despike_size = despike_size
        logger.info(
            'Precomputing teleseismic spreading curve '
            f'for {phase} wave, source depth {source_depth_in_km:g} km')
        self._precompute(min_delta_deg, max_delta_deg, step_deg)
        logger.info(
            f'Precomputed teleseismic spreading curve '
            f'({self.delta_grid.size} points, '
            f'{self.delta_grid[0]:g}°–{self.delta_grid[-1]:g}°)')

    def _precompute(self, min_delta_deg, max_delta_deg, step_deg):
        """
        Compute the spreading curve on the angular-distance grid.

        The computation follows these steps:

        1. Select the direct branch at each grid point.
        2. Compute ``d i_h / d Δ`` with a Savitzky-Golay derivative.
        3. Compute the Okal (1992) ``g(Δ)`` and ``G(Δ) = a / g(Δ)``,
           using the density and velocity at the source and at a reference
           receiver located at the surface.
        4. Despike and smooth ``G`` in log space.

        The reference receiver (surface) is used so that the curve can be
        precomputed independently of the station depth: the actual station
        depth only changes the ``ρ_r c_r`` factor, which is applied at call
        time as a simple multiplicative correction (see :meth:`__call__`).
        """
        delta_grid = np.arange(min_delta_deg, max_delta_deg, step_deg)
        branch = _select_direct_branch(
            delta_grid, self.source_depth_in_km, _PHASE_LISTS[self.phase])
        self.delta_grid = branch['distance']
        if self.delta_grid.size == 0:
            logger.warning(
                'No arrivals found: teleseismic spreading curve is empty')
            self.g_ref = np.array([])
            self._ref_medium_factor = np.nan
            return
        takeoff = branch['takeoff_angle']
        incident = branch['incident_angle']
        # Derivative of the takeoff angle with respect to the angular
        # distance: d i_h / d Δ (dimensionless, see _takeoff_derivative).
        dtdd = _takeoff_derivative(
            takeoff, self.delta_grid, self.smooth_window, self.smooth_poly)
        # Density and velocity at the source.
        rho_h, c_h = _medium_properties(self.source_depth_in_km, self.phase)
        # Density and velocity at the reference receiver (surface).
        rho_r0, c_r0 = _medium_properties(0.0, self.phase)
        self._ref_medium_factor = rho_r0 * c_r0
        # Okal (1992), eq. 4:
        #   g(Δ)² = (ρ_h c_h) / (ρ_r c_r) * sin(i_h) / sin(Δ)
        #           * 1 / cos(i_r) * |d i_h / d Δ|
        #   G(Δ) = a / g(Δ)
        delta_rad = np.deg2rad(self.delta_grid)
        with np.errstate(divide='ignore', invalid='ignore'):
            g = np.sqrt(
                (rho_h * c_h) / self._ref_medium_factor
                * np.sin(np.deg2rad(takeoff)) / np.sin(delta_rad)
                * 1.0 / np.cos(np.deg2rad(incident))
                * np.abs(dtdd)
            )
            g_ref = EARTH_RADIUS_M / g
        g_ref = np.where(np.isfinite(g_ref) & (g_ref > 0), g_ref, np.nan)
        # Despike and smooth in log space.
        self.g_ref = _smooth_log(
            g_ref, self.smooth_window, self.smooth_poly, self.despike_size)

    def __call__(self, angular_distance, station_depth_in_km=0.0):
        """
        Return the spreading coefficient at ``angular_distance``.

        The precomputed curve (valid for a receiver at the surface) is
        interpolated at the requested distance(s) and corrected for the
        density and velocity at the station depth:

            G(Δ) = G_ref(Δ) * sqrt( (ρ_r c_r) / (ρ_r0 c_r0) )

        where ``ρ_r0 c_r0`` is the factor at the surface.  For a station at
        the surface the correction is exactly 1.

        :param angular_distance: Angular distance (degrees), scalar or array.
        :type angular_distance: float or numpy.ndarray
        :param station_depth_in_km: Station depth (km).
        :type station_depth_in_km: float
        :return: Spreading coefficient (m).
        :rtype: float or numpy.ndarray
        """
        if self.g_ref.size == 0:
            return np.nan
        angular_distance = np.asarray(angular_distance, dtype=float)
        is_scalar = angular_distance.ndim == 0
        if is_scalar:
            angular_distance = angular_distance.reshape(1)
        # np.interp clips the distances outside the grid to the edge values.
        g_ref = np.interp(
            angular_distance, self.delta_grid, self.g_ref)
        # Correction for the actual station depth.
        rho_r, c_r = _medium_properties(station_depth_in_km, self.phase)
        station_factor = np.sqrt(rho_r * c_r / self._ref_medium_factor)
        g = g_ref * station_factor
        return g[0] if is_scalar else g


def geom_spreading_teleseismic(
        angular_distance, source_depth_in_km, station_depth_in_km, phase,
        min_delta_deg=0.5, max_delta_deg=100.0, step_deg=0.1,
        smooth_window=5, smooth_poly=2, despike_size=10):
    """
    Calculate the geometrical spreading coefficient for teleseismic body
    waves, following Okal (1992), eq. 4.

    The spreading curve is precomputed once for the given
    (source depth, phase) pair and cached; subsequent calls only interpolate
    the curve at the requested distance.  This is much faster than computing
    a ray trace for every station.

    :param angular_distance: Angular distance (degrees), scalar or array.
    :type angular_distance: float or numpy.ndarray
    :param source_depth_in_km: Source depth (km).
    :type source_depth_in_km: float
    :param station_depth_in_km: Station depth (km).
    :type station_depth_in_km: float
    :param phase: Wave type (``'P'`` or ``'S'``).
    :type phase: str
    :param min_delta_deg: Minimum angular distance of the grid (degrees).
    :type min_delta_deg: float
    :param max_delta_deg: Maximum angular distance of the grid (degrees).
    :type max_delta_deg: float
    :param step_deg: Grid step (degrees).
    :type step_deg: float
    :param smooth_window: Savitzky-Golay window for the derivative and the
        smoothing.
    :type smooth_window: int
    :param smooth_poly: Savitzky-Golay polynomial order.
    :type smooth_poly: int
    :param despike_size: Median filter size for despiking; 1 disables.
    :type despike_size: int
    :return: Geometrical spreading correction (in m).
    :rtype: float or numpy.ndarray
    """
    if phase not in _PHASE_LISTS:
        raise ValueError(f'Invalid phase: {phase}')
    key = (
        source_depth_in_km, phase, min_delta_deg, max_delta_deg, step_deg,
        smooth_window, smooth_poly, despike_size)
    spreader = _spreading_cache.get(key)
    if spreader is None:
        spreader = TeleseismicSpreading(
            source_depth_in_km, phase,
            min_delta_deg=min_delta_deg, max_delta_deg=max_delta_deg,
            step_deg=step_deg, smooth_window=smooth_window,
            smooth_poly=smooth_poly, despike_size=despike_size)
        _spreading_cache[key] = spreader
    return spreader(angular_distance, station_depth_in_km)
