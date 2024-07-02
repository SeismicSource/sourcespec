# -*- coding: utf-8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Functions to compute geometrical spreading coefficients.

:copyright:
    2012-2026 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
import logging
import numpy as np
from obspy.taup import TauPyModel
from .ssp_util import MediumProperties
model = TauPyModel(model='iasp91')
v_model = model.model.s_mod.v_mod
logger = logging.getLogger(__name__.rsplit('.', maxsplit=1)[-1])


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


def geom_spread_r_power_n_segmented(hypo_dist_in_km, exponents,
                                    hinge_distances):
    """
    Geometrical spreading function defined as piecewise continuous powerlaw,
    as defined in Boore (2003), eq. 9

    :param hypo_dist_in_km: Hypocentral distance (km).
    :type hypo_dist_in_km: float
    :param exponents: Exponents for different powerlaw segments
    :type exponents: numpy.ndarray
    :param hinge_distances: Distances defining start of powerlaw segments
    :type hinge_distances: numpy.ndarray
    :return: Geometrical spreading correction (for distance in m)
    :rtype: float
    """
    if np.isscalar(hypo_dist_in_km):
        hypo_dist_in_km = np.array([hypo_dist_in_km], dtype='float')
        is_scalar = True
    else:
        hypo_dist_in_km = np.asarray(hypo_dist_in_km, dtype='float')
        is_scalar = False
    hinge_distances = np.asarray(hinge_distances)
    Rref = hinge_distances[0]
    exponents = -np.asarray(exponents)
    # Do not allow distances less than Rref
    hypo_dist_in_km = np.maximum(Rref, hypo_dist_in_km)
    Zhinges = (hinge_distances[:-1] / hinge_distances[1:]) ** exponents[:-1]
    Zhinges = np.cumprod(Zhinges)
    R0, p0 = hinge_distances[0], exponents[0]
    Z = (R0 / hypo_dist_in_km) ** p0
    for n in range(1, len(hinge_distances)):
        Rn, pn = hinge_distances[n], exponents[n]
        idxs = hypo_dist_in_km > Rn
        Z[idxs] = Zhinges[n-1] * ((Rn / hypo_dist_in_km[idxs]) ** pn)
    # Convert spreading correction to metric distance
    Z *= 1e3
    if is_scalar:
        Z = Z[0]
    return Z


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
    freqs = np.atleast_1d(freqs).astype(float)
    exponent = np.ones_like(freqs)
    low_freq = freqs <= 0.2
    mid_freq = np.logical_and(freqs > 0.2, freqs <= 0.25)
    high_freq = freqs >= 0.25
    exponent[low_freq] = 0.5
    exponent[mid_freq] = 0.5 + 2 * np.log10(5 * freqs[mid_freq])
    exponent[high_freq] = 0.7
    coeff = cutoff_dist * (dist / cutoff_dist)**exponent
    return coeff[0] if len(coeff) == 1 else coeff


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
    # Don't need to specify coordinates, since we use a spherically symmetric
    # Earth; don't need to specify a config object, since we use the global
    # model (iasp91)
    medium_properties_source = MediumProperties(
        0, 0, source_depth_in_km, None)
    medium_properties_station = MediumProperties(
        0, 0, station_depth_in_km, None)
    if phase == 'P':
        v_source = medium_properties_source.get_from_taup('vp')
        v_station = medium_properties_station.get_from_taup('vp')
        phase_list = ['p', 'P', 'pP', 'sP']
    elif phase == 'S':
        v_source = medium_properties_source.get_from_taup('vs')
        v_station = medium_properties_station.get_from_taup('vs')
        phase_list = ['s', 'S', 'sS', 'pS']
    else:
        raise ValueError(f'Invalid phase: {phase}')
    rho_source = medium_properties_source.get_from_taup('rho')
    rho_station = medium_properties_station.get_from_taup('rho')
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
    # As a first guess, we set the aperture to half the angular distance
    aperture = angular_distance/2
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


# ------------------------- Plotting functions -------------------------

LINEWIDTH = 2.5
ZORDER = 0


def _plot_r_power_n(ax, x_data, hypo_dists, model_params):
    default_params = {'exponent': 1}
    params = model_params or default_params
    coeff = geom_spread_r_power_n(hypo_dists, **params)
    label = f'r^{params["exponent"]}'
    ax.plot(
        x_data, coeff,
        label=label, linewidth=LINEWIDTH, zorder=ZORDER)
    print(f'Done with r_power_n, params: {params}')


def _plot_r_power_n_segmented(ax, x_data, hypo_dists, model_params):
    default_params = {
        'exponents': np.array((1, 0, 0.5)),
        'hinge_distances': np.array((1, 70, 130)),
    }
    params = model_params or default_params
    coeff = geom_spread_r_power_n_segmented(hypo_dists, **params)
    exponents_str = ', '.join(map(str, params['exponents']))
    hinges_str = ', '.join(map(str, params['hinge_distances']))
    label = (
        f'Segmented r^n (exponents: {exponents_str}, '
        f'hinges: {hinges_str} km)'
    )
    ax.plot(
        x_data, coeff,
        label=label, linewidth=LINEWIDTH, zorder=ZORDER)
    print(f'Done with r_power_n_segmented, params: {params}')


def _plot_boatwright(ax, x_data, hypo_dists, model_params):
    default_params = {
        'cutoff_dist_in_km': 200,
        'freqs': np.array((0.))
    }
    params = model_params or default_params
    # set params['freqs'] to 0, if not provided
    params['freqs'] = params.get('freqs', 0.)
    freqs = np.atleast_1d(params['freqs'])
    if len(freqs) > 1:
        raise ValueError(
            'Multiple frequencies not supported for plotting '
            'Boatwright model'
        )
    coeff = [
        geom_spread_boatwright(d, **params)
        for d in hypo_dists
    ]
    label = (
        f'Boatwright (cutoff: {params["cutoff_dist_in_km"]} km, '
        f'freqs: {params["freqs"]} Hz)'
    )
    ax.plot(
        x_data, coeff,
        label=label, linewidth=LINEWIDTH, zorder=ZORDER)
    print(f'Done with boatwright, params: {params}')


def _plot_teleseismic(ax, x_data, angular_dist, source_depth, model_params):
    default_params = {'phase': 'P', 'station_depth_in_km': 0}
    params = model_params or default_params
    # set params['station_depth_in_km'] to 0, if not provided
    params['station_depth_in_km'] = (
        params.get('station_depth_in_km', 0))
    gst = np.vectorize(geom_spread_teleseismic)
    coeff = gst(angular_dist, source_depth, **params)
    label = (
        f'Teleseismic {params["phase"]} '
        f'(station depth: {params["station_depth_in_km"]} km)'
    )
    ax.plot(
        x_data, coeff,
        label=label, linewidth=LINEWIDTH, zorder=ZORDER)
    print(f'Done with teleseismic, params: {params}')


def plot_geom_spread_models(source_depth, epi_dists, models=None,
                            plottype='linlog', xaxis='epi_dist',
                            return_fig=False):
    """
    Plot geometrical spreading coefficients for different models.

    :param source_depth: Source depth (km).
    :type source_depth: float
    :param epi_dists: Epicentral distances (km).
    :type epi_dists: numpy.ndarray

    :param models: List of tuples where each tuple contains:

        - The model name.
        - A dictionary of model parameters (or ``None`` if the model has no
          parameters or to use defaults).

        Valid model names are:

        - ``'r_power_n'``: rⁿ geometrical spreading.
        - ``'r_power_n_segmented'``: Piecewise continuous power law,
          as defined in Boore (2003), eq. 9.
        - ``'boatwright'``: Boatwright et al. (2002).
        - ``'teleseismic'``: Teleseismic (Okal, 1992).

        If ``None``, defaults to ``[('r_power_n', {'exponent': 1})]``.
    :type models: list[tuple[str, dict or None]]

    :param plottype: Type of plot scaling. Options:

        - ``'linlin'``
        - ``'linlog'``
        - ``'loglin'``
        - ``'loglog'``
    :type plottype: str

    :param xaxis: X-axis type. Options:

        - ``'epi_dist'``: Epicentral distance.
        - ``'hypo_dist'``: Hypocentral distance.
    :type xaxis: str

    :param return_fig: If True, return the figure instead of showing it.
    :type return_fig: bool

    :return: Matplotlib figure, if ``return_fig`` is True.
    :rtype: matplotlib.figure.Figure

    :examples:

    >>> import numpy as np
    >>> source_depth = 10
    >>> epi_dists = np.arange(10, 1000, 10)
    >>> models = [
    >>>     ('r_power_n', {'exponent': 1}),
    >>>     ('boatwright', {'cutoff_dist_in_km': 100}),
    >>>     ('boatwright', {'cutoff_dist_in_km': 200, 'freqs': 0}),
    >>>     ('boatwright', {'cutoff_dist_in_km': 200, 'freqs': 10}),
    >>>     ('teleseismic', {'phase': 'P'}),
    >>>     ('teleseismic', {'phase': 'S'}),
    >>> ]
    >>> plot_geom_spread_models(
    >>>     source_depth, epi_dists, models=models,
    >>>     plottype='linlog', xaxis='epi_dist')

    >>> import numpy as np
    >>> source_depth = 0
    >>> epi_dists = np.arange(1, 1000, 1)
    >>> models = [
    >>>     ('r_power_n', {'exponent': 1}),
    >>>     ('r_power_n_segmented', None),
    >>> ]
    >>> plot_geom_spread_models(
    >>>     source_depth, epi_dists, models=models,
    >>>     plottype='loglog', xaxis='hypo_dist')
    """
    # Lazy import matplotlib, since it's only needed for this function
    # pylint: disable=import-outside-toplevel
    import matplotlib.pyplot as plt

    # source depth must be a float
    source_depth = float(source_depth)
    # epi_dists must be a numpy array
    epi_dists = np.atleast_1d(epi_dists)
    hypo_dists = np.sqrt(epi_dists**2 + source_depth**2)
    angular_dist = epi_dists / 111.111

    fig, ax = plt.subplots(figsize=(10, 6))

    if xaxis == 'epi_dist':
        x_data = epi_dists
        x_label = 'Epicentral distance (km)'
    elif xaxis == 'hypo_dist':
        x_data = hypo_dists
        x_label = 'Hypocentral distance (km)'
    else:
        raise ValueError(f'Invalid xaxis: {xaxis}')
    if plottype == 'linlin':
        ax.set_xscale('linear')
        ax.set_yscale('linear')
    elif plottype == 'linlog':
        ax.set_xscale('linear')
        ax.set_yscale('log')
    elif plottype == 'loglin':
        ax.set_xscale('log')
        ax.set_yscale('linear')
    elif plottype == 'loglog':
        ax.set_xscale('log')
        ax.set_yscale('log')
    else:
        raise ValueError(f'Invalid plottype: {plottype}')

    # Set default model if none is provided
    if models is None:
        models = [('r_power_n', {'exponent': 1})]

    print('Computing and plotting geometrical spreading coefficients...')
    # pylint: disable=global-statement, invalid-name
    global ZORDER
    for model_name, model_params in models:
        ZORDER += 10
        if model_name == 'r_power_n':
            _plot_r_power_n(ax, x_data, hypo_dists, model_params)
        elif model_name == 'r_power_n_segmented':
            _plot_r_power_n_segmented(ax, x_data, hypo_dists, model_params)
        elif model_name == 'boatwright':
            _plot_boatwright(ax, x_data, hypo_dists, model_params)
        elif model_name == 'teleseismic':
            _plot_teleseismic(ax, x_data, angular_dist, source_depth,
                              model_params)
        else:
            raise ValueError(f'Invalid model name: {model_name}')

    ax.minorticks_on()
    ax.grid(True, which='major', linestyle='-', linewidth=0.75)
    ax.grid(True, which='minor', linestyle='--', linewidth=0.5)
    ax.legend()
    ax.set_xlabel(x_label)
    ax.set_ylabel('Geometrical spreading coefficient')
    ax.set_title(f'Source depth: {source_depth} km')

    if return_fig:
        return fig
    plt.show()
