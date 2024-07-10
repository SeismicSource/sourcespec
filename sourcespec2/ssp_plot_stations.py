# -*- coding: utf8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Station plotting routine.

:copyright:
    2018-2024 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
import os
import logging
import contextlib
import warnings
import numpy as np
import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader
import cartopy.feature as cfeature
from obspy.imaging.beachball import beach
from pyproj import Geod
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import colors
import matplotlib.patheffects as PathEffects
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from .config import config
from .adjustText import adjust_text
from .cached_tiler import CachedTiler
from .map_tiles import (
    EsriHillshade,
    EsriHillshadeDark,
    EsriOcean,
    EsriImagery,
    StamenTerrain,
)
from .savefig import savefig
from ._version import get_versions
logger = logging.getLogger(__name__.rsplit('.', maxsplit=1)[-1])
# Reduce logging level for Matplotlib to avoid DEBUG messages
mpl_logger = logging.getLogger('matplotlib')
mpl_logger.setLevel(logging.WARNING)
# The following code fails with Sphinx and with Shapely<1,8, so we need to
# ignore any exception
with contextlib.suppress(Exception):
    # Ignore Shapely deprecation warnings, since they depend on Cartopy
    from shapely.errors import ShapelyDeprecationWarning
    warnings.filterwarnings('ignore', category=ShapelyDeprecationWarning)
TILER = {
    'hillshade': EsriHillshade,
    'hillshade_dark': EsriHillshadeDark,
    'ocean': EsriOcean,
    'satellite': EsriImagery,
    'stamen_terrain': StamenTerrain,
}
ZORDER_BASEMAP = 1
ZORDER_TILES = 2
ZORDER_COASTLINES = 5
ZORDER_CIRCLES = 10
ZORDER_EPICENTER = 20
ZORDER_STATIONS = 20
ZORDER_STATION_TEXTS = 30


# TODO: add table with values


def _round_to_base(x, base):
    """Round to base."""
    sign_x = np.sign(x)
    x = np.abs(x)
    if x <= base / 2:
        return 0
    round_x = 0
    while round_x == 0:
        round_x = int(base * np.ceil(float(x) / base))
        base = base / 2.
    return round_x*sign_x


# Source: https://stackoverflow.com/a/20528097
def _shiftedColorMap(cmap, start=0, midpoint=0.5, stop=1.0, name='shifted'):
    """
    Offset the "center" of a colormap.

    Useful for data with a negative min and positive max and you want the
    middle of the colormap's dynamic range to be at zero.

    Input
    -----
      cmap : The matplotlib colormap to be altered
      start : Offset from lowest point in the colormap's range.
          Defaults to 0.0 (no lower offset). Should be between
          0.0 and `midpoint`.
      midpoint : The new center of the colormap. Defaults to
          0.5 (no shift). Should be between 0.0 and 1.0. In
          general, this should be  1 - vmax / (vmax + abs(vmin))
          For example if your data range from -15.0 to +5.0 and
          you want the center of the colormap at 0.0, `midpoint`
          should be set to  1 - 5/(5 + 15)) or 0.75
      stop : Offset from highest point in the colormap's range.
          Defaults to 1.0 (no upper offset). Should be between
          `midpoint` and 1.0.
    """
    cdict = {
        'red': [],
        'green': [],
        'blue': [],
        'alpha': []
    }

    # regular index to compute the colors
    reg_index = np.linspace(start, stop, 257)

    # shifted index to match the data
    shift_index = np.hstack([
        np.linspace(0.0, midpoint, 128, endpoint=False),
        np.linspace(midpoint, 1.0, 129, endpoint=True)
    ])

    for ri, si in zip(reg_index, shift_index):
        r, g, b, a = cmap(ri)

        cdict['red'].append((si, r, r))
        cdict['green'].append((si, g, g))
        cdict['blue'].append((si, b, b))
        cdict['alpha'].append((si, a, a))

    return colors.LinearSegmentedColormap(name, cdict)


def _get_circles_step(maxdist, ncircles):
    """
    Return the step for the circles to be plotted.

    :param maxdist: maximum distance from the epicenter.
    :type maxdist: float
    :param ncircles: number of circles to plot.
    :type ncircles: int

    :return: step for the circles.
    :rtype: float
    """
    step = _round_to_base(maxdist / ncircles, base=5)
    if step == 0:
        step = 1
    if maxdist < 1:  # 1 km
        step = 0.2
    return step


def _plot_circles(ax, evlo, evla, distances):
    """
    Plot circles around the epicenter.

    :param ax: axes to plot the circles.
    :type ax: matplotlib.axes.Axes
    :param evlo: epicenter longitude.
    :type evlo: float
    :param evla: epicenter latitude.
    :type evla: float
    :param distances: distances from the epicenter.
    :type distances: list of float

    :return: list of texts with circle labels.
    :rtype: list of matplotlib.text.Text
    """
    geodetic_transform = ccrs.PlateCarree()
    g = Geod(ellps='WGS84')
    texts = []
    maxdist = np.max(distances)
    for dist in distances:
        azimuths = np.arange(0, 360, 1)
        circle = np.array(
            [g.fwd(evlo, evla, az, dist * 1e3)[:2] for az in azimuths]
        )
        p0 = circle[np.argmax(circle[:, 1])]
        ax.plot(
            circle[:, 0], circle[:, 1],
            color='#777777', linestyle='--',
            transform=geodetic_transform,
            zorder=ZORDER_CIRCLES)
        dist_text = (
            f'{int(dist * 1000)} m' if maxdist < 1
            else
            f'{int(dist)} km'
        )
        t = ax.text(
            p0[0], p0[1], dist_text, size=8, weight='bold',
            verticalalignment='center',
            horizontalalignment='center',
            clip_on=True,
            transform=geodetic_transform, zorder=ZORDER_CIRCLES)
        t.set_path_effects([
            PathEffects.Stroke(linewidth=0.8, foreground='white'),
            PathEffects.Normal()
        ])
        texts.append(t)
    return texts


def _plot_epicenter_as_beachball(ax, event):
    """
    Plot the epicenter as a beachball.

    :param ax: axes to plot the beachball.
    :type ax: matplotlib.axes.Axes
    :param event: event object.
    :type event: :class:`~sourcespec.ssp_event.SSPEvent`
    """
    geodetic_transform = ccrs.PlateCarree()
    fm = event.focal_mechanism
    # TODO: draw full moment tensor, if available
    # The following three lines will raise an exception if the focal mechanism
    # is not valid (e.g., values are None)
    strike = float(fm.strike)
    dip = float(fm.dip)
    rake = float(fm.rake)
    hypo = event.hypocenter
    evlo = hypo.longitude.value_in_deg
    evla = hypo.latitude.value_in_deg
    xy = ax.projection.transform_point(evlo, evla, geodetic_transform)
    # compute beachball width from map extent
    # Note: in previous versionos of SourceSpec, we used the argument axes=ax
    # to let beach() adapth the width to the axes size. However, this does not
    # work anymore with recent versions of Matplotlib, since the
    # IdentityTransform() used by beach() does not provide anymore a
    # translate() method.
    # TODO: this should be fixed in ObsPy 1.4.1,
    # see https://github.com/obspy/obspy/issues/2887
    # Change this again when bumping dependencies to ObsPy >= 1.4.1
    extent = ax.get_extent()
    xlen = extent[1] - extent[0]
    ylen = extent[3] - extent[2]
    width = min(xlen, ylen) / 15
    meca = beach(
        (strike, dip, rake),
        xy=xy,
        width=width,
        linewidth=1,
        facecolor='k',
        zorder=ZORDER_EPICENTER,
    )
    ax.add_collection(meca)


def _plot_epicenter_as_star(ax, event):
    """
    Plot the epicenter as a star.

    :param ax: axes to plot the star.
    :type ax: matplotlib.axes.Axes
    :param event: event object.
    :type event: :class:`~sourcespec.ssp_event.SSPEvent`
    """
    geodetic_transform = ccrs.PlateCarree()
    hypo = event.hypocenter
    evlo = hypo.longitude.value_in_deg
    evla = hypo.latitude.value_in_deg
    ax.plot(
        evlo, evla, marker='*', markersize=20,
        markeredgewidth=1, markeredgecolor='white',
        color='k', transform=geodetic_transform,
        zorder=ZORDER_EPICENTER
    )


def _add_event_info(event, ax):
    """
    Add event information as plot title.

    :param event: event object.
    :type event: :class:`~sourcespec.ssp_event.SSPEvent`
    :param ax: axes to plot the event information.
    :type ax: matplotlib.axes.Axes
    """
    evid = event.event_id
    hypo = event.hypocenter
    evlo = hypo.longitude.value_in_deg
    evla = hypo.latitude.value_in_deg
    evdp = hypo.depth.value_in_km
    textstr = f'{event.name} — ' if event.name else ''
    textstr += (
        f'evid: {evid}\n'
        f'lon: {evlo:.3f} lat: {evla:.3f} depth: {evdp:.1f} km'
    )
    with contextlib.suppress(AttributeError):
        textstr += f' time: {hypo.origin_time.format_iris_web_service()}'
    ax.text(
        0., 1.05, textstr, fontsize=10,
        ha='left', va='top', linespacing=1.5, transform=ax.transAxes)


def _add_tiles(ax, tiler, alpha=1):
    """
    Add map tiles to basemap.

    :param ax: axes to plot the tiles.
    :type ax: matplotlib.axes.Axes
    :param tiler: tiler object.
    :type tiler: :class:`~sourcespec.cached_tiler.CachedTiler`
    :param alpha: transparency of the tiles.
    :type alpha: float
    """
    if config.plot_map_tiles_zoom_level:
        tile_zoom_level = config.plot_map_tiles_zoom_level
    else:
        tile_zoom_level = 12 if ax.maxdiagonal <= 100 else 8
        logger.info(f'Map zoom level autoset to: {tile_zoom_level}')
    while True:
        if tile_zoom_level == 0:
            logger.warning('No map tiles found. Map will be blank.')
            break
        ax.add_image(tiler, tile_zoom_level, alpha=alpha, zorder=ZORDER_TILES)
        try:
            # draw tiles to check if they exist
            ax.get_figure().canvas.draw()
            break
        except ValueError:
            logger.warning(
                f'No map tiles found for zoom level {tile_zoom_level}. '
                f'Trying zoom level {tile_zoom_level-1}')
            ax.img_factories = []
            tile_zoom_level -= 1


def _add_coastlines(ax):
    """
    Add coastlines and borders to basemap.

    :param ax: axes to plot the coastlines.
    :type ax: matplotlib.axes.Axes
    """
    if config.plot_coastline_resolution == 'no_coastline':
        return
    # add coastlines from GSHHS
    res_map = {
        'full': 'f',
        'high': 'h',
        'intermediate': 'i',
        'low': 'l',
        'crude': 'c'
    }
    inv_res_map = {v: k for k, v in res_map.items()}
    if ax.global_projection:
        coastline_resolution = 'c'
    elif config.plot_coastline_resolution:
        coastline_resolution = res_map[config.plot_coastline_resolution]
    else:
        coastline_resolution = 'h' if ax.maxdiagonal <= 100 else 'i'
        logger.info(
            'Coastline resolution autoset to: '
            f'{inv_res_map[coastline_resolution]}'
        )
    shpfile = shpreader.gshhs(coastline_resolution)
    shp = shpreader.Reader(shpfile)
    with warnings.catch_warnings():
        # silent a warning on:
        # "Shapefile shape has invalid polygon: no exterior rings found"
        warnings.simplefilter('ignore')
        ax.add_geometries(
            shp.geometries(), ccrs.PlateCarree(),
            edgecolor='black', facecolor='none',
            zorder=ZORDER_COASTLINES)
    ax.add_feature(
        cfeature.BORDERS, edgecolor='black', facecolor='none',
        zorder=ZORDER_COASTLINES)


def _add_gridlines_to_orthographic_axes(ax, bounding_box):
    """
    Add gridlines to global orthographic GeoAxes.

    :param ax: axes to plot the gridlines.
    :type ax: matplotlib.axes.Axes
    :param bounding_box: bounding box of the map.
    :type bounding_box: list of float

    .. note::
    We need to compute gridlines manually, since Cartopy has occasional
    bugs with gridlines in global projections.

    The error is the following:
        shapely.errors.GEOSException:
        IllegalArgumentException: point array must contain 0 or >1 elements
    """
    lonmin, lonmax, latmin, latmax = bounding_box
    lonmin_grid = _round_to_base(lonmin, base=10)
    lonmax_grid = _round_to_base(lonmax, base=10)
    latmin_grid = _round_to_base(latmin, base=10)
    latmax_grid = _round_to_base(latmax, base=10)
    # longitude gridlines, every 30 degrees
    xticks = np.arange(lonmin_grid, lonmax_grid + 1, 30)
    # latitude gridlines, every 20 degrees
    if latmin_grid < 0 < latmax_grid:
        # make sure 0 is in the yticks
        yticks_neg = np.arange(latmin_grid, 0, 20)
        yticks_pos = np.arange(0, latmax_grid + 1, 20)
        yticks = np.unique(np.concatenate((yticks_neg, yticks_pos)))
    else:
        yticks = np.arange(latmin_grid, latmax_grid + 1, 20)
    ax.gridlines(
        draw_labels=True, color='#777777', linestyle='--',
        xlocs=xticks, ylocs=yticks)


def _make_geoaxes_mercator(fig, buonding_box, maxdiagonal):
    """
    Create a GeoAxes with Mercator projection and optionally add map tiles.

    :param fig: figure to add the axes.
    :type fig: matplotlib.figure.Figure
    :param buonding_box: bounding box of the map.
    :type buonding_box: list of float
    :param maxdiagonal: maximum diagonal of the map.
    :type maxdiagonal: float

    :return: axes with Mercator projection.
    :rtype: matplotlib.axes.Axes
    """
    land_10m = cfeature.NaturalEarthFeature(
        'physical', 'land', '10m',
        edgecolor='face',
        facecolor=cfeature.COLORS['land'])
    ocean_10m = cfeature.NaturalEarthFeature(
        'physical', 'ocean', '10m',
        edgecolor='face',
        facecolor=cfeature.COLORS['water'])
    map_style = config.plot_map_style
    api_key = config.plot_map_api_key
    if map_style == 'no_basemap':
        ax = fig.add_subplot(111, projection=ccrs.Mercator())
        ax.add_feature(land_10m, zorder=ZORDER_BASEMAP)
        ax.add_feature(ocean_10m, zorder=ZORDER_BASEMAP)
    else:
        tile_dir = 'maptiles'
        tiler = CachedTiler(
            TILER[map_style](apikey=api_key),
            tile_dir
        )
        ax = fig.add_subplot(111, projection=tiler.crs)
    ax.set_extent(buonding_box, crs=ccrs.Geodetic())
    ax.global_projection = False
    ax.maxdiagonal = maxdiagonal
    if map_style != 'no_basemap':
        if map_style in ['hillshade', 'hillshade_dark']:
            # add a sea mask to the hillshade map
            ax.add_feature(ocean_10m, zorder=ZORDER_TILES+1)
        _add_tiles(ax, tiler)
    if map_style in ['hillshade', 'hillshade_dark', 'ocean', 'satellite']:
        ax.attribution_text = 'Map powered by Esri and Natural Earth'
    elif map_style == 'stamen_terrain':
        ax.attribution_text = 'Map powered by Stamen Design and Natural Earth'
    else:
        ax.attribution_text = 'Map powered by Natural Earth'
    ax.gridlines(draw_labels=True, color='#777777', linestyle='--')
    return ax


def _make_geoaxes_orthographic(fig, evlo, evla, bounding_box):
    """
    Create a GeoAxes with global Orthographic projection.

    The basemap is the Earth stock image from Cartopy.

    :param fig: figure to add the axes.
    :type fig: matplotlib.figure.Figure
    :param evlo: epicenter longitude.
    :type evlo: float
    :param evla: epicenter latitude.
    :type evla: float
    :param bounding_box: bounding box of the map.
    :type bounding_box: list of float

    :return: axes with Orthographic projection.
    :rtype: matplotlib.axes.Axes
    """
    _projection = ccrs.Orthographic(
        central_longitude=evlo, central_latitude=evla)
    ax = fig.add_subplot(111, projection=_projection)
    ax.global_projection = True
    ax.stock_img()
    ax.attribution_text = 'Map powered by Natural Earth'
    _add_gridlines_to_orthographic_axes(ax, bounding_box)
    return ax


def _make_basemap(maxdist):
    """
    Create basemap with tiles, coastlines, hypocenter
    and distance circles.

    :param maxdist: maximum distance from the epicenter.
    :type maxdist: float

    :return: ax0 (invisible axis for title and footer),
        ax (GeoAxes with basemap),
        circle_texts (list of texts with circle labels).
    :rtype: tuple of matplotlib.axes.Axes, matplotlib.axes.Axes,
        list of matplotlib.text.Text
    """
    g = Geod(ellps='WGS84')
    event = config.event
    evlo = event.hypocenter.longitude.value_in_deg
    evla = event.hypocenter.latitude.value_in_deg
    ncircles = 5
    circles_step = _get_circles_step(maxdist, ncircles)
    # compute bounding box, adding 0.2 to ncircles to have some padding
    maxdist = (ncircles + 0.2) * circles_step
    lonmin = g.fwd(evlo, evla, 270, maxdist * 1e3)[0]
    lonmax = g.fwd(evlo, evla, 90, maxdist * 1e3)[0]
    latmin = g.fwd(evlo, evla, 180, maxdist * 1e3)[1]
    latmax = g.fwd(evlo, evla, 0, maxdist * 1e3)[1]
    bounding_box = [lonmin, lonmax, latmin, latmax]
    # maxdiagonal is the distance between the two corners of the map
    # (i.e., the diagonal of the bounding box)
    maxdiagonal = g.inv(lonmin, latmin, lonmax, latmax)[2] / 1e3
    # Reduce dpi for vector formats, since the only raster are the maptiles
    if config.plot_save_format in ['pdf', 'pdf_multipage', 'svg']:
        figsize = (8.5, 8.5)
        dpi = 72
    else:
        figsize = (7.5, 7.5)
        dpi = 200
    fig = plt.figure(figsize=figsize, dpi=dpi)
    # Create an invisible axis and use it for title and footer
    ax0 = fig.add_subplot(111, label='ax0')
    ax0.set_axis_off()
    _add_event_info(config.event, ax0)
    if maxdist < 3000:  # km
        ax = _make_geoaxes_mercator(fig, bounding_box, maxdiagonal)
    else:
        ax = _make_geoaxes_orthographic(fig, evlo, evla, bounding_box)
    _add_coastlines(ax)
    circles_distances = np.arange(1, ncircles + 1) * circles_step
    circle_texts = _plot_circles(ax, evlo, evla, circles_distances)
    try:
        _plot_epicenter_as_beachball(ax, event)
    except Exception:
        _plot_epicenter_as_star(ax, event)
    return ax0, ax, circle_texts


def _add_footer(ax, attribution_text=None):
    """
    Add code and author information at the figure footer.

    :param ax: axes to plot the footer.
    :type ax: matplotlib.axes.Axes
    :param attribution_text: additional attribution text.
    :type attribution_text: str
    """
    textstr = (
        f'SourceSpec v{get_versions()["version"]} '
        f'- {config.end_of_run.strftime("%Y-%m-%d %H:%M:%S")} '
        f'{config.end_of_run_tz} '
    )
    textstr2 = ''
    if config.author_name is not None:
        textstr2 += config.author_name
    elif config.author_email is not None:
        textstr2 += config.author_email
    if config.agency_short_name is not None:
        if textstr2 != '':
            textstr2 += ' - '
        textstr2 += config.agency_short_name
    elif config.agency_full_name is not None:
        if textstr2 != '':
            textstr2 += ' - '
        textstr2 += config.agency_full_name
    if attribution_text is not None:
        if textstr2 != '':
            textstr2 += ' - '
        textstr2 += attribution_text
    if textstr2 != '':
        textstr = f'{textstr}\n{textstr2} '
    ax.text(
        1.12, 0, textstr, fontsize=8, linespacing=1.5,
        ha='right', va='top', transform=ax.transAxes)


def _add_main_title(ax, vname, vmean, verr):
    """
    Add to the figure the main title with the value and its error.

    :param ax: axes to plot the title.
    :type ax: matplotlib.axes.Axes
    :param vname: name of the value.
    :type vname: str
    :param vmean: mean value.
    :type vmean: float
    :param verr: error value.
    :type verr: float or tuple of float
    """
    verr_minus, verr_plus = _get_verr_minus_plus(verr)
    if vname == 'fc':
        textstr = (
            f'fc {vmean:.3f} ± {verr_minus:.3f} Hz'
            if np.isclose(verr_minus, verr_plus, atol=1e-3)
            else f'fc {vmean:.3f} [- {verr_minus:.3f}, + {verr_plus:.3f}] Hz'
        )
    elif vname == 'mag':
        textstr = (
            f'Mw {vmean:.2f} ± {verr_minus:.2f}'
            if np.isclose(verr_minus, verr_plus, atol=1e-2)
            else f'Mw {vmean:.2f} [- {verr_minus:.2f}, + {verr_plus:.2f}]'
        )
    ax.text(
        0., 1.09, textstr, fontsize=14,
        ha='left', va='top', transform=ax.transAxes)


def _contrast_color(color):
    """
    Return the best contrasting color, either black or white.

    Source: https://stackoverflow.com/a/3943023/2021880

    :param color: color in RGB format.
    :type color: tuple of float

    :return: 'black' or 'white'.
    :rtype: str
    """
    R, G, B = [
        c / 12.92 if c <= 0.03928
        else ((c + 0.055) / 1.055)**2.4
        for c in color[:3]]
    L = 0.2126 * R + 0.7152 * G + 0.0722 * B  # luminance
    return 'black' if L > 0.179 else 'white'


def _get_verr_minus_plus(verr):
    """
    Return the minus and plus error values for the given error value.

    :param verr: error value.
    :type verr: float or tuple of float

    :return: minus and plus error values.
    :rtype: tuple of float
    """
    if isinstance(verr, tuple):
        verr_minus, verr_plus = verr
    else:
        verr_minus = verr_plus = verr
    return verr_minus, verr_plus


def _get_cmap_and_norm(values, outliers, vname, vmean, verr):
    """
    Return the colormap and normalization for the given values.

    :param values: values to plot.
    :type values: numpy.ndarray
    :param outliers: boolean array with the outliers.
    :type outliers: numpy.ndarray
    :param vname: name of the value.
    :type vname: str
    :param vmean: mean value.
    :type vmean: float
    :param verr: error value.
    :type verr: float or tuple of float

    :return: colormap, normalization and colorbar extension.
    :rtype: matplotlib.colors.Colormap, matplotlib.colors.Normalize, str
    """
    verr_minus, verr_plus = _get_verr_minus_plus(verr)
    values_no_outliers = values[~outliers]
    if vname == 'mag':
        vmax = np.max(np.abs(values_no_outliers - vmean))
        vmin = -vmax
        vmax += vmean
        vmin += vmean
        if np.isclose(vmax, vmin):
            vmax = vmean + 0.1
            vmin = vmean - 0.1
        cbar_extend = 'neither'
        cmap = cm.Spectral_r  # pylint: disable=no-member
    elif vname == 'fc':
        vmin = np.min(values_no_outliers)
        vmax = np.max(values_no_outliers)
        cbar_extend = 'neither'
        # limit colorbar to ±3sigma
        if vmax > vmean + 3 * verr_plus:
            vmax = vmean + 3 * verr_plus
            cbar_extend = 'max'
        if vmin < vmean - 3 * verr_minus:
            vmin = vmean - 3 * verr_minus
            cbar_extend = 'both' if cbar_extend == 'max' else 'min'
        if vmax == vmin:
            vmax = vmean + 1
            vmin = vmean - 1
        midpoint = (vmean - vmin) / (vmax - vmin)
        cmap = _shiftedColorMap(
            cm.PRGn, midpoint=midpoint)  # pylint: disable=no-member
    # Corner case, when errorbars are larger than vmin or vmax
    if vmin > vmean - verr_minus:
        vmin = vmean - (verr_minus * 1.1)
    if vmax < vmean + verr_plus:
        vmax = vmean + (verr_plus * 1.1)
    norm = colors.Normalize(vmin=vmin, vmax=vmax)
    return cmap, norm, cbar_extend


def _plot_stations_scatter(ax, lonlat_dist, st_ids, values, cmap, norm):
    """
    Plot the stations as scatter points on the map.

    :param ax: axes to plot the stations.
    :type ax: matplotlib.axes.Axes
    :param lonlat_dist: array with the station coordinates and distances.
    :type lonlat_dist: numpy.ndarray
    :param st_ids: list with the station IDs.
    :type st_ids: list of str
    :param values: values to plot.
    :type values: numpy.ndarray
    :param cmap: colormap.
    :type cmap: matplotlib.colors.Colormap
    :param norm: normalization.
    :type norm: matplotlib.colors.Normalize

    :return: list of texts with station labels.
    :rtype: list of matplotlib.text.Text
    """
    trans = ccrs.PlateCarree()
    lonlat = lonlat_dist[:, :2]
    ax.scatter(
        lonlat[:, 0], lonlat[:, 1],
        marker='^', s=100,
        color=cmap(norm(values)), edgecolor='k',
        zorder=ZORDER_STATIONS, transform=trans)
    texts = []
    if config.plot_station_names_on_map:
        for _lonlat, _statid in zip(lonlat, st_ids):
            _statid = '.'.join(_statid.split('.')[:2])
            _statid = f'  {_statid}'
            station_text_size = config.plot_station_text_size
            t = ax.text(
                _lonlat[0], _lonlat[1], _statid,
                size=station_text_size, weight='bold',
                va='center', zorder=ZORDER_STATION_TEXTS, transform=trans)
            t.set_path_effects([
                PathEffects.Stroke(linewidth=0.8, foreground='white'),
                PathEffects.Normal()
            ])
            texts.append(t)
    return texts


def _adjust_text_labels(station_texts, circle_texts, ax):
    """
    Adjust the text labels so that they do not overlap.

    :param station_texts: list of texts with station labels.
    :type station_texts: list of matplotlib.text.Text
    :param circle_texts: list of texts with circle labels.
    :type circle_texts: list of matplotlib.text.Text
    :param ax: axes to plot the text labels.
    :type ax: matplotlib.axes.Axes
    """
    if not station_texts:
        return
    # store original text positions and texts
    text_pos = [(t.get_position(), t.get_text()) for t in station_texts]
    # compute mean position
    x_pos_mean = np.mean([p[0][0] for p in text_pos])
    y_pos_mean = np.mean([p[0][1] for p in text_pos])
    # first adjust text labels relatively to each other
    adjust_text(station_texts, ax=ax, maxshift=1e3)
    # then, try to stay away from circle texts
    adjust_text(station_texts, add_objects=circle_texts, ax=ax, maxshift=1e3)
    # check if some text labels are too far away from the mean position
    # (i.e., bug in adjust_text) and move them back to their original position
    for t in station_texts:
        txt = t.get_text()
        x_pos, y_pos = t.get_position()
        delta_x = np.abs(x_pos - x_pos_mean)
        delta_y = np.abs(y_pos - y_pos_mean)
        if delta_x > 100 or delta_y > 100:
            x_pos, y_pos = [tp[0] for tp in text_pos if tp[1] == txt][0]
            t.set_position((x_pos, y_pos))


def _add_colorbar(ax, cmap, norm, vmean, verr, vname, cbar_extend):
    """
    Add a colorbar to the given axes.

    :param ax: axes to plot the colorbar.
    :type ax: matplotlib.axes.Axes
    :param cmap: colormap.
    :type cmap: matplotlib.colors.Colormap
    :param norm: normalization.
    :type norm: matplotlib.colors.Normalize
    :param vmean: mean value.
    :type vmean: float
    :param verr: error value.
    :type verr: float or tuple of float
    :param vname: name of the value.
    :type vname: str
    :param cbar_extend: colorbar extension.
    :type cbar_extend: str
    """
    ax_divider = make_axes_locatable(ax)
    cax = ax_divider.append_axes(
        'right', size='6%', pad='15%', axes_class=plt.Axes)
    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    fig = ax.get_figure()
    fig.colorbar(sm, cax=cax, extend=cbar_extend)
    cax.get_yaxis().set_visible(True)
    cax.axhline(vmean, lw=2, color='black')
    linestyle = (0, (2, 1))
    linewidth = 1.5
    verr_minus, verr_plus = _get_verr_minus_plus(verr)
    color = _contrast_color(cmap(norm(vmean - verr_minus)))
    cax.axhline(
        vmean - verr_minus, lw=linewidth, linestyle=linestyle, color=color)
    color = _contrast_color(cmap(norm(vmean + verr_plus)))
    cax.axhline(
        vmean + verr_plus, lw=linewidth, linestyle=linestyle, color=color)
    if vname == 'fc':
        cm_label = 'Corner Frequency (Hz)'
    elif vname == 'mag':
        cm_label = 'Magnitude'
    cax.set_ylabel(cm_label)


def _savefig(fig, vname):
    """
    Save the figure to a file.

    :param fig: figure to save.
    :type fig: matplotlib.figure.Figure
    :param vname: name of the value.
    :type vname: str
    """
    evid = config.event.event_id
    figfile_base = os.path.join(config.options.outdir, evid)
    figfile_base += f'.map_{vname}.'
    fmt = config.plot_save_format
    if fmt == 'pdf_multipage':
        fmt = 'pdf'
    figfile = figfile_base + fmt
    if config.plot_show:
        plt.show()
    if config.plot_save:
        savefig(fig, figfile, fmt, quantize_colors=False, bbox_inches='tight')
        if vname == 'mag':
            logger.info(f'Station-magnitude map saved to: {figfile}')
        elif vname == 'fc':
            logger.info(f'Station-corner_freq map saved to: {figfile}')
        config.figures['station_maps'].append(figfile)


def _make_station_map(
        lonlat_dist, st_ids, values, outliers, vmean, verr, vname):
    """
    Make a map of stations with the given values.

    :param lonlat_dist: array with the station coordinates and distances.
    :type lonlat_dist: numpy.ndarray
    :param st_ids: list with the station IDs.
    :type st_ids: list of str
    :param values: values to plot.
    :type values: numpy.ndarray
    :param outliers: boolean array with the outliers.
    :type outliers: numpy.ndarray
    :param vmean: mean value.
    :type vmean: float
    :param verr: error value.
    :type verr: float or tuple of float
    :param vname: name of the value.
    :type vname: str
    """
    maxdist = np.max(lonlat_dist[:, 2])
    ax0, ax, circle_texts = _make_basemap(maxdist)
    _add_main_title(ax0, vname, vmean, verr)
    cmap, norm, cbar_extend = _get_cmap_and_norm(
        values, outliers, vname, vmean, verr)
    station_texts = _plot_stations_scatter(
        ax, lonlat_dist, st_ids, values, cmap, norm)
    _add_colorbar(ax, cmap, norm, vmean, verr, vname, cbar_extend)
    _add_footer(ax0, ax.attribution_text)
    _adjust_text_labels(station_texts, circle_texts, ax)
    _savefig(ax0.get_figure(), vname)


def _spread_overlapping_stations(lonlat_dist, min_dlonlat=1e-3, spread=0.03):
    """
    Spread overlapping stations horizontally and vertically.

    :param lonlat_dist: array with the station coordinates and distances.
    :type lonlat_dist: numpy.ndarray
    :param min_dlonlat: minimum distance between stations.
    :type min_dlonlat: float
    :param spread: spread distance.
    :type spread: float

    :return: array with the station coordinates and distances.
    :rtype: numpy.ndarray
    """
    dlonlat = np.diff(lonlat_dist[:, :2], axis=0)
    dlonlat = np.sum(dlonlat**2, axis=1)**(0.5)
    # find indexes of overlapping stations
    idx = np.where(dlonlat < min_dlonlat)[0]
    # find consecutive indexes (i.e., more than two stations overlapping)
    didx = np.diff(idx)
    idx_idx = np.where(didx == 1)[0]
    idx_consec = idx[idx_idx + 1]
    # spread stations horizontally
    lonlat_dist[idx] += np.array((-spread, 0, 0))
    lonlat_dist[idx + 1] += np.array((spread, 0, 0))
    # further spread vertically a third station
    lonlat_dist[idx_consec] += np.array((0, spread, 0))
    return lonlat_dist


def plot_stations(sspec_output):
    """
    Plot station map, color coded by magnitude or fc.

    :param sspec_output: SourceSpecOutput object
    :type sspec_output: ssp_data_types.SourceSpecOutput
    """
    # Check config, if we need to plot at all
    if not config.plot_show and not config.plot_save:
        return
    matplotlib.rcParams['pdf.fonttype'] = 42  # to edit text in Illustrator
    stationpar = sspec_output.station_parameters
    st_ids = sorted(stationpar.keys())
    lonlat_dist = np.array([
        (stationpar[k]['longitude'], stationpar[k]['latitude'],
         stationpar[k]['epi_dist_in_km'])
        for k in st_ids])
    maxdist = np.max(lonlat_dist[:, 2])
    # empirical rule to define overlapping station spread
    spread = maxdist / 4e3
    _spread_overlapping_stations(lonlat_dist, min_dlonlat=1e-3, spread=spread)
    summary_values = sspec_output.reference_values()
    summary_uncertainties = sspec_output.reference_uncertainties()
    mag = np.array([stationpar[k]['Mw'].value for k in st_ids])
    mag_outliers = np.array([stationpar[k]['Mw'].outlier for k in st_ids])
    summary_mag = summary_values['Mw']
    summary_mag_err = summary_uncertainties['Mw']
    _make_station_map(
        lonlat_dist, st_ids,
        mag, mag_outliers, summary_mag, summary_mag_err, 'mag')
    fc = np.array([stationpar[k]['fc'].value for k in st_ids])
    fc_outliers = np.array([stationpar[k]['fc'].outlier for k in st_ids])
    summary_fc = summary_values['fc']
    summary_fc_err = summary_uncertainties['fc']
    _make_station_map(
        lonlat_dist, st_ids,
        fc, fc_outliers, summary_fc, summary_fc_err, 'fc')
