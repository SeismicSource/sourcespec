# -*- coding: utf8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Station plotting routine.

:copyright:
    2018-2023 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
import os
import logging
import contextlib
import numpy as np
import warnings
import cartopy.crs as ccrs
import cartopy.io.img_tiles as cimgt
import cartopy.io.shapereader as shpreader
import cartopy.feature as cfeature
from obspy.imaging.beachball import beach
from sourcespec.adjustText import adjust_text
from pyproj import Geod
from sourcespec.cached_tiler import CachedTiler
from sourcespec.savefig import savefig
from sourcespec._version import get_versions
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib.patheffects as PathEffects
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
logger = logging.getLogger(__name__.split('.')[-1])
# Reduce logging level for Matplotlib to avoid DEBUG messages
mpl_logger = logging.getLogger('matplotlib')
mpl_logger.setLevel(logging.WARNING)
# The following code fails with Sphinx and with Shapely<1,8, so we need to
# ignore any exception
with contextlib.suppress(Exception):
    # Ignore Shapely deprecation warnings, since they depend on Cartopy
    from shapely.errors import ShapelyDeprecationWarning
    warnings.filterwarnings('ignore', category=ShapelyDeprecationWarning)


# TODO:
# add table with values


def _round(x, base=5):
    """Round to base."""
    if x <= base / 2:
        return 0
    round_x = 0
    while round_x == 0:
        round_x = int(base * round(float(x) / base))
        base = base / 2.
    return round_x


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


def _plot_circles(ax, evlo, evla, maxdist, ncircles=5):
    geodetic_transform = ccrs.PlateCarree()
    g = Geod(ellps='WGS84')
    step = _round(maxdist / ncircles)
    if step == 0:
        step = 1
    if maxdist < 1:  # 1 km
        ncircles = 5
        step = 0.2
    texts = []
    for dist in np.arange(step, maxdist + step, step):
        azimuths = np.arange(0, 360, 1)
        circle = np.array(
            [g.fwd(evlo, evla, az, dist * 1e3)[:2] for az in azimuths]
        )
        p0 = circle[np.argmax(circle[:, 1])]
        ax.plot(
            circle[:, 0], circle[:, 1],
            color='#777777', linestyle='--',
            transform=geodetic_transform)
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
            transform=geodetic_transform, zorder=10)
        t.set_path_effects([
            PathEffects.Stroke(linewidth=0.8, foreground='white'),
            PathEffects.Normal()
        ])
        texts.append(t)
    return texts


def _plot_epicenter_as_beachball(ax, event):
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
        zorder=10,
    )
    ax.add_collection(meca)


def _plot_epicenter_as_star(ax, event):
    geodetic_transform = ccrs.PlateCarree()
    hypo = event.hypocenter
    evlo = hypo.longitude.value_in_deg
    evla = hypo.latitude.value_in_deg
    ax.plot(
        evlo, evla, marker='*', markersize=20,
        markeredgewidth=1, markeredgecolor='white',
        color='k', transform=geodetic_transform,
        zorder=10
    )


def _add_event_info(event, ax):
    """Add event information as plot title."""
    evid = event.event_id
    hypo = event.hypocenter
    evlo = hypo.longitude.value_in_deg
    evla = hypo.latitude.value_in_deg
    evdp = hypo.depth.value_in_km
    textstr = (
        f'evid: {evid} \nlon: {evlo:.3f} lat: {evla:.3f} '
        f'depth: {evdp:.1f} km'
    )
    with contextlib.suppress(AttributeError):
        textstr += f' time: {hypo.origin_time.format_iris_web_service()}'
    ax.text(
        0., 1.05, textstr, fontsize=10,
        ha='left', va='top', linespacing=1.5, transform=ax.transAxes)


def _add_tiles(config, ax, stamen_terrain):
    """Add map tiles to basemap."""
    if config.plot_map_tiles_zoom_level:
        tile_zoom_level = config.plot_map_tiles_zoom_level
    else:
        tile_zoom_level = 12 if ax.maxdiagonal <= 100 else 8
        logger.info(f'Map zoom level autoset to: {tile_zoom_level}')
    while True:
        if tile_zoom_level == 0:
            logger.warning('No map tiles found. Map will be blank.')
            break
        ax.add_image(stamen_terrain, tile_zoom_level)
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


def _add_coastlines(config, ax):
    """Add coastlines and borders to basemap."""
    # add coastlines from GSHHS
    res_map = {
        'full': 'f',
        'high': 'h',
        'intermediate': 'i',
        'low': 'l',
        'crude': 'c'
    }
    inv_res_map = {v: k for k, v in res_map.items()}
    if config.plot_coastline_resolution:
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
            edgecolor='black', facecolor='none')
    ax.add_feature(cfeature.BORDERS, edgecolor='black', facecolor='none')


def _make_basemap(config, maxdist):
    """Create basemap with tiles, coastlines, hypocenter
    and distance circles."""
    g = Geod(ellps='WGS84')
    # increase bounding box for large maxdist,
    # to account for deformation in Mercator projection
    mult = 1.1 if maxdist < 500 else 1.5
    maxdiagonal = maxdist * (2 ** 0.5) * mult
    event = config.event
    hypo = event.hypocenter
    evlo = hypo.longitude.value_in_deg
    evla = hypo.latitude.value_in_deg
    lonmax, latmax, _ = g.fwd(evlo, evla, 45, maxdiagonal * 1000.)
    lonmin = 2 * evlo - lonmax
    latmin = 2 * evla - latmax
    tile_dir = 'maptiles'
    stamen_terrain = CachedTiler(cimgt.Stamen('terrain-background'), tile_dir)
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
    ax = fig.add_subplot(111, projection=stamen_terrain.crs)
    ax.set_extent([lonmin, lonmax, latmin, latmax], crs=ccrs.Geodetic())
    ax.maxdiagonal = maxdiagonal
    _add_tiles(config, ax, stamen_terrain)
    _add_coastlines(config, ax)
    ax.gridlines(draw_labels=True, color='#777777', linestyle='--')
    circle_texts = _plot_circles(ax, evlo, evla, maxdist, 5)
    try:
        _plot_epicenter_as_beachball(ax, event)
    except Exception:
        _plot_epicenter_as_star(ax, event)
    return ax0, ax, circle_texts


def _add_footer(config, ax):
    """
    Add code and author information at the figure footer.
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
    if textstr2 != '':
        textstr = f'{textstr}\n{textstr2} '
    ax.text(
        1.12, 0, textstr, fontsize=8, linespacing=1.5,
        ha='right', va='top', transform=ax.transAxes)


def _add_main_title(config, ax, vname, vmean, verr):
    """
    Add to the figure the main title with the value and its error.
    """
    verr_minus, verr_plus = _get_verr_minus_plus(verr)
    textstr = f'{config.options.evname} — ' if config.options.evname else ''
    if vname == 'fc':
        textstr += (
            f'fc {vmean:.3f} ± {verr_minus:.3f} Hz'
            if np.isclose(verr_minus, verr_plus, atol=1e-3)
            else f'fc {vmean:.3f} [- {verr_minus:.3f}, + {verr_plus:.3f}] Hz'
        )
    elif vname == 'mag':
        textstr += (
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
    """
    if isinstance(verr, tuple):
        verr_minus, verr_plus = verr
    else:
        verr_minus = verr_plus = verr
    return verr_minus, verr_plus


def _get_cmap_and_norm(values, outliers, vname, vmean, verr):
    """
    Return the colormap and normalization for the given values.
    """
    verr_minus, verr_plus = _get_verr_minus_plus(verr)
    values_no_outliers = values[~outliers]
    if vname == 'mag':
        vmax = np.max(np.abs(values_no_outliers - vmean))
        vmin = -vmax
        vmax += vmean
        vmin += vmean
        if vmax == vmin:
            vmax = vmean + 0.5
            vmin = vmean - 0.5
        cbar_extend = 'neither'
        cmap = cm.Spectral_r
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
        cmap = _shiftedColorMap(cm.PRGn, midpoint=midpoint)
    # Corner case, when errorbars are larger than vmin or vmax
    if vmin > vmean - verr_minus:
        vmin = vmean - (verr_minus * 1.1)
    if vmax < vmean + verr_plus:
        vmax = vmean + (verr_plus * 1.1)
    norm = colors.Normalize(vmin=vmin, vmax=vmax)
    return cmap, norm, cbar_extend


def _plot_stations_scatter(
        config, ax, lonlat_dist, st_ids, values, cmap, norm, circle_texts):
    """
    Plot the stations as scatter points on the map.
    """
    trans = ccrs.PlateCarree()
    lonlat = lonlat_dist[:, :2]
    ax.scatter(
        lonlat[:, 0], lonlat[:, 1],
        marker='^', s=100,
        color=cmap(norm(values)), edgecolor='k',
        zorder=99, transform=trans)
    if config.plot_station_names_on_map:
        texts = []
        for _lonlat, _statid in zip(lonlat, st_ids):
            _statid = '.'.join(_statid.split('.')[:2])
            _statid = f'  {_statid}'
            station_text_size = config.plot_station_text_size
            t = ax.text(
                _lonlat[0], _lonlat[1], _statid,
                size=station_text_size, weight='bold',
                va='center', zorder=999, transform=trans)
            t.set_path_effects([
                PathEffects.Stroke(linewidth=0.8, foreground='white'),
                PathEffects.Normal()
            ])
            texts.append(t)
        # first adjust text labels relatively to each other
        adjust_text(texts, ax=ax, maxshift=1e3)
        # then, try to stay away from circle texts
        adjust_text(texts, add_objects=circle_texts, ax=ax, maxshift=1e3)


def _add_colorbar(
        ax, cmap, norm, vmean, verr, vname, cbar_extend):
    """
    Add a colorbar to the given axes.
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


def _savefig(config, fig, vname):
    """
    Save the figure to a file.
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
        config, lonlat_dist, st_ids, values, outliers, vmean, verr, vname):
    """
    Make a map of stations with the given values.
    """
    maxdist = np.max(lonlat_dist[:, 2])
    ax0, ax, circle_texts = _make_basemap(config, maxdist)
    _add_main_title(config, ax0, vname, vmean, verr)
    cmap, norm, cbar_extend = _get_cmap_and_norm(
        values, outliers, vname, vmean, verr)
    _plot_stations_scatter(
        config, ax, lonlat_dist, st_ids, values, cmap, norm, circle_texts)
    _add_colorbar(ax, cmap, norm, vmean, verr, vname, cbar_extend)
    _add_footer(config, ax0)
    _savefig(config, ax0.get_figure(), vname)


def _spread_overlapping_stations(lonlat_dist, min_dlonlat=1e-3, spread=0.03):
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


def plot_stations(config, sspec_output):
    """
    Plot station map, color coded by magnitude or fc.

    :param config: Configuration object
    :type config: config.Config
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
        config, lonlat_dist, st_ids,
        mag, mag_outliers, summary_mag, summary_mag_err, 'mag')
    fc = np.array([stationpar[k]['fc'].value for k in st_ids])
    fc_outliers = np.array([stationpar[k]['fc'].outlier for k in st_ids])
    summary_fc = summary_values['fc']
    summary_fc_err = summary_uncertainties['fc']
    _make_station_map(
        config, lonlat_dist, st_ids,
        fc, fc_outliers, summary_fc, summary_fc_err, 'fc')
