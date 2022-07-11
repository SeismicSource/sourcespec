# -*- coding: utf8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Station plotting routine.

:copyright:
    2018-2022 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
import os
import logging
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
matplotlib.rcParams['pdf.fonttype'] = 42  # to edit text in Illustrator
# Reduce logging level for Matplotlib to avoid DEBUG messages
mpl_logger = logging.getLogger('matplotlib')
mpl_logger.setLevel(logging.WARNING)


# TODO:
# add table with values


def _round(x, base=5):
    """Round to base."""
    if x <= base/2:
        return 0
    round_x = 0
    while round_x == 0:
        round_x = int(base * round(float(x)/base))
        base = base/2.
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

    newcmap = colors.LinearSegmentedColormap(name, cdict)
    return newcmap


def _plot_circles(ax, evlon, evlat, maxdist, ncircles=5):
    geodetic_transform = ccrs.PlateCarree()
    g = Geod(ellps='WGS84')
    step = _round(maxdist/ncircles)
    if step == 0:
        step = 1
    if maxdist < 1:  # 1 km
        ncircles = 5
        step = 0.2
    texts = []
    for dist in np.arange(step, maxdist+step, step):
        azimuths = np.arange(0, 360, 1)
        circle = np.array(
            [g.fwd(evlon, evlat, az, dist*1e3)[0:2] for az in azimuths]
        )
        p0 = circle[np.argmax(circle[:, 1])]
        ax.plot(circle[:, 0], circle[:, 1],
                color='#777777', linestyle='--',
                transform=geodetic_transform)
        if maxdist < 1:
            dist_text = '{} m'.format(int(dist*1000))
        else:
            dist_text = '{} km'.format(int(dist))
        t = ax.text(p0[0], p0[1], dist_text, size=8, weight='bold',
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


def _plot_hypo(ax, hypo):
    geodetic_transform = ccrs.PlateCarree()
    try:
        strike = hypo.strike
        dip = hypo.dip
        rake = hypo.rake
        xy = ax.projection.transform_point(
            hypo.longitude, hypo.latitude, geodetic_transform)
        meca = beach(
            (strike, dip, rake),
            xy=xy,
            width=30,
            linewidth=1,
            facecolor='k',
            zorder=10,
            axes=ax
        )
        ax.add_collection(meca)
    except AttributeError:
        # plot hypocenter as a star
        ax.plot(
            hypo.longitude, hypo.latitude, marker='*', markersize=20,
            markeredgewidth=1, markeredgecolor='white',
            color='k', transform=geodetic_transform,
            zorder=10
        )


def _make_basemap(config, maxdist):
    g = Geod(ellps='WGS84')
    hypo = config.hypo
    maxdiagonal = maxdist*(2**0.5)*1.10
    lonmax, latmax, _ = g.fwd(
        hypo.longitude, hypo.latitude, 45, maxdiagonal*1000.)
    lonmin = 2*hypo.longitude - lonmax
    latmin = 2*hypo.latitude - latmax
    tile_dir = 'maptiles'
    stamen_terrain = CachedTiler(cimgt.Stamen('terrain-background'), tile_dir)
    # Create a GeoAxes
    figsize = (7.5, 7.5)
    fig = plt.figure(figsize=figsize, dpi=200)
    ax = fig.add_subplot(111, projection=stamen_terrain.crs)
    # Add event information as a title
    textstr = 'evid: {} \nlon: {:.3f} lat: {:.3f} depth: {:.1f} km'
    textstr = textstr.format(
        hypo.evid, hypo.longitude, hypo.latitude, hypo.depth)
    try:
        textstr += ' time: {}'.format(
            hypo.origin_time.format_iris_web_service())
    except AttributeError:
        pass
    ax.text(0., 1.15, textstr, fontsize=10,
            ha='left', va='top', linespacing=1.5, transform=ax.transAxes)
    trans = ccrs.Geodetic()
    ax.set_extent([lonmin, lonmax, latmin, latmax], crs=trans)
    if config.plot_map_tiles_zoom_level:
        tile_zoom_level = config.plot_map_tiles_zoom_level
    else:
        if maxdiagonal <= 100:
            tile_zoom_level = 12
        else:
            tile_zoom_level = 8
        logger.info('Map zoom level autoset to: {}'.format(tile_zoom_level))
    # Add the image twice, so that the CachedTiler has time to cache
    # (avoids white tiles on map)
    ax.add_image(stamen_terrain, tile_zoom_level)
    ax.add_image(stamen_terrain, tile_zoom_level)
    ax.gridlines(draw_labels=True, color='#777777', linestyle='--')
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
        if maxdiagonal <= 100:
            coastline_resolution = 'h'
        else:
            coastline_resolution = 'i'
        logger.info(
            'Coastline resolution autoset to: {}'.format(
                inv_res_map[coastline_resolution])
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
    circle_texts = _plot_circles(ax, hypo.longitude, hypo.latitude, maxdist, 5)
    _plot_hypo(ax, hypo)
    return ax, circle_texts


def _contrast_color(color):
    """
    Return the best contrasting color, either black or white.

    Source: https://stackoverflow.com/a/3943023/2021880
    """
    R, G, B = [
        c/12.92 if c <= 0.03928 else ((c+0.055)/1.055)**2.4 for c in color[:3]]
    L = 0.2126*R + 0.7152*G + 0.0722*B  # luminance
    if L > 0.179:
        return 'black'
    else:
        return 'white'


def _plot_stations(config, lonlat_dist, st_ids, values, vmean, verr, vname):
    maxdist = np.max(lonlat_dist[:, 2])
    ax, circle_texts = _make_basemap(config, maxdist)

    if config.options.evname is not None:
        textstr = '{} — '.format(config.options.evname)
    else:
        textstr = ''
    if vname == 'mag':
        verr_minus = verr_plus = verr
        textstr += 'Mw {:.2f} ± {:.2f}'.format(vmean, verr)
    elif vname == 'fc':
        verr_minus, verr_plus = verr
        textstr += 'fc {:.3f} [- {:.3f}, + {:.3f}] Hz'.format(
            vmean, verr_minus, verr_plus)
    ax.text(0., 1.22, textstr, fontsize=14,
            ha='left', va='top', transform=ax.transAxes)

    if vname == 'mag':
        vmax = np.max(np.abs(values-vmean))
        vmin = -vmax
        vmax += vmean
        vmin += vmean
        if vmax == vmin:
            vmax = vmean+0.5
            vmin = vmean-0.5
        cbar_extend = 'neither'
        cmap = cm.Spectral_r
    elif vname == 'fc':
        vmin = np.min(values)
        vmax = np.max(values)
        cbar_extend = 'neither'
        # limit colorbar to ±3sigma
        if vmax > vmean+3*verr_plus:
            vmax = vmean+3*verr_plus
            cbar_extend = 'max'
        if vmin < vmean-3*verr_minus:
            vmin = vmean-3*verr_minus
            if cbar_extend == 'max':
                cbar_extend = 'both'
            else:
                cbar_extend = 'min'
        if vmax == vmin:
            vmax = vmean+1
            vmin = vmean-1
        midpoint = (vmean-vmin)/(vmax-vmin)
        cmap = _shiftedColorMap(cm.PRGn, midpoint=midpoint)
    # Corner case, when errorbars are larger than vmin or vmax
    if vmin > vmean-verr_minus:
        vmin = vmean-(verr_minus*1.1)
    if vmax < vmean+verr_plus:
        vmax = vmean+(verr_plus*1.1)
    norm = colors.Normalize(vmin=vmin, vmax=vmax)

    trans = ccrs.PlateCarree()
    lonlat = lonlat_dist[:, :2]
    ax.scatter(
        lonlat[:, 0], lonlat[:, 1],
        marker='^', s=100,
        color=cmap(norm(values)), edgecolor='k',
        zorder=99, transform=trans)
    if config.plot_station_names_on_map:
        texts = []
        for ll, id in zip(lonlat, st_ids):
            id = '.'.join(id.split('.')[:2])
            id = '  ' + id
            station_text_size = config.plot_station_text_size
            t = ax.text(
                ll[0], ll[1], id, size=station_text_size, weight='bold',
                va='center', zorder=999, transform=trans)
            t.set_path_effects([
                PathEffects.Stroke(linewidth=0.8, foreground='white'),
                PathEffects.Normal()
            ])
            texts.append(t)

    # Add a colorbar
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
    color = _contrast_color(cmap(norm(vmean-verr_minus)))
    cax.axhline(
        vmean-verr_minus, lw=linewidth, linestyle=linestyle, color=color)
    color = _contrast_color(cmap(norm(vmean+verr_plus)))
    cax.axhline(
        vmean+verr_plus, lw=linewidth, linestyle=linestyle, color=color)
    if vname == 'mag':
        cm_label = 'Magnitude'
    elif vname == 'fc':
        cm_label = 'Corner Frequency (Hz)'
    cax.set_ylabel(cm_label)
    # Add code and author information at the figure bottom
    textstr = 'SourceSpec v{} '.format(get_versions()['version'])
    textstr += '– {} {} '.format(
        config.end_of_run.strftime('%Y-%m-%d %H:%M:%S'),
        config.end_of_run_tz)
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
        textstr = '{}\n{} '.format(textstr, textstr2)
    cax.text(1., -0.1, textstr, fontsize=8, linespacing=1.5,
             ha='right', va='top', transform=cax.transAxes)

    if config.plot_station_names_on_map:
        # first adjust text labels relatively to each other
        adjust_text(texts, ax=ax, maxshift=1e3)
        # then, try to stay away from circle texts
        adjust_text(texts, add_objects=circle_texts, ax=ax, maxshift=1e3)

    evid = config.hypo.evid
    figfile_base = os.path.join(config.options.outdir, evid)
    figfile_base += '.map_{}.'.format(vname)
    fmt = config.plot_save_format
    if fmt == 'pdf_multipage':
        fmt = 'pdf'
    figfile = figfile_base + fmt
    if config.plot_show:
        plt.show()
    if config.plot_save:
        savefig(fig, figfile, fmt, quantize_colors=False, bbox_inches='tight')
        if vname == 'mag':
            logger.info('Station-magnitude map saved to: ' + figfile)
        elif vname == 'fc':
            logger.info('Station-corner_freq map saved to: ' + figfile)
        config.figures['station_maps'].append(figfile)


def _spread_overlapping_stations(lonlat_dist, min_dlonlat=1e-3, spread=0.03):
    dlonlat = np.diff(lonlat_dist[:, :2], axis=0)
    dlonlat = np.sum(dlonlat**2, axis=1)**(0.5)
    # find indexes of overlapping stations
    idx = np.where(dlonlat < min_dlonlat)[0]
    # find consecutive indexes (i.e., more than two stations overlapping)
    didx = np.diff(idx)
    idx_idx = np.where(didx == 1)[0]
    idx_consec = idx[idx_idx+1]
    # spread stations horizontally
    lonlat_dist[idx] += np.array((-spread, 0, 0))
    lonlat_dist[idx+1] += np.array((spread, 0, 0))
    # further spread vertically a third station
    lonlat_dist[idx_consec] += np.array((0, spread, 0))
    return lonlat_dist


def plot_stations(config, sourcepar):
    """Plot station map, color coded by magnitude or fc."""
    # Check config, if we need to plot at all
    if not config.plot_show and not config.plot_save:
        return
    stationpar = sourcepar.station_parameters
    st_ids = sorted(stationpar.keys())
    lonlat_dist = np.array([
        (stationpar[k]['lon'], stationpar[k]['lat'], stationpar[k]['epi_dist'])
        for k in st_ids])
    maxdist = np.max(lonlat_dist[:, 2])
    # empirical rule to define overlapping station spread
    spread = maxdist/4e3
    _spread_overlapping_stations(lonlat_dist, min_dlonlat=1e-3, spread=spread)
    mag = np.array([stationpar[k]['Mw'] for k in st_ids])
    magmean = sourcepar.means_weight['Mw']
    magerr = sourcepar.errors_weight['Mw']
    _plot_stations(config, lonlat_dist, st_ids, mag, magmean, magerr, 'mag')
    fc = np.array([stationpar[k]['fc'] for k in st_ids])
    fcmean = sourcepar.means_weight['fc']
    fcerr = sourcepar.errors_weight['fc']
    _plot_stations(config, lonlat_dist, st_ids, fc, fcmean, fcerr, 'fc')
