# -*- coding: utf8 -*-
"""
Station plotting routine.

:copyright:
    2018-2019 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement, Version 2.1
    (http://www.cecill.info/index.en.html)
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import os
import logging
import numpy as np
import cartopy.crs as ccrs
import cartopy.io.img_tiles as cimgt
import cartopy.feature as cfeature
from pyproj import Geod
from sourcespec.cached_tiler import CachedTiler
from sourcespec.ssp_version import get_git_version

# TODO:
# add table with values


def _import_mpl(config):
    import matplotlib
    matplotlib.rcParams['pdf.fonttype'] = 42  # to edit text in Illustrator
    # Reduce logging level for Matplotlib to avoid DEBUG messages
    mpl_logger = logging.getLogger('matplotlib')
    mpl_logger.setLevel(logging.WARNING)
    if config.PLOT_SHOW:
        global plt
        import matplotlib.pyplot as plt
    else:
        global Figure
        global FigureCanvasAgg
        from matplotlib.figure import Figure
        from matplotlib.backends.backend_agg import FigureCanvasAgg
    global cm
    global colors
    global transforms
    global patches
    global PathEffects
    global make_axes_locatable
    import matplotlib.cm as cm
    import matplotlib.colors as colors
    import matplotlib.transforms as transforms
    import matplotlib.patches as patches
    import matplotlib.patheffects as PathEffects
    from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable


def _round(x, base=5):
    """Round to base."""
    round_x = 0
    while round_x == 0:
        round_x = int(base * round(float(x)/base))
        base = base/2.
    return round_x


def _plot_circles(ax, evlon, evlat, maxdist, ncircles=5):
    geodetic_transform = ccrs.PlateCarree()
    g = Geod(ellps='WGS84')
    ax.plot(evlon, evlat, marker='*', markersize=20,
            markeredgewidth=1, markeredgecolor='white',
            color='k', transform=geodetic_transform,
            zorder=10)
    step = _round(maxdist/ncircles)
    for dist in np.arange(step, maxdist+step, step):
        azimuths = np.arange(0, 360, 1)
        circle = np.array(
            [g.fwd(evlon, evlat, az, dist*1e3)[0:2] for az in azimuths]
        )
        p0 = circle[np.argmax(circle[:, 1])]
        ax.plot(circle[:, 0], circle[:, 1],
                color='#777777', linestyle='--',
                transform=geodetic_transform)
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


def plot_stations(config, sourcepar, st=None):
    """Plot station map, color coded by magnitude or fc."""
    _import_mpl(config)
    st_ids = [
        k for k in sorted(sourcepar.keys())
        if k not in ['means', 'errors', 'means_weight', 'errors_weight']]
    lonlat = np.array([
        (sourcepar[k]['lon'], sourcepar[k]['lat']) for k in st_ids])
    epi_dist = np.array([sourcepar[k]['epi_dist'] for k in st_ids])
    if st is not None:
        lonlat_all = np.array([
            (tr.stats.coords.longitude, tr.stats.coords.latitude)
            for tr in st
        ])
        epi_dist = np.array([tr.stats.epi_dist for tr in st])
    maxdist = np.max(epi_dist)
    maxdiagonal = maxdist*(2**0.5)*1.10
    mag = np.array([sourcepar[k]['Mw'] for k in st_ids])
    magmean = sourcepar['means_weight']['Mw']
    magerr = sourcepar['errors_weight']['Mw']
    # fc = np.array([sourcepar[k]['fc'] for k in st_ids])
    g = Geod(ellps='WGS84')
    hypo = config.hypo
    lonmax, latmax, _ = g.fwd(
        hypo.longitude, hypo.latitude, 45, maxdiagonal*1000.)
    lonmin = 2*hypo.longitude - lonmax
    latmin = 2*hypo.latitude - latmax

    tile_dir = 'maptiles'
    stamen_terrain = CachedTiler(cimgt.Stamen('terrain-background'), tile_dir)
    # Create a GeoAxes
    figsize = (10, 10)
    if config.PLOT_SHOW:
        fig = plt.figure(figsize=figsize)
    else:
        fig = Figure(figsize=figsize)
    ax = fig.add_subplot(111, projection=stamen_terrain.crs)
    # Add event information as a title
    textstr = 'evid: {} \nlon: {:.3f} lat: {:.3f} ' +\
              'depth: {:.1f} km time: {}'
    textstr = textstr.format(
        hypo.evid, hypo.longitude, hypo.latitude, hypo.depth,
        hypo.origin_time.format_iris_web_service())
    ax.text(0., 1.15, textstr, fontsize=10,
            ha='left', va='top', linespacing=1.5, transform=ax.transAxes)
    if config.options.evname is not None:
        textstr = '{} — '.format(config.options.evname)
    else:
        textstr = ''
    textstr += 'Mw {:.1f} ± {:.1f}'.format(magmean, magerr)
    ax.text(0., 1.22, textstr, fontsize=14,
            ha='left', va='top', transform=ax.transAxes)
    trans = ccrs.Geodetic()
    ax.set_extent([lonmin, lonmax, latmin, latmax])
    if maxdiagonal <= 100:
        tile_zoom_level = 12
    else:
        tile_zoom_level = 8
    ax.add_image(stamen_terrain, tile_zoom_level)
    ax.gridlines(draw_labels=True, color='#777777', linestyle='--')
    countries = cfeature.NaturalEarthFeature(
        category='cultural',
        name='admin_0_countries',
        scale='10m',
        facecolor='none')
    ax.add_feature(countries, edgecolor='k')
    _plot_circles(ax, hypo.longitude, hypo.latitude, maxdist, 5)

    vmax = np.max(np.abs(mag-magmean))
    vmin = -vmax
    vmax += magmean
    vmin += magmean
    norm = colors.Normalize(vmin=vmin, vmax=vmax)
    cmap = cm.Spectral_r

    if st is not None:
        ax.scatter(
            lonlat_all[:, 0], lonlat_all[:, 1],
            marker='^', s=100,
            color='k', edgecolor='k',
            zorder=88, transform=trans)
    ax.scatter(
        lonlat[:, 0], lonlat[:, 1],
        marker='^', s=100,
        color=cmap(norm(mag)), edgecolor='k',
        zorder=99, transform=trans)

    # Add a colorbar
    ax_divider = make_axes_locatable(ax)
    cax = ax_divider.append_axes('right', size='100%',
                                 pad='-30%', aspect=15.,
                                 map_projection=stamen_terrain.crs)
    cax.background_patch.set_visible(False)
    cax.outline_patch.set_visible(False)
    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    fig.colorbar(sm, cax=cax)
    cax.get_yaxis().set_visible(True)
    cax.axhline(magmean, color='k')
    cax.axhline(magmean-magerr, linestyle=':', color='k')
    cax.axhline(magmean+magerr, linestyle=':', color='k')
    cm_label = 'Magnitude'
    cax.set_ylabel(cm_label)
    # Add code information at the figure bottom
    textstr = 'SourceSpec v{} '.format(get_git_version())
    textstr += '– {} {}\n'.format(
        config.end_of_run.strftime('%Y-%m-%d %H:%M:%S'),
        config.end_of_run_tz)
    cax.text(1., -0.1, textstr, fontsize=8,
             ha='right', va='top', transform=cax.transAxes)

    evid = config.hypo.evid
    figfile_base = os.path.join(config.options.outdir, evid + '.map_mag.')
    fmt = config.PLOT_SAVE_FORMAT
    if fmt == 'pdf_multipage':
        fmt = 'pdf'
    figfile = figfile_base + fmt
    if config.PLOT_SHOW:
        plt.show()
    if config.PLOT_SAVE:
        if config.PLOT_SHOW:
            fig.savefig(figfile, bbox_inches='tight')
        else:
            canvas = FigureCanvasAgg(fig)
            canvas.print_figure(figfile, bbox_inches='tight')
        logging.info('Station-magnitude map saved to: ' + figfile)
