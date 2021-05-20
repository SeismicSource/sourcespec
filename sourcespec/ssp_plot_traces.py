# -*- coding: utf-8 -*-
"""
Trace plotting routine.

:copyright:
    2015-2021 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement, Version 2.1
    (http://www.cecill.info/index.en.html)
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import os
import math
import logging
from sourcespec._version import get_versions
logger = logging.getLogger(__name__.split('.')[-1])

import matplotlib
# Reduce logging level for Matplotlib to avoid DEBUG messages
mpl_logger = logging.getLogger('matplotlib')
mpl_logger.setLevel(logging.WARNING)
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.transforms as transforms
import matplotlib.patches as patches
import matplotlib.patheffects as PathEffects
from matplotlib.ticker import ScalarFormatter as sf
class ScalarFormatter(sf):  #NOQA
    def _set_format(self, vmin=None, vmax=None):
        self.format = '%1.1f'


phase_label_pos = {'P': 0.9, 'S': 0.93}
phase_label_color = {'P': 'black', 'S': 'black'}


def _nplots(config, st, maxlines, ncols):
    """Determine the number of lines and columns of the plot."""
    # Remove the channel letter to determine the number of plots
    if config.plot_traces_ignored:
        nplots = len(set(tr.id[:-1] for tr in st))
    else:
        nplots = len(set(tr.id[:-1] for tr in st if not tr.stats.ignore))
    nplots = len(set(tr.id[:-1] for tr in st))
    nlines = int(math.ceil(nplots/ncols))
    if nlines > maxlines:
        nlines = maxlines
    if nplots < ncols:
        ncols = 1
    return nlines, ncols


def _make_fig(config, nlines, ncols):
    if nlines <= 3:
        figsize = (16, 9)
    else:
        figsize = (16, 18)
    fig = plt.figure(figsize=figsize)
    # Create an invisible axis and use it for title and footer
    ax0 = fig.add_subplot(111, label='ax0')
    ax0.set_axis_off()
    # Add event information as a title
    hypo = config.hypo
    textstr = 'evid: {} lon: {:.3f} lat: {:.3f} depth: {:.1f} km'
    textstr = textstr.format(
        hypo.evid, hypo.longitude, hypo.latitude, hypo.depth)
    try:
        textstr += ' time: {}'.format(
            hypo.origin_time.format_iris_web_service())
    except AttributeError:
        pass
    ax0.text(0., 1.06, textstr, fontsize=12,
             ha='left', va='top', transform=ax0.transAxes)
    if config.options.evname is not None:
        textstr = config.options.evname
        ax0.text(0., 1.1, textstr, fontsize=14,
                 ha='left', va='top', transform=ax0.transAxes)
    # Add code information at the figure bottom
    textstr = 'SourceSpec v{} '.format(get_versions()['version'])
    ax0.text(1., -0.1, textstr, fontsize=10,
             ha='right', va='top', transform=ax0.transAxes)
    axes = []
    for n in range(nlines*ncols):
        plotn = n+1
        if plotn == 1:
            ax = fig.add_subplot(nlines, ncols, plotn)
        else:
            ax = fig.add_subplot(nlines, ncols, plotn, sharex=axes[0])
        ax.grid(True, which='both', linestyle='solid', color='#DDDDDD',
                zorder=0)
        ax.set_axisbelow(True)
        ax.xaxis.set_tick_params(which='both', labelbottom=False)
        ax.yaxis.set_tick_params(which='both', labelleft=True)
        ax.tick_params(width=2)  # FIXME: ticks are below grid lines!
        ax.tick_params(labelsize=8)
        ax.yaxis.offsetText.set_fontsize(8)
        yfmt = ScalarFormatter()
        yfmt.set_powerlimits((-1, 1))
        ax.yaxis.set_major_formatter(yfmt)
        axes.append(ax)
    fig.subplots_adjust(hspace=.1, wspace=.20)
    return fig, axes


def _savefig(config, figures, async_plotter):
    evid = config.hypo.evid
    figfile_base = os.path.join(config.options.outdir, evid + '.traces.')
    fmt = config.PLOT_SAVE_FORMAT
    pad_inches = matplotlib.rcParams['savefig.pad_inches']
    bbox = figures[0].get_tightbbox(figures[0].canvas.get_renderer())
    bbox = bbox.padded(pad_inches)
    if fmt == 'pdf_multipage':
        figfile = figfile_base + 'pdf'
        with PdfPages(figfile) as pdf:
            for fig in figures:
                pdf.savefig(fig, bbox_inches=bbox)
                if not config.PLOT_SHOW:
                    plt.close(fig)
        logger.info('Trace plots saved to: ' + figfile)
        return
    for n, fig in enumerate(figures):
        if len(figures) == 1:
            figfile = figfile_base + fmt
        else:
            figfile = figfile_base + '{:02d}.{}'.format(n, fmt)
        if config.PLOT_SHOW or (async_plotter is None):
            fig.savefig(figfile, bbox_inches=bbox)
        else:
            async_plotter.save(fig, figfile, bbox_inches=bbox)
        logger.info('Trace plots saved to: ' + figfile)
        if not config.PLOT_SHOW:
            plt.close(fig)


def _plot_trace(config, trace, ntraces, tmax,
                ax, ax_text, trans, trans3, path_effects):
    t1 = (trace.stats.arrivals['P'][1] - config.pre_p_time)
    t2 = (trace.stats.arrivals['S'][1] + 3 * config.win_length)
    trace.trim(starttime=t1, endtime=t2, pad=True, fill_value=0)
    orientation = trace.stats.channel[-1]
    if orientation in config.vertical_channel_codes:
        color = 'purple'
    if orientation in config.horizontal_channel_codes_1:
        color = 'green'
        if ntraces > 1:
            trace.data = (trace.data / tmax - 1) * tmax
    if orientation in config.horizontal_channel_codes_2:
        color = 'blue'
        if ntraces > 1:
            trace.data = (trace.data / tmax + 1) * tmax
    # dim out ignored traces
    if trace.stats.ignore:
        alpha = 0.3
    else:
        alpha = 1.0
    ax.plot(trace.times(), trace, color=color,
            alpha=alpha, zorder=20, rasterized=True)
    ax.text(0.05, trace.data.mean(), trace.stats.channel,
            fontsize=8, color=color, transform=trans3, zorder=22,
            path_effects=path_effects)
    _text = 'S/N: {:.1f}'.format(trace.stats.sn_ratio)
    ax.text(0.95, trace.data.mean(), _text, ha='right',
            fontsize=8, color=color, transform=trans3, zorder=22,
            path_effects=path_effects)
    for phase in 'P', 'S':
        a = trace.stats.arrivals[phase][1] - trace.stats.starttime
        text = trace.stats.arrivals[phase][0]
        ax.axvline(a, linestyle='--',
                   color=phase_label_color[phase], zorder=21)
        ax.text(a, phase_label_pos[phase], text,
                fontsize=8, transform=trans,
                zorder=22, path_effects=path_effects)
    # Noise window
    try:
        N1 = trace.stats.arrivals['N1'][1] - trace.stats.starttime
        N2 = trace.stats.arrivals['N2'][1] - trace.stats.starttime
        rect = patches.Rectangle((N1, 0), width=N2-N1, height=1,
                                 transform=trans, color='#eeeeee',
                                 alpha=0.5, zorder=-1)
        ax.add_patch(rect)
    except KeyError:
        pass
    # S-wave window
    S1 = trace.stats.arrivals['S1'][1] - trace.stats.starttime
    S2 = trace.stats.arrivals['S2'][1] - trace.stats.starttime
    rect = patches.Rectangle((S1, 0), width=S2-S1, height=1,
                             transform=trans, color='yellow',
                             alpha=0.5, zorder=-1)
    ax.add_patch(rect)
    if not ax_text:
        text_y = 0.1
        color = 'black'
        id_no_channel = '.'.join(trace.id.split('.')[:-1])
        ax_text = '%s %s %.1f km (%.1f km)' %\
                  (id_no_channel,
                   trace.stats.instrtype,
                   trace.stats.hypo_dist,
                   trace.stats.epi_dist)
        ax.text(0.05, text_y, ax_text,
                fontsize=8,
                horizontalalignment='left',
                verticalalignment='bottom',
                color=color,
                transform=ax.transAxes,
                zorder=50,
                path_effects=path_effects)
        ax_text = True


def _add_labels(axes, plotn, ncols):
    """Add xlabels to the last row of plots."""
    # A row has "ncols" plots: the last row is from `plotn-ncols` to `plotn`
    n0 = plotn-ncols if plotn-ncols > 0 else 0
    for ax in axes[n0:plotn]:
        ax.xaxis.set_tick_params(which='both', labelbottom=True)
        ax.set_xlabel('Time (s)', fontsize=8)


def plot_traces(config, st, spec_st=None, ncols=4, block=True,
                async_plotter=None):
    """
    Plot displacement traces.

    Display to screen and/or save to file.
    """
    # Check config, if we need to plot at all
    if not config.PLOT_SHOW and not config.PLOT_SAVE:
        return

    matplotlib.rcParams['pdf.fonttype'] = 42  # to edit text in Illustrator

    nlines, ncols = _nplots(config, st, config.plot_traces_maxrows, ncols)
    figures = []
    fig, axes = _make_fig(config, nlines, ncols)
    figures.append(fig)

    # Path effect to contour text in white
    path_effects = [PathEffects.withStroke(linewidth=3, foreground="white")]

    # Plot!
    plotn = 0
    if config.plot_traces_ignored:
        stalist = sorted(set((tr.stats.hypo_dist, tr.id[:-1]) for tr in st))
    else:
        stalist = sorted(set(
            (tr.stats.hypo_dist, tr.id[:-1]) for tr in st
            if not tr.stats.ignore
        ))
    for t in stalist:
        plotn += 1
        _, traceid = t
        # 'code' is band+instrument code
        network, station, location, code = traceid.split('.')
        st_sel = st.select(network=network, station=station, location=location)
        if plotn > nlines*ncols:
            _add_labels(axes, plotn-1, ncols)
            fig, axes = _make_fig(config, nlines, ncols)
            figures.append(fig)
            plotn = 1
        ax_text = False
        ax = axes[plotn-1]
        instrtype = [t.stats.instrtype for t in st_sel.traces
                     if t.stats.channel[:-1] == code][0]
        if instrtype == 'acc':
            ax.set_ylabel(u'Acceleration (m/s²)', fontsize=8, labelpad=0)
        elif instrtype == 'shortp' or instrtype == 'broadb':
            ax.set_ylabel('Velocity (m/s)', fontsize=8, labelpad=0)
        # Custom transformation for plotting phase labels:
        # x coords are data, y coords are axes
        trans = transforms.blended_transform_factory(ax.transData,
                                                     ax.transAxes)
        trans2 = transforms.blended_transform_factory(ax.transAxes,
                                                      ax.transData)
        trans3 = transforms.offset_copy(trans2, fig=fig, x=0, y=0.1)

        maxes = [abs(t.max()) for t in st_sel.traces
                 if t.stats.channel[:-1] == code]
        ntraces = len(maxes)
        tmax = max(maxes)
        for trace in st_sel.traces:
            if trace.stats.channel[:-1] != code:
                continue
            if not config.plot_traces_ignored and trace.stats.ignore:
                continue
            _plot_trace(
                config, trace, ntraces, tmax, ax, ax_text,
                trans, trans3, path_effects)

    # Add lables for the last figure
    _add_labels(axes, plotn, ncols)
    # Turn off the unused axes
    for ax in axes[plotn:]:
        ax.set_axis_off()

    if config.PLOT_SHOW:
        plt.show(block=block)
    if config.PLOT_SAVE:
        _savefig(config, figures, async_plotter)
