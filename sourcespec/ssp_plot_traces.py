# -*- coding: utf-8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Trace plotting routine.

:copyright:
    2015-2023 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
import os
import contextlib
import numpy as np
import logging
from sourcespec.savefig import savefig
from sourcespec._version import get_versions
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.transforms as transforms
import matplotlib.patches as patches
import matplotlib.patheffects as PathEffects
from matplotlib.ticker import ScalarFormatter as sf
logger = logging.getLogger(__name__.split('.')[-1])
# Reduce logging level for Matplotlib to avoid DEBUG messages
mpl_logger = logging.getLogger('matplotlib')
mpl_logger.setLevel(logging.WARNING)


class ScalarFormatter(sf):  #NOQA
    def _set_format(self, vmin=None, vmax=None):
        self.format = '%1.1f'


phase_label_pos = {'P': 0.9, 'S': 0.93}
phase_label_color = {'P': 'black', 'S': 'black'}


def _nplots(config, st, maxlines, ncols):
    """Determine the number of lines and columns of the plot."""
    # Remove the channel letter to determine the number of plots
    if config.plot_traces_ignored:
        nplots = len({tr.id[:-1] for tr in st})
    else:
        nplots = len({tr.id[:-1] for tr in st if not tr.stats.ignore})
    nplots = len({tr.id[:-1] for tr in st})
    nlines = int(np.ceil(nplots/ncols))
    nlines = min(nlines, maxlines)
    if nplots < ncols:
        ncols = 1
    return nlines, ncols


def _make_fig(config, nlines, ncols):
    figsize = (16, 9) if nlines <= 3 else (16, 18)
    # high dpi needed to rasterize png
    # vector formats (pdf, svg) do not have rasters
    dpi = 300 if config.plot_save_format == 'png' else 72
    fig = plt.figure(figsize=figsize, dpi=dpi)
    # Create an invisible axis and use it for title and footer
    ax0 = fig.add_subplot(111, label='ax0')
    ax0.set_axis_off()
    # Add event information as a title
    hypo = config.hypo
    textstr = (
        f'evid: {hypo.evid} lon: {hypo.longitude:.3f} '
        f'lat: {hypo.latitude:.3f} depth: {hypo.depth:.1f} km'
    )
    with contextlib.suppress(AttributeError):
        textstr += f' time: {hypo.origin_time.format_iris_web_service()}'
    ax0.text(0., 1.06, textstr, fontsize=12,
             ha='left', va='top', transform=ax0.transAxes)
    if config.options.evname is not None:
        textstr = config.options.evname
        ax0.text(0., 1.1, textstr, fontsize=14,
                 ha='left', va='top', transform=ax0.transAxes)
    # Add code and author information at the figure bottom
    textstr = f'SourceSpec v{get_versions()["version"]} '
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
    ax0.text(1., -0.1, textstr, fontsize=10, linespacing=1.5,
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


def _savefig(config, figures):
    evid = config.hypo.evid
    figfile_base = os.path.join(config.options.outdir, f'{evid}.traces.')
    fmt = config.plot_save_format
    pad_inches = matplotlib.rcParams['savefig.pad_inches']
    bbox = figures[0].get_tightbbox(figures[0].canvas.get_renderer())
    bbox = bbox.padded(pad_inches)
    nfigures = len(figures)
    if nfigures == 1 or fmt == 'pdf_multipage':
        if fmt == 'pdf_multipage':
            figfile = f'{figfile_base}pdf'
            pdf = PdfPages(figfile)
        else:
            figfile = figfile_base + fmt
        figfiles = [figfile, ]
    else:
        figfiles = [f'{figfile_base}{n:02d}.{fmt}' for n in range(nfigures)]
    for n in range(nfigures):
        if fmt == 'pdf_multipage':
            pdf.savefig(figures[n], bbox_inches=bbox)
        else:
            savefig(figures[n], figfiles[n], fmt, bbox_inches=bbox)
        if not config.plot_show:
            plt.close(figures[n])
    for figfile in figfiles:
        logger.info(f'Trace plots saved to: {figfile}')
    config.figures['traces'] += figfiles
    if fmt == 'pdf_multipage':
        pdf.close()


def _plot_min_max(ax, x_vals, y_vals, linewidth, color, alpha, zorder):
    """Quick and dirty plot using less points. Useful for vector plotting."""
    ax_width_in_pixels = int(np.ceil(ax.bbox.width))
    nsamples = len(x_vals)
    samples_per_pixel = int(np.ceil(nsamples/ax_width_in_pixels))
    # find the closest multiple of samples_per_pixel (rounded down)
    nsamples -= nsamples % samples_per_pixel
    # resample x_vals
    x_vals = x_vals[:nsamples:samples_per_pixel]
    # reshape y_vals so that each row has a number of elements equal to
    # samples_per_pixel
    y_vals = y_vals[:nsamples].reshape(-1, samples_per_pixel)
    # find min and max for each row
    y_min = y_vals.min(axis=1)
    y_max = y_vals.max(axis=1)
    # double the number of elements in x_vals
    dx = x_vals[1] - x_vals[0]
    x_vals = np.column_stack((x_vals, x_vals+dx/2)).flatten()
    # alternate mins and maxs in y_vals
    y_vals = np.column_stack((y_min, y_max)).flatten()
    ax.plot(
        x_vals, y_vals, linewidth=linewidth, color=color, alpha=alpha,
        zorder=zorder)


def _freq_string(freq):
    """Return a string representing the rounded frequency."""
    if 1e-2 <= freq <= 1e2:
        # int or float notation for frequencies between 0.01 and 100
        int_freq = int(round(freq))
        return (
            f'{int_freq}' if np.abs(int_freq - freq) < 1e-1
            else f'{freq:.1f}'
        )
    else:
        # scientific notation for frequencies outside of 0.01 and 100
        freq_str = f'{freq:.1e}'
        m, n = map(float, freq_str.split('e'))
        n = int(n)
        int_m = int(round(m))
        return (
            f'{int_m}e{n}' if np.abs(int_m - m) < 1e-1
            else f'{m:.1f}e{n}'
        )


def _plot_trace(config, trace, ntraces, tmax,
                ax, ax_text, trans, trans3, path_effects):
    # Origin and height to draw vertical patches for noise and signal windows
    rectangle_patch_origin = 0
    rectangle_patch_height = 1
    orientation = trace.stats.channel[-1]
    if orientation in config.vertical_channel_codes:
        color = 'purple'
        if ntraces > 1:
            rectangle_patch_origin = 1./3
            rectangle_patch_height = 1./3
    if orientation in config.horizontal_channel_codes_1:
        color = 'green'
        if ntraces > 1:
            trace.data = (trace.data / tmax - 1) * tmax
            rectangle_patch_origin = 0
            rectangle_patch_height = 1./3
    if orientation in config.horizontal_channel_codes_2:
        color = 'blue'
        if ntraces > 1:
            trace.data = (trace.data / tmax + 1) * tmax
            rectangle_patch_origin = 2./3
            rectangle_patch_height = 1./3
    # dim out ignored traces
    alpha = 0.3 if trace.stats.ignore else 1.0
    if config.plot_save_format == 'png':
        ax.plot(trace.times(), trace, linewidth=1, color=color,
                alpha=alpha, zorder=20, rasterized=True)
    else:
        # reduce the number of points to plot for vector formats
        _plot_min_max(
            ax, trace.times(), trace.data, linewidth=1, color=color,
            alpha=alpha, zorder=20)
    ax.text(0.05, trace.data.mean(), trace.stats.channel,
            fontsize=8, color=color, transform=trans3, zorder=22,
            path_effects=path_effects)
    _text = f'S/N: {trace.stats.sn_ratio:.1f}'
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
    with contextlib.suppress(KeyError):
        N1 = trace.stats.arrivals['N1'][1] - trace.stats.starttime
        N2 = trace.stats.arrivals['N2'][1] - trace.stats.starttime
        rect = patches.Rectangle(
            (N1, rectangle_patch_origin),
            width=N2-N1, height=rectangle_patch_height,
            transform=trans, color='#eeeeee',
            zorder=-1)
        ax.add_patch(rect)
    # Signal window
    if config.wave_type[0] == 'S':
        t1 = trace.stats.arrivals['S1'][1] - trace.stats.starttime
        t2 = trace.stats.arrivals['S2'][1] - trace.stats.starttime
    elif config.wave_type[0] == 'P':
        t1 = trace.stats.arrivals['P1'][1] - trace.stats.starttime
        t2 = trace.stats.arrivals['P2'][1] - trace.stats.starttime
    rect = patches.Rectangle(
        (t1, rectangle_patch_origin),
        width=t2-t1, height=rectangle_patch_height,
        transform=trans, color='yellow',
        alpha=0.5, zorder=-1)
    ax.add_patch(rect)
    if not ax_text:
        text_y = 0.01
        color = 'black'
        id_no_channel = '.'.join(trace.id.split('.')[:-1])
        fmin_str = _freq_string(trace.stats.filter.freqmin)
        fmax_str = _freq_string(trace.stats.filter.freqmax)
        ax_text = (
            f'{id_no_channel} {trace.stats.instrtype} '
            f'{trace.stats.hypo_dist:.1f} km ({trace.stats.epi_dist:.1f} km)\n'
            f'filter: {fmin_str} - {fmax_str} Hz'
        )
        ax.text(0.05, text_y, ax_text,
                fontsize=8,
                horizontalalignment='left',
                verticalalignment='bottom',
                color=color,
                transform=ax.transAxes,
                zorder=50,
                path_effects=path_effects)
        ax_text = True
    # Reason why trace is ignored
    if trace.stats.ignore:
        _text = trace.stats.ignore_reason
        ax.text(0.5, trace.data.mean(), _text, ha='center',
                fontsize=8, color=color, transform=trans3, zorder=22,
                path_effects=path_effects)


def _add_labels(axes, plotn, ncols):
    """Add xlabels to the last row of plots."""
    # A row has "ncols" plots: the last row is from `plotn-ncols` to `plotn`
    n0 = max(plotn-ncols, 0)
    for ax in axes[n0:plotn]:
        ax.xaxis.set_tick_params(which='both', labelbottom=True)
        ax.set_xlabel('Time (s)', fontsize=8)


def _set_ylim(axes):
    """Set symmetric ylim."""
    for ax in axes:
        ylim = ax.get_ylim()
        ymax = np.max(np.abs(ylim))
        ax.set_ylim(-ymax, ymax)


def _trim_traces(config, st):
    for trace in st:
        t1 = (trace.stats.arrivals['P'][1] - config.noise_pre_time)
        t2 = (trace.stats.arrivals['S'][1] + 3 * config.win_length)
        trace.trim(starttime=t1, endtime=t2)


def plot_traces(config, st, ncols=None, block=True):
    """
    Plot traces in the original instrument unit (velocity or acceleration).

    Display to screen and/or save to file.
    """
    # Check config, if we need to plot at all
    if not config.plot_show and not config.plot_save:
        return

    matplotlib.rcParams['pdf.fonttype'] = 42  # to edit text in Illustrator

    if ncols is None:
        ntr = len({t.id[:-1] for t in st})
        ncols = 4 if ntr > 6 else 3

    nlines, ncols = _nplots(config, st, config.plot_traces_maxrows, ncols)
    fig, axes = _make_fig(config, nlines, ncols)
    figures = [fig]
    # Path effect to contour text in white
    path_effects = [PathEffects.withStroke(linewidth=3, foreground='white')]

    # Plot!
    plotn = 0
    if config.plot_traces_ignored:
        stalist = sorted({(tr.stats.hypo_dist, tr.id[:-1]) for tr in st})
    else:
        stalist = sorted(
            {
                (tr.stats.hypo_dist, tr.id[:-1])
                for tr in st
                if not tr.stats.ignore
            }
        )
    ax_text = False
    for t in stalist:
        plotn += 1
        _, traceid = t
        # 'code' is band+instrument code
        network, station, location, code = traceid.split('.')
        st_sel = st.select(network=network, station=station, location=location)
        if plotn > nlines*ncols:
            _set_ylim(axes)
            _add_labels(axes, plotn-1, ncols)
            fig, axes = _make_fig(config, nlines, ncols)
            figures.append(fig)
            plotn = 1
        ax = axes[plotn-1]
        if config.trace_units == 'auto':
            instrtype = [t.stats.instrtype for t in st_sel.traces
                         if t.stats.channel[:-1] == code][0]
        else:
            instrtype = config.trace_units
        if instrtype == 'acc':
            ax.set_ylabel(u'Acceleration (m/s²)', fontsize=8, labelpad=0)
        elif instrtype in ['broadb', 'shortp', 'vel']:
            ax.set_ylabel('Velocity (m/s)', fontsize=8, labelpad=0)
        elif instrtype in ['disp']:
            ax.set_ylabel('Displacement (m)', fontsize=8, labelpad=0)
        # Custom transformation for plotting phase labels:
        # x coords are data, y coords are axes
        trans = transforms.blended_transform_factory(ax.transData,
                                                     ax.transAxes)
        trans2 = transforms.blended_transform_factory(ax.transAxes,
                                                      ax.transData)
        trans3 = transforms.offset_copy(trans2, fig=fig, x=0, y=0.1)

        _trim_traces(config, st_sel)
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

    _set_ylim(axes)
    # Add labels for the last figure
    _add_labels(axes, plotn, ncols)
    # Turn off the unused axes
    for ax in axes[plotn:]:
        ax.set_axis_off()

    if config.plot_show:
        plt.show(block=block)
    if config.plot_save:
        _savefig(config, figures)
