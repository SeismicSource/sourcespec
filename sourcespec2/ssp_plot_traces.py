# -*- coding: utf-8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Trace plotting routine.

:copyright:
    2015-2026 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
import os
import contextlib
import logging
import numpy as np
from obspy.core import Stream
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import transforms
from matplotlib import patches
import matplotlib.patheffects as PathEffects
from matplotlib.ticker import ScalarFormatter as sf
from .savefig import savefig
from ._version import get_versions
logger = logging.getLogger(__name__.rsplit('.', maxsplit=1)[-1])
# Reduce logging level for Matplotlib to avoid DEBUG messages
logging.getLogger('matplotlib').setLevel(logging.WARNING)


class ScalarFormatter(sf):
    """A ScalarFormatter with a custom format."""
    def _set_format(self, vmin=None, vmax=None):
        # pylint: disable=unused-argument
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
    nlines = int(np.ceil(nplots / ncols))
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
    evid = config.event.event_id
    hypo = config.event.hypocenter
    ev_lon = hypo.longitude.value_in_deg
    ev_lat = hypo.latitude.value_in_deg
    ev_depth = hypo.depth.value_in_km
    textstr = f'{config.options.evname} — ' if config.options.evname else ''
    textstr += (
        f'evid: {evid} '
        f'lon: {ev_lon:.3f} lat: {ev_lat:.3f} depth: {ev_depth:.1f} km'
    )
    with contextlib.suppress(AttributeError):
        textstr += f' time: {hypo.origin_time.format_iris_web_service()}'
    ax0.text(0., 1.06, textstr, fontsize=12,
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
    for n in range(nlines * ncols):
        plotn = n + 1
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


# Keep track of saved figure numbers to avoid saving the same figure twice
SAVED_FIGURE_CODES = []
# Bounding box for saving figures
BBOX = None


def _savefig(config, figures, suffix, force_numbering=False):
    global BBOX  # pylint: disable=global-statement
    evid = config.event.event_id
    figfile_base = os.path.join(config.options.outdir, f'{evid}.traces.')
    if suffix is not None:
        figfile_base += f'{suffix}.'
    fmt = config.plot_save_format
    if BBOX is None:
        # draw the figure without rendering it,
        # so that we can get the bounding box
        figures[0].draw_without_rendering()
        BBOX = figures[0].get_tightbbox(figures[0].canvas.get_renderer())
        pad_inches = matplotlib.rcParams['savefig.pad_inches']
        BBOX = BBOX.padded(pad_inches)
    nfigures = len(figures)
    if (nfigures == 1 or fmt == 'pdf_multipage') and not force_numbering:
        if fmt == 'pdf_multipage':
            figfile = f'{figfile_base}pdf'
            pdf = PdfPages(figfile)
        else:
            figfile = figfile_base + fmt
        figfiles = [figfile, ]
    else:
        figfiles = [f'{figfile_base}{n:02d}.{fmt}' for n in range(nfigures)]
    for n in range(nfigures):
        figcode = f'{suffix}.{n}' if suffix is not None else str(n)
        if figcode in SAVED_FIGURE_CODES:
            continue
        if fmt == 'pdf_multipage':
            pdf.savefig(figures[n], bbox_inches=BBOX)
        else:
            savefig(figures[n], figfiles[n], fmt, bbox_inches=BBOX)
        if not config.plot_show:
            plt.close(figures[n])
            # dereference the figure to free up memory
            figures[n] = None
        SAVED_FIGURE_CODES.append(figcode)
        figkey = f'traces_{suffix}' if suffix is not None else 'traces'
        config.figures[figkey].append(figfiles[n])
        logger.info(f'Trace plots saved to: {figfiles[n]}')
    if fmt == 'pdf_multipage':
        pdf.close()


def _plot_min_max(ax, x_vals, y_vals, linewidth, color, alpha, zorder):
    """Quick and dirty plot using less points. Useful for vector plotting."""
    ax_width_in_pixels = int(np.ceil(ax.bbox.width))
    nsamples = len(x_vals)
    samples_per_pixel = int(np.ceil(nsamples / ax_width_in_pixels))
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
    x_vals = np.column_stack((x_vals, x_vals + dx / 2)).flatten()
    # alternate mins and maxs in y_vals
    y_vals = np.column_stack((y_min, y_max)).flatten()
    ax.plot(
        x_vals, y_vals, linewidth=linewidth, color=color, alpha=alpha,
        zorder=zorder)


def _freq_string(freq):
    """Return a string representing the rounded frequency."""
    # float .2f notation for frequencies between 0.01 and 0.1
    if 1e-2 <= freq < 1e-1:
        return f'{freq:.2f}'
    # int or float .1f notation for frequencies between 0.1 and 100
    if 1e-1 <= freq <= 1e2:
        int_freq = int(round(freq))
        return (
            f'{int_freq}' if np.abs(int_freq - freq) < 1e-1
            else f'{freq:.1f}'
        )
    # scientific notation for frequencies outside of 0.01 and 100
    freq_str = f'{freq:.1e}'
    m, n = map(float, freq_str.split('e'))
    n = int(n)
    int_m = int(round(m))
    return (
        f'{int_m}e{n}' if np.abs(int_m - m) < 1e-1
        else f'{m:.1f}e{n}'
    )


def _plot_trace(config, trace, ntraces, tmax, ax, trans, trans3, path_effects):
    # Origin and height to draw vertical patches for noise and signal windows
    rectangle_patch_origin = 0
    rectangle_patch_height = 1
    orientation = trace.stats.channel[-1]
    if orientation in config.vertical_channel_codes:
        color = 'purple'
        if ntraces > 1:
            rectangle_patch_origin = 1. / 3
            rectangle_patch_height = 1. / 3
    # shift trace data for the horizontal components to avoid overlap
    # avoid division by zero
    tmax = 1 if tmax == 0 else tmax
    if orientation in config.horizontal_channel_codes_1:
        color = 'green'
        if ntraces > 1:
            trace.data = (trace.data / tmax - 1) * tmax
            rectangle_patch_origin = 0
            rectangle_patch_height = 1. / 3
    if orientation in config.horizontal_channel_codes_2:
        color = 'blue'
        if ntraces > 1:
            trace.data = (trace.data / tmax + 1) * tmax
            rectangle_patch_origin = 2. / 3
            rectangle_patch_height = 1. / 3
    # dim out ignored traces
    alpha = 0.3 if getattr(trace.stats, 'ignore', False) else 1.0
    times = trace.times() + trace.stats.time_offset
    if config.plot_save_format == 'png':
        ax.plot(
            times, trace.data, linewidth=1, color=color,
            alpha=alpha, zorder=20, rasterized=True)
    else:
        # reduce the number of points to plot for vector formats
        _plot_min_max(
            ax, times, trace.data, linewidth=1, color=color,
            alpha=alpha, zorder=20)
    # Calculate y position for text labels
    # trace.data can be a masked array, so handle it appropriately
    if np.ma.isMaskedArray(trace.data):
        data_mean = trace.data.mean()
        # Masked array mean() returns a masked scalar if all values are masked
        data_mean = 0.0 if np.ma.is_masked(data_mean) else float(data_mean)
    else:
        data_mean = float(np.nanmean(trace.data))
    if not np.isfinite(data_mean):
        # Use 0 as fallback if all data is NaN or empty
        data_mean = 0.0
    ax.text(0.05, data_mean, trace.stats.channel,
            fontsize=8, color=color, transform=trans3, zorder=22,
            path_effects=path_effects)
    with contextlib.suppress(AttributeError):
        sn_ratio = trace.stats.sn_ratio
        _text = f'S/N: {sn_ratio:.1f}'
        ax.text(0.95, data_mean, _text, ha='right',
                fontsize=8, color=color, transform=trans3, zorder=22,
                path_effects=path_effects)
    starttime = trace.stats.starttime - trace.stats.time_offset
    for phase in 'P', 'S':
        try:
            a = trace.stats.arrivals[phase][1] - starttime
        except KeyError:
            continue
        text = trace.stats.arrivals[phase][0]
        ax.axvline(a, linestyle='--',
                   color=phase_label_color[phase], zorder=21)
        ax.text(a, phase_label_pos[phase], text,
                fontsize=8, transform=trans,
                zorder=22, path_effects=path_effects)
    # Noise window
    with contextlib.suppress(KeyError):
        N1 = trace.stats.arrivals['N1'][1] - starttime
        N2 = trace.stats.arrivals['N2'][1] - starttime
        rect = patches.Rectangle(
            (N1, rectangle_patch_origin),
            width=(N2 - N1), height=rectangle_patch_height,
            transform=trans, color='#eeeeee',
            zorder=-1)
        ax.add_patch(rect)
        if config.wave_type[0] == 'S':
            t1 = trace.stats.arrivals['S1'][1] - starttime
            t2 = trace.stats.arrivals['S2'][1] - starttime
        elif config.wave_type[0] == 'P':
            t1 = trace.stats.arrivals['P1'][1] - starttime
            t2 = trace.stats.arrivals['P2'][1] - starttime
        rect = patches.Rectangle(
            (t1, rectangle_patch_origin),
            width=(t2 - t1), height=rectangle_patch_height,
            transform=trans, color='yellow',
            alpha=0.5, zorder=-1)
        ax.add_patch(rect)
    # Reason why trace is ignored
    if getattr(trace.stats, 'ignore', False):
        _text = trace.stats.ignore_reason
        color = 'black'
        ax.text(
            0.5, data_mean, _text, ha='center',
            fontsize=8, color=color, transform=trans3, zorder=22,
            path_effects=path_effects)
    _add_station_info_text(trace, ax, path_effects)


def _add_station_info_text(trace, ax, path_effects):
    with contextlib.suppress(AttributeError):
        if ax.has_station_info_text:
            return
    text_y = 0.01
    color = 'black'
    id_no_channel = '.'.join(trace.id.split('.')[:-1])
    station_info_text = (
        f'{id_no_channel} {trace.stats.instrtype} '
        f'{trace.stats.hypo_dist:.1f} km ({trace.stats.epi_dist:.1f} km)'
    )
    if getattr(trace.stats, 'processed', False):
        fmin_str = _freq_string(trace.stats.filter.freqmin)
        fmax_str = _freq_string(trace.stats.filter.freqmax)
        station_info_text += f'\nfilter: {fmin_str} - {fmax_str} Hz'
    else:
        station_info_text += '\nraw trace'
    ax.text(
        0.05, text_y, station_info_text, fontsize=8,
        horizontalalignment='left', verticalalignment='bottom', color=color,
        transform=ax.transAxes, zorder=50, path_effects=path_effects)
    ax.has_station_info_text = True


def _add_labels(axes, plotn, ncols):
    """Add xlabels to the last row of plots."""
    # A row has "ncols" plots: the last row is from `plotn-ncols` to `plotn`
    n0 = max(plotn - ncols, 0)
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
    """
    Trim traces and compute time offsets for plotting.

    Trimming starts at noise window start (N1) and ends at signal window end
    (S2) plus 2x window length. Pads with masked values if window exceeds trace
    bounds.
    Sets time_offset attribute relative to earliest trace start for alignment.
    """
    for trace in st:
        try:
            t1 = trace.stats.arrivals['N1'][1]
            t2 = trace.stats.arrivals['S2'][1] + 2 * config.win_length
            if t2 < t1:
                # this can happen for raw traces
                raise KeyError
        except KeyError:
            t1 = trace.stats.starttime
            t2 = trace.stats.endtime
        # Trim trace to the specified time window. If the requested window
        # extends beyond the trace boundaries, pad with masked values to
        # maintain proper time alignment without plotting spurious data.
        trace.trim(starttime=t1, endtime=t2, pad=True, fill_value=None)
    # compute time offset for correctly aligning traces when plotting
    min_starttime = min(tr.stats.starttime for tr in st)
    for trace in st:
        trace.stats.time_offset = trace.stats.starttime - min_starttime


def _get_ylabel(config, st_sel, processed):
    if config.correct_instrumental_response and not processed:
        return 'Counts'
    if config.trace_units == 'auto':
        instrtype = [t.stats.instrtype for t in st_sel.traces]
        if len(set(instrtype)) > 1:
            raise ValueError(
                'All traces with the same band+instrument code must have the '
                'same instrument type')
        instrtype = instrtype[0]
    else:
        instrtype = config.trace_units
    if instrtype in ['acc']:
        return 'Acceleration (m/s²)'
    if instrtype in ['broadb', 'shortp', 'vel']:
        return 'Velocity (m/s)'
    if instrtype in ['disp']:
        return 'Displacement (m)'
    raise ValueError(f'Unknown instrument type: {instrtype}')


def plot_traces(config, st, ncols=None, block=True, suffix=None):
    """
    Plot raw (counts) or processed traces (instrument units and filtered).

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
        stalist = sorted({
            (getattr(tr.stats, 'hypo_dist', 0), tr.id[:-1])
            for tr in st
        })
    else:
        stalist = sorted({
            (getattr(tr.stats, 'hypo_dist', 0), tr.id[:-1])
            for tr in st
            if not tr.stats.ignore
        })
    for _, traceid in stalist:
        # select traces with same band+instrument code
        network, station, location, code = traceid.split('.')
        st_sel = st.select(
            network=network, station=station, location=location,
            channel=f'{code}*')
        if not config.plot_traces_ignored:
            st_sel = Stream(tr for tr in st_sel if not tr.stats.ignore)
        processed = [getattr(tr.stats, 'processed', False) for tr in st_sel]
        if len(set(processed)) > 1:
            raise ValueError(
                'All traces with the same band+instrument code must have the '
                'same processed status')
        processed = processed[0]
        if not st_sel:
            continue
        plotn += 1
        if plotn > nlines * ncols:
            _set_ylim(axes)
            _add_labels(axes, plotn - 1, ncols)
            if (
                config.plot_save_asap and
                config.plot_save and not config.plot_show and
                config.plot_save_format != 'pdf_multipage'
            ):
                # save figure here to free up memory
                _savefig(config, figures, suffix, force_numbering=True)
            fig, axes = _make_fig(config, nlines, ncols)
            figures.append(fig)
            plotn = 1
        ax = axes[plotn - 1]
        ylabel = _get_ylabel(config, st_sel, processed)
        ax.set_ylabel(ylabel, fontsize=8, labelpad=0)
        # Custom transformation for plotting phase labels:
        # x coords are data, y coords are axes
        trans =\
            transforms.blended_transform_factory(ax.transData, ax.transAxes)
        trans2 =\
            transforms.blended_transform_factory(ax.transAxes, ax.transData)
        trans3 = transforms.offset_copy(trans2, fig=fig, x=0, y=0.1)

        _trim_traces(config, st_sel)
        max_values = [abs(tr.max()) for tr in st_sel if len(tr.data)]
        ntraces = len(max_values)
        if ntraces == 0:
            continue
        tmax = max(max_values)
        for trace in st_sel:
            if len(trace.data) == 0:
                continue
            _plot_trace(
                config, trace, ntraces, tmax, ax, trans, trans3, path_effects)

    _set_ylim(axes)
    # Add labels for the last figure
    _add_labels(axes, plotn, ncols)
    # Turn off the unused axes
    for ax in axes[plotn:]:
        ax.set_axis_off()

    if config.plot_show:
        plt.show(block=block)
    if config.plot_save:
        _savefig(config, figures, suffix)
