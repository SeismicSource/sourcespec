# -*- coding: utf-8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Spectral plotting routine.

:copyright:
    2012 Claudio Satriano <satriano@ipgp.fr>

    2013-2014 Claudio Satriano <satriano@ipgp.fr>,
              Emanuela Matrullo <matrullo@geologie.ens.fr>

    2015-2023 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
import os
import math
import logging
import warnings
import contextlib
from collections import defaultdict
from obspy import Stream
from sourcespec.ssp_util import spec_minmax, moment_to_mag, mag_to_moment
from sourcespec.savefig import savefig
from sourcespec._version import get_versions
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.patheffects as PathEffects
logger = logging.getLogger(__name__.split('.')[-1])
# Reduce logging level for Matplotlib to avoid DEBUG messages
mpl_logger = logging.getLogger('matplotlib')
mpl_logger.setLevel(logging.WARNING)


synth_colors = [
    '#201F1F',
    '#94F75B',
    '#3EC2AA',
    '#FECC38',
    '#FC4384',
]


class PlotParams():

    def __init__(self):
        self.plot_type = None
        self.stack_plots = False
        self.nlines = 0
        self.ncols = 0
        self.freq_minmax = None
        self.moment_minmax = None
        self.mag_minmax = None
        self.plotn = 0
        self.figures = []
        self.axes = None
        self.ax0 = None

    def set_plot_params(self, config, spec_st, specnoise_st):
        """Determine the number of plots and axes min and max."""
        nplots = 0
        moment_minmax = None
        freq_minmax = None
        _spec_st = (
            spec_st
            if config.plot_spectra_ignored
            else Stream(sp for sp in spec_st if not sp.stats.ignore)
        )
        specids = {'.'.join(sp.id.split('.')[:-1]) for sp in _spec_st}
        for specid in specids:
            network, station, location = specid.split('.')
            spec_st_sel = _spec_st.select(
                network=network, station=station, location=location)
            if specnoise_st:
                specnoise_sel = specnoise_st.select(
                    network=network, station=station, location=location)
                if specnoise_sel[0].data.sum() != 0:
                    spec_st_sel += specnoise_sel
            for spec in spec_st_sel:
                moment_minmax, freq_minmax = spec_minmax(
                    spec.data, spec.get_freq(), moment_minmax, freq_minmax)
            # increment the number of plots,
            # based on the number of uniqs band+instrument codes
            nplots += len({x.stats.channel[:-1] for x in spec_st_sel})
        nlines = int(math.ceil(nplots / self.ncols))
        maxlines = config.plot_spectra_maxrows
        self.nlines = min(nlines, maxlines)
        if self.plot_type != 'weight':
            moment_minmax[1] *= 10
            self.mag_minmax = moment_to_mag(moment_minmax)
        self.freq_minmax = freq_minmax
        self.moment_minmax = moment_minmax


def _make_fig(config, plot_params):
    nlines = plot_params.nlines
    ncols = plot_params.ncols
    stack_plots = plot_params.stack_plots
    figsize = (16, 9) if nlines <= 3 or stack_plots else (16, 18)
    dpi = 100 if config.plot_show else 300
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
    textstr = (
        f'evid: {evid} lon: {ev_lon:.3f} lat: {ev_lat:.3f} '
        f'depth: {ev_depth:.1f} km'
    )
    with contextlib.suppress(AttributeError):
        textstr += f' time: {hypo.origin_time.format_iris_web_service()}'
    ax0.text(
        0., 1.06, textstr, fontsize=12,
        ha='left', va='top', transform=ax0.transAxes)
    if config.options.evname is not None:
        textstr = config.options.evname
        ax0.text(
            0., 1.1, textstr, fontsize=14,
            ha='left', va='top', transform=ax0.transAxes)
    # Add code and author information at the figure bottom
    textstr = f'SourceSpec v{get_versions()["version"]} '
    if not stack_plots:
        textstr += (
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
    ax0.text(
        1., -0.1, textstr, fontsize=10, linespacing=1.5,
        ha='right', va='top', transform=ax0.transAxes)
    axes = []
    for n in range(nlines * ncols):
        plotn = n + 1
        # ax1 has moment units (or weight)
        if plotn == 1:
            ax = (
                fig.add_subplot(1, 1, 1, label='ax')
                if stack_plots
                else fig.add_subplot(nlines, ncols, plotn)
            )
        elif not stack_plots:
            ax = fig.add_subplot(
                nlines, ncols, plotn,
                sharex=axes[0][0], sharey=axes[0][0])
        ax.set_xlim(plot_params.freq_minmax)
        with warnings.catch_warnings():
            # ignore warnings when ymin == ymax (ex., no weighting)
            warnings.simplefilter('ignore')
            ax.set_ylim(plot_params.moment_minmax)
        ax.grid(
            True, which='both', linestyle='solid', color='#DDDDDD', zorder=0)
        ax.set_axisbelow(True)
        ax.xaxis.set_tick_params(which='both', labelbottom=False)
        ax.yaxis.set_tick_params(which='both', labelleft=False)
        ax.tick_params(width=2)  # FIXME: ticks are below grid lines!
        if plot_params.plot_type == 'weight':
            ax2 = None
        elif ((stack_plots and plotn == 1) or not stack_plots):
            ax2 = ax.twinx()
            # ax2 has magnitude units
            ax2.set_ylim(plot_params.mag_minmax)
            ax2.yaxis.set_tick_params(
                which='both', labelright=False, width=0)
        ax.has_station_text = False
        axes.append((ax, ax2))
    fig.subplots_adjust(hspace=.025, wspace=.03)
    plot_params.figures.append(fig)
    plot_params.axes = axes
    plot_params.ax0 = ax0
    plot_params.plotn = 0


def _savefig(config, plottype, figures):
    evid = config.event.event_id
    if plottype == 'regular':
        suffix = '.ssp.'
        message = 'Spectral'
    elif plottype == 'weight':
        suffix = '.sspweight.'
        message = 'Weight'
    figfile_base = os.path.join(config.options.outdir, evid + suffix)
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
        logger.info(f'{message} plots saved to: {figfile}')
    config.figures[f'spectra_{plottype}'] += figfiles
    if fmt == 'pdf_multipage':
        pdf.close()


def _add_labels(plot_params):
    """
    Add xlabels to the last row plots.

    Add ylabels to the first and last columns.
    """
    plotn = plot_params.plotn
    ncols = plot_params.ncols
    plot_type = plot_params.plot_type
    axes = plot_params.axes
    # A row has "ncols" plots: the last row is from `plotn-ncols` to `plotn`
    n0 = max(plotn - ncols, 0)
    for ax, ax2 in axes[n0:plotn]:
        ax.xaxis.set_tick_params(which='both', labelbottom=True)
        ax.set_xlabel('Frequency (Hz)')
    # Show the y-labels only for the first column
    for i in range(0, len(axes) + ncols, ncols):
        try:
            ax = axes[i][0]
        except IndexError:
            continue
        try:
            # for ax2 we take the last column
            ax2 = axes[i - 1][1]
        except IndexError:
            continue
        ax.yaxis.set_tick_params(which='both', labelleft=True)
        ax.set_ylabel('Weight')
        if plot_type != 'weight':
            ax.set_ylabel('Seismic moment (N·m)')
            if ax2:
                ax2.yaxis.set_tick_params(
                    which='both', labelright=True, pad=0, width=2)
                ax2.set_ylabel('Magnitude')
    # still some work to do on the last plot
    ax, ax2 = axes[plotn - 1]
    if ax2:
        ax2.yaxis.set_tick_params(
            which='both', labelright=True, pad=0, width=2)
        ax2.set_ylabel('Magnitude')


def _color_lines(config, orientation, plotn, stack_plots):
    if orientation in config.vertical_channel_codes:
        color = 'purple'
        linestyle = 'solid'
        linewidth = 1
    if orientation in config.horizontal_channel_codes_1:
        color = 'green'
        linestyle = 'solid'
        linewidth = 1
    if orientation in config.horizontal_channel_codes_2:
        color = 'blue'
        linestyle = 'solid'
        linewidth = 1
    if orientation == 'H':
        # root sum of squares spectrum, corrected
        color = 'red'
        linestyle = 'solid'
        linewidth = 2
    elif orientation == 'S':
        # synthetic spectrum
        color = (
            synth_colors[(plotn - 1) % len(synth_colors)]
            if stack_plots
            else 'black'
        )
        linestyle = 'solid'
        linewidth = 2
    elif orientation == 'h':
        # root sum of squares spectrum, uncorrected
        color = 'chocolate'
        linestyle = 'solid'
        linewidth = 1
    elif orientation == 's':
        # synthetic spectrum, no attenuation
        color = 'gray'
        linestyle = 'solid'
        linewidth = 1
    elif orientation == 't':
        # synthetic spectrum, no corner frequency
        color = 'gray'
        linestyle = 'dashed'
        linewidth = 1
    return color, linestyle, linewidth


def _add_legend(config, plot_params, spec_st, specnoise_st):
    stack_plots = plot_params.stack_plots
    plot_type = plot_params.plot_type
    ax0 = plot_params.ax0
    # check the available channel codes
    channel_codes = {s.stats.channel[-1] for s in spec_st}
    ncol0 = 0
    handles0 = []
    if 'H' in channel_codes:
        ncol0 += 1
        orientation = 'H'
        color, linestyle, linewidth =\
            _color_lines(config, orientation, 0, stack_plots)
        linewidth = 2
        if plot_type == 'weight':
            label = 'Weight'
        else:
            label = 'Root sum of squares'
            if 'h' in channel_codes:
                label += ', corr.'
        _h, = ax0.plot(range(2), linestyle=linestyle, linewidth=linewidth,
                       color=color, label=label)
        handles0.append(_h)
    if 'h' in channel_codes:
        ncol0 += 1
        orientation = 'h'
        color, linestyle, linewidth =\
            _color_lines(config, orientation, 0, stack_plots)
        linewidth = 2
        label = 'RSS, uncorr.'
        _h, = ax0.plot(range(2), linestyle=linestyle, linewidth=linewidth,
                       color=color, label=label)
        handles0.append(_h)
    Z_codes = sorted(
        c for c in channel_codes if c in config.vertical_channel_codes)
    if len(Z_codes) > 0:
        ncol0 += 1
        orientation = Z_codes[0]
        color, linestyle, linewidth =\
            _color_lines(config, orientation, 0, stack_plots)
        linewidth = 2
        label = ', '.join(Z_codes)
        _h, = ax0.plot(range(2), linestyle=linestyle, linewidth=linewidth,
                       color=color, label=label)
        handles0.append(_h)
    H1_codes = sorted(
        c for c in channel_codes if c in config.horizontal_channel_codes_1)
    if len(H1_codes) > 0:
        ncol0 += 1
        orientation = H1_codes[0]
        color, linestyle, linewidth =\
            _color_lines(config, orientation, 0, stack_plots)
        linewidth = 2
        label = ', '.join(H1_codes)
        _h, = ax0.plot(range(2), linestyle=linestyle, linewidth=linewidth,
                       color=color, label=label)
        handles0.append(_h)
    H2_codes = sorted(
        c for c in channel_codes if c in config.horizontal_channel_codes_2)
    if len(H2_codes) > 0:
        ncol0 += 1
        orientation = H2_codes[0]
        color, linestyle, linewidth =\
            _color_lines(config, orientation, 0, stack_plots)
        linewidth = 2
        label = ', '.join(H2_codes)
        _h, = ax0.plot(range(2), linestyle=linestyle, linewidth=linewidth,
                       color=color, label=label)
        handles0.append(_h)
    if specnoise_st:
        ncol0 += 1
        linewidth = 2
        color = 'gray'
        linestyle = ':'
        label = 'Noise'
        _h, = ax0.plot(range(2), linestyle=linestyle, linewidth=linewidth,
                       color=color, label=label)
        handles0.append(_h)
    # Create a second axis for a second legend
    ax1 = ax0.get_figure().add_subplot(111, label='ax1', zorder=-1)
    ax1.set_axis_off()
    ncol1 = 0
    handles1 = []
    if 'S' in channel_codes:
        ncol1 += 1
        orientation = 'S'
        color, linestyle, linewidth =\
            _color_lines(config, orientation, 0, stack_plots)
        linewidth = 2
        label = 'Brune fit'
        _h, = ax1.plot(range(2), linestyle=linestyle, linewidth=linewidth,
                       color=color, label=label)
        handles1.append(_h)
    if 's' in channel_codes:
        ncol1 += 1
        orientation = 's'
        color, linestyle, linewidth =\
            _color_lines(config, orientation, 0, stack_plots)
        linewidth = 2
        label = 'Brune fit no attenuation'
        _h, = ax1.plot(range(2), linestyle=linestyle, linewidth=linewidth,
                       color=color, label=label)
        handles1.append(_h)
    if 't' in channel_codes:
        ncol1 += 1
        orientation = 't'
        color, linestyle, linewidth =\
            _color_lines(config, orientation, 0, stack_plots)
        linewidth = 2
        label = 'Brune fit no fc'
        _h, = ax1.plot(range(2), linestyle=linestyle, linewidth=linewidth,
                       color=color, label=label)
        handles1.append(_h)
    # Put the two legends on the two axes
    legend0_y = -0.111 if handles0 and handles1 else -0.127
    legend1_y = -0.147 if handles0 and handles1 else -0.127
    if handles0:
        ax0.legend(handles=handles0, bbox_to_anchor=(0, legend0_y),
                   loc='lower left', borderaxespad=0, ncol=ncol0)
    if handles1:
        ax1.legend(handles=handles1, bbox_to_anchor=(0, legend1_y),
                   loc='lower left', borderaxespad=0, ncol=ncol1)
    # Make lines invisible
    for h in handles0 + handles1:
        h.set_visible(False)


def _snratio_text(spec, ax, color, path_effects):
    global snratio_text_ypos
    if spec.stats.spectral_snratio is None:
        return
    snratio_text = f'S/N: {spec.stats.spectral_snratio:.1f}'
    ax.text(
        0.95, snratio_text_ypos, snratio_text, ha='right', va='top',
        fontsize=8, color=color, path_effects=path_effects,
        transform=ax.transAxes, zorder=20)
    snratio_text_ypos -= 0.05


def _station_text(spec, ax, color, path_effects, stack_plots):
    global station_text_ypos
    station_text = f'{spec.id[:-1]} {spec.stats.instrtype}'
    if not stack_plots:
        color = 'black'
    station_text += (
        f'\n{spec.stats.hypo_dist:.1f} km ({spec.stats.epi_dist:.1f} km)')
    if not ax.has_station_text or stack_plots:
        ax.text(
            0.05, station_text_ypos, station_text,
            horizontalalignment='left',
            verticalalignment='bottom',
            color=color,
            transform=ax.transAxes,
            zorder=50,
            path_effects=path_effects)
        ax.has_station_text = True


def _params_text(spec, ax, color, path_effects, stack_plots):
    global station_text_ypos
    if stack_plots:
        params_text_ypos = station_text_ypos - 0.04
    else:
        params_text_ypos = 0.03
        color = 'black'
    fc = spec.stats.par['fc']
    Mw = spec.stats.par['Mw']
    Mo = mag_to_moment(Mw)
    t_star = spec.stats.par['t_star']
    if 'par_err' in spec.stats.keys():
        fc_err_left, fc_err_right = spec.stats.par_err['fc']
        Mw_err_left, Mw_err_right = spec.stats.par_err['Mw']
        t_star_err_left, t_star_err_right =\
            spec.stats.par_err['t_star']
    else:
        fc_err_left = fc_err_right = 0.
        Mw_err_left = Mw_err_right = 0.
        t_star_err_left = t_star_err_right = 0.
    Mo_text = f'Mo: {Mo:.2g}N·m'
    Mw_err_text = (
        f'±{Mw_err_left:.2f}'
        if round(Mw_err_left, 2) == round(Mw_err_right, 2)
        else f'[-{Mw_err_left:.2f},+{Mw_err_right:.2f}]'
    )
    Mw_text = f'Mw: {Mw:.2f}{Mw_err_text}'
    fc_err_text = (
        f'±{fc_err_left:.2f}Hz'
        if round(fc_err_left, 2) == round(fc_err_right, 2)
        else f'[-{fc_err_left:.2f},+{fc_err_right:.2f}]Hz'
    )
    fc_text = f'fc: {fc:.2f}{fc_err_text}'
    t_star_err_text = (
        f'±{t_star_err_left:.2f}s'
        if round(t_star_err_left, 2) == round(t_star_err_right, 2)
        else f'[-{t_star_err_left:.2f},+{t_star_err_right:.2f}]s'
    )
    t_star_text = f't*: {t_star:.2f}{t_star_err_text}'
    if len(fc_text + t_star_text) > 38:
        sep = '\n'
        params_text_ypos -= 0.01
        station_text_ypos += 0.06
    else:
        sep = ' '
    params_text = f'{Mo_text} {Mw_text}\n{fc_text}{sep}{t_star_text}'
    ax.text(
        0.05, params_text_ypos, params_text,
        horizontalalignment='left',
        verticalalignment='bottom',
        color=color,
        fontsize=9,
        transform=ax.transAxes,
        zorder=50,
        path_effects=path_effects)


def _plot_fc_and_mw(spec, ax, ax2):
    fc = spec.stats.par['fc']
    Mw = spec.stats.par['Mw']
    if 'par_err' in spec.stats.keys():
        fc_err_left, fc_err_right = spec.stats.par_err['fc']
        fc_min = fc - fc_err_left
        if fc_min < 0:
            fc_min = 0.01
        ax.axvspan(
            fc_min, fc + fc_err_right, color='#bbbbbb', alpha=0.3, zorder=1)
        Mw_err_left, Mw_err_right = spec.stats.par_err['Mw']
        ax2.axhspan(
            Mw - Mw_err_left, Mw + Mw_err_right, color='#bbbbbb',
            alpha=0.3, zorder=1)
    ax.axvline(fc, color='#999999', linewidth=2., zorder=1)


def _plot_spec(config, plot_params, spec, spec_noise):
    """Plot one spectrum (and its associated noise)."""
    plotn = plot_params.plotn
    plot_type = plot_params.plot_type
    stack_plots = plot_params.stack_plots
    axes = plot_params.axes
    ax, ax2 = axes[plotn - 1]
    orientation = spec.stats.channel[-1]
    # Path effect to contour text in white
    path_effects = [PathEffects.withStroke(linewidth=3, foreground='white')]
    color, linestyle, linewidth =\
        _color_lines(config, orientation, plotn, stack_plots)
    # dim out ignored spectra
    alpha = 0.3 if spec.stats.ignore else 1.0
    if plot_type == 'regular':
        zorder = defaultdict(lambda: 20, {'S': 21, 'H': 22})
        ax.loglog(
            spec.freq_log, spec.data_log, color=color, alpha=alpha,
            linestyle=linestyle, linewidth=linewidth,
            zorder=zorder[orientation])
        special_orientations = ['S', 's', 't', 'H', 'h']
        # Write spectral S/N for regular Z,N,E components
        if orientation not in special_orientations:
            _snratio_text(spec, ax, color, path_effects)
        # Plot fc and Mw if a synthetic spectrum S is available
        if orientation == 'S':
            _plot_fc_and_mw(spec, ax, ax2)
    elif plot_type == 'weight':
        ax.semilogx(
            spec.get_freq(), spec.data, color=color, alpha=alpha,
            zorder=20)
    else:
        raise ValueError(f'Unknown plot type: {plot_type}')

    if spec_noise:
        ax.loglog(
            spec_noise.get_freq(), spec_noise.data,
            linestyle=':', linewidth=linewidth,
            color=color, alpha=alpha, zorder=20)
    if orientation == 'S':
        _params_text(spec, ax, color, path_effects, stack_plots)
    # station_text must be written after params_text, since it might move up
    # if params_text is too tall
    _station_text(spec, ax, color, path_effects, stack_plots)


def _plot_specid(config, plot_params, specid, spec_st, specnoise_st):
    """Plot all spectra having the same specid."""
    plotn = plot_params.plotn + 1
    nlines = plot_params.nlines
    ncols = plot_params.ncols
    # 'code' is band+instrument code
    network, station, location, code = specid.split('.')
    spec_st_sel = spec_st.select(
        network=network, station=station, location=location)
    if plotn > nlines * ncols:
        # Add labels and legend before making a new figure
        _add_labels(plot_params)
        _add_legend(config, plot_params, spec_st, specnoise_st)
        _make_fig(config, plot_params)
        plotn = 1
    plot_params.plotn = plotn
    special_orientations = ['S', 's', 't', 'H', 'h']
    orientations = [sp.stats.channel[-1] for sp in spec_st_sel]
    # compute the number of instrument components (N, Z, E, 1, 2, ...)
    ncomponents = len(
        [o for o in orientations if o not in special_orientations])
    global snratio_text_ypos
    global station_text_ypos
    snratio_text_ypos = 0.95
    if plot_params.stack_plots:
        station_text_ypos = 0.05 + 0.10 * (plotn - 1)
    else:
        station_text_ypos = 0.15
    # sort specs by orientation letter, so that synthetic is plotted first
    # sort_order[...] gives 10 if the key is not present
    sort_order = defaultdict(lambda: 10, {'S': 0, 's': 1, 't': 2, 'H': 99})
    spec_st_sel_sort = sorted(
        spec_st_sel, key=lambda s: sort_order[s.stats.channel[-1]])
    for spec in spec_st_sel_sort:
        if spec.stats.channel[:-1] != code:
            continue
        if not config.plot_spectra_ignored and spec.stats.ignore:
            continue
        # If there is only one component, do not plot the 'H' spectrum
        # which would coincide with the component spectrum
        # (but if the spectrum is corrected, e.g., the 'h' component is
        # available, then plot it)
        orientation = spec.stats.channel[-1]
        if (not plot_params.stack_plots and ncomponents == 1
                and orientation == 'H' and 'h' not in orientations):
            continue
        spec_noise = None
        if orientation not in ['S', 's', 't', 'h']:
            specid = spec.get_id()
            with contextlib.suppress(Exception):
                spec_noise = specnoise_st.select(id=specid)[0]
        _plot_spec(config, plot_params, spec, spec_noise)


def plot_spectra(config, spec_st, specnoise_st=None, ncols=None,
                 stack_plots=False, plot_type='regular'):
    """
    Plot spectra for signal and noise.

    Display to screen and/or save to file.
    """
    # Check config, if we need to plot at all
    if not config.plot_show and not config.plot_save:
        return

    matplotlib.rcParams['pdf.fonttype'] = 42  # to edit text in Illustrator

    if ncols is None:
        nspec = len({s.id[:-1] for s in spec_st})
        ncols = 4 if nspec > 6 else 3

    plot_params = PlotParams()
    plot_params.plot_type = plot_type
    plot_params.stack_plots = stack_plots
    plot_params.ncols = ncols
    plot_params.set_plot_params(config, spec_st, specnoise_st)
    _make_fig(config, plot_params)

    # Plot!
    if config.plot_spectra_ignored:
        stalist = sorted({(sp.stats.hypo_dist, sp.id[:-1]) for sp in spec_st})
    else:
        stalist = sorted({
            (sp.stats.hypo_dist, sp.id[:-1])
            for sp in spec_st
            if not sp.stats.ignore
        })
    for _, specid in stalist:
        _plot_specid(config, plot_params, specid, spec_st, specnoise_st)

    # Add labels and legend for the last figure
    _add_labels(plot_params)
    _add_legend(config, plot_params, spec_st, specnoise_st)
    # Turn off the unused axes
    for ax, ax2 in plot_params.axes[plot_params.plotn:]:
        ax.set_axis_off()
        if ax2:
            ax2.set_axis_off()

    if config.plot_show:
        plt.show()
    if config.plot_save:
        _savefig(config, plot_type, plot_params.figures)
