# -*- coding: utf-8 -*-
"""
Spectral plotting routine.

:copyright:
    2012 Claudio Satriano <satriano@ipgp.fr>

    2013-2014 Claudio Satriano <satriano@ipgp.fr>,
              Emanuela Matrullo <matrullo@geologie.ens.fr>

    2015-2017 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement, Version 2.1
    (http://www.cecill.info/index.en.html)
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import os
import math
import logging
from sourcespec.ssp_util import spec_minmax, moment_to_mag, mag_to_moment

synth_colors = [
    '#201F1F',
    '#94F75B',
    '#3EC2AA',
    '#FECC38',
    '#FC4384',
]


def plot_spectra(config, spec_st, specnoise_st=None, ncols=4,
                 stack_plots=False, plottype='regular', async_plotter=None):
    """
    Plot spectra for signal and noise.

    Display to screen and/or save to file.
    """
    # Check config, if we need to plot at all
    if not config.PLOT_SHOW and not config.PLOT_SAVE:
        return
    import matplotlib
    matplotlib.rcParams['pdf.fonttype'] = 42  # to edit text in Illustrator
    if config.PLOT_SHOW:
        import matplotlib.pyplot as plt
    else:
        from matplotlib.figure import Figure
        from matplotlib.backends.backend_agg import FigureCanvasAgg
    import matplotlib.patheffects as PathEffects

    # Determine the number of plots and axes min and max:
    nplots = 0
    moment_minmax = None
    freq_minmax = None
    for station in set(x.stats.station for x in spec_st.traces):
        spec_st_sel = spec_st.select(station=station)
        if specnoise_st:
            specnoise_sel = specnoise_st.select(station=station)
            spec_st_sel += specnoise_sel
        for spec in spec_st_sel.traces:
            moment_minmax, freq_minmax =\
                spec_minmax(spec.data, spec.get_freq(),
                            moment_minmax, freq_minmax)
        # 'code' is band+instrument code
        for code in set(x.stats.channel[0:2] for x in spec_st_sel):
            nplots += 1
    nlines = int(math.ceil(nplots/ncols))
    if plottype != 'weight':
        moment_minmax[1] *= 10
        mag_minmax = moment_to_mag(moment_minmax)

    # OK, now we can plot!
    if nlines <= 3 or stack_plots:
        figsize = (16, 9)
    else:
        figsize = (16, 18)
    if config.PLOT_SHOW:
        fig = plt.figure(figsize=figsize)
    else:
        fig = Figure(figsize=figsize)
    fig.subplots_adjust(hspace=.025, wspace=.03)

    # Path effect to contour text in white
    path_effects = [PathEffects.withStroke(linewidth=3, foreground="white")]

    # Plot!
    axes = []
    plotn = 0
    stalist = [s[1] for s in sorted(set((t.stats.hypo_dist, t.stats.station)
                                        for t in spec_st))]
    for station in stalist:
        spec_st_sel = spec_st.select(station=station)
        # 'code' is band+instrument code
        for code in sorted(set(x.stats.channel[0:2] for x in spec_st_sel)):
            plotn += 1
            ax_text = False

            # ax1 has moment units (or weight)
            if plotn == 1:
                if stack_plots:
                    ax = fig.add_subplot(1, 1, 1)
                else:
                    ax = fig.add_subplot(nlines, ncols, plotn)
            else:
                if not stack_plots:
                    ax = fig.add_subplot(nlines, ncols, plotn,
                                         sharex=axes[0][0], sharey=axes[0][0])
            ax.set_xlim(freq_minmax)
            ax.set_ylim(moment_minmax)
            ax.grid(True, which='both', linestyle='solid', color='#DDDDDD',
                    zorder=0)
            ax.set_axisbelow(True)
            [t.set_visible(False) for t in ax.get_xticklabels()]
            [t.set_visible(False) for t in ax.get_yticklabels()]
            ax.tick_params(width=2)  # FIXME: ticks are below grid lines!

            # ax2 has magnitude units
            if plottype != 'weight':
                if ((stack_plots and plotn == 1) or not stack_plots):
                    ax2 = ax.twinx()
                    ax2.set_ylim(mag_minmax)
                    [t.set_visible(False) for t in ax2.get_xticklabels()]
                    for tick in ax2.yaxis.get_major_ticks():
                        tick.set_pad(-2)
                        tick.label2.set_horizontalalignment('right')
                    ax2.yaxis.set_tick_params(width=0)
            else:
                ax2 = None
            axes.append((ax, ax2))

            for spec in spec_st_sel.traces:
                if spec.stats.channel[0:2] != code:
                    continue
                orientation = spec.stats.channel[2]
                if orientation in ['Z', '1']:
                    color = 'purple'
                    linestyle = 'solid'
                    linewidth = 1
                if orientation in ['N', '2', 'R']:
                    color = 'green'
                    linestyle = 'solid'
                    linewidth = 1
                if orientation in ['E', '3', 'T']:
                    color = 'blue'
                    linestyle = 'solid'
                    linewidth = 1
                if orientation == 'H':
                    color = 'red'
                    linestyle = 'solid'
                    linewidth = 1
                if orientation == 'S':
                    if stack_plots:
                        color = synth_colors[(plotn-1) % len(synth_colors)]
                    else:
                        color = 'black'
                    linestyle = 'solid'
                    linewidth = 2
                if orientation == 's':
                    color = 'gray'
                    linestyle = 'solid'
                    linewidth = 1
                if orientation == 't':
                    color = 'gray'
                    linestyle = 'dashed'
                    linewidth = 1
                if plottype in ['regular', 'noise']:
                    ax.loglog(spec.get_freq(), spec.data, color=color,
                              linestyle=linestyle, linewidth=linewidth,
                              zorder=20)
                    if orientation == 'S':
                        fc = spec.stats.par['fc']
                        fc_err = spec.stats.par_err['fc']
                        fc_min = fc-fc_err
                        if fc_min < 0:
                            fc_min = 0.01
                        ax.axvspan(fc_min, fc+fc_err, color='#bbbbbb',
                                   alpha=0.3, zorder=1)
                        ax.axvline(fc, color='#999999',
                                   linewidth=2., zorder=1)
                        Mw = spec.stats.par['Mw']
                        Mw_err = spec.stats.par_err['Mw']
                        ax2.axhspan(Mw-Mw_err, Mw+Mw_err, color='#bbbbbb',
                                    alpha=0.3, zorder=1)
                elif plottype == 'weight':
                    ax.semilogx(spec.get_freq(), spec.data, color=color,
                                zorder=20)
                else:
                    raise ValueError('Unknown plot type: %s' % plottype)
                # leg = ax.legend(('N', 'E', 'H'), 'lower right')

                if specnoise_st:
                    if spec.stats.channel[2] != 'S':
                        specid = spec.get_id()
                        try:
                            sp_noise = specnoise_st.select(id=specid)[0]
                        except IndexError:
                            continue
                        orientation = sp_noise.stats.channel[2]
                        if orientation in ['Z', '1']:
                            color = 'purple'
                        if orientation in ['N', '2']:
                            color = 'green'
                        if orientation in ['E', '3']:
                            color = 'blue'
                        if orientation == 'H':
                            color = 'red'
                        ax.loglog(sp_noise.get_freq(), sp_noise.data,
                                  linestyle=':', linewidth=2.,
                                  color=color, zorder=20)

                if not ax_text:
                    ax_text = '%s %s' % (spec.id[0:-1], spec.stats.instrtype)
                    if stack_plots:
                        text_y = 0.05 + (plotn-1) * 0.05
                    else:
                        text_y = 0.15
                        color = 'black'
                        ax_text += '\n%.1f km (%.1f km)' % (
                                                        spec.stats.hypo_dist,
                                                        spec.stats.epi_dist)
                    ax.text(0.05, text_y, ax_text,
                            horizontalalignment='left',
                            verticalalignment='bottom',
                            color=color,
                            transform=ax.transAxes,
                            zorder=50,
                            path_effects=path_effects)
                    ax_text = True

                if orientation == 'S':
                    if stack_plots:
                        text_y2 = text_y - 0.02
                    else:
                        text_y2 = 0.04
                        color = 'black'
                    fc = spec.stats.par['fc']
                    fc_err = spec.stats.par_err['fc']
                    Mw = spec.stats.par['Mw']
                    Mw_err = spec.stats.par_err['Mw']
                    Mo = mag_to_moment(Mw)
                    t_star = spec.stats.par['t_star']
                    t_star_err = spec.stats.par_err['t_star']
                    ax.text(0.05, text_y2,
                            'Mo: %.2g Mw: %.2f±%.2f\n'
                            'fc: %.2f±%.2f Hz t*: %.2f±%.2fs' %
                            (Mo, Mw, Mw_err, fc, fc_err, t_star, t_star_err),
                            horizontalalignment='left',
                            verticalalignment='bottom',
                            color=color,
                            fontsize=9,
                            transform=ax.transAxes,
                            zorder=50,
                            path_effects=path_effects)

    # Show the x-labels only for the last row
    for ax, ax2 in axes[-ncols:]:
        [t.set_visible(True) for t in ax.get_xticklabels()]
        ax.set_xlabel('Frequency (Hz)')
    # Show the y-labels only for the first column
    for i in range(0, len(axes)+ncols, ncols):
        try:
            ax = axes[i][0]
        except IndexError:
            continue
        try:
            ax2 = axes[i-1][1]
        except IndexError:
            continue
        [t.set_visible(True) for t in ax.get_yticklabels()]
        ax.set_ylabel('Weight')
        if plottype != 'weight':
            ax.set_ylabel('Seismic moment (Nm)')
            if ax2:
                [t.set_visible(True) for t in ax2.get_yticklabels()]
                ax2.set_ylabel('Magnitude')

    if config.PLOT_SHOW:
        plt.show()
    if config.PLOT_SAVE:
        # TODO: improve this:
        evid = spec_st.traces[0].stats.hypo.evid
        if plottype == 'regular':
            suffix = '.ssp.'
            message = 'Spectral'
        elif plottype == 'noise':
            suffix = '.sspnoise.'
            message = 'Noise'
        elif plottype == 'weight':
            suffix = '.sspweight.'
            message = 'Weight'
        figurefile = os.path.join(config.options.outdir, evid + suffix +
                                  config.PLOT_SAVE_FORMAT)
        if config.PLOT_SHOW:
            fig.savefig(figurefile, bbox_inches='tight')
        else:
            canvas = FigureCanvasAgg(fig)
            if async_plotter is not None:
                async_plotter.save(canvas, figurefile, bbox_inches='tight')
            else:
                canvas.print_figure(figurefile, bbox_inches='tight')
        logging.info(message + ' plots saved to: ' + figurefile)
    if not config.PLOT_SHOW:
        fig.clf()
