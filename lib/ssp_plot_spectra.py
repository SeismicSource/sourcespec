# -*- coding: utf-8 -*-
# ssp_plot_spectra.py
#
# (c) 2012 Claudio Satriano <satriano@ipgp.fr>
# (c) 2013-2014 Claudio Satriano <satriano@ipgp.fr>,
#               Emanuela Matrullo <matrullo@geologie.ens.fr>
# (c) 2015-2016 Claudio Satriano <satriano@ipgp.fr>
"""Spectral plotting routine."""
from __future__ import division
import os
import math
import logging
from ssp_util import spec_minmax, moment_to_mag, mag_to_moment
from ssp_setup import unload_matplotlib

synth_colors = [
    '#201F1F',
    '#94F75B',
    '#3EC2AA',
    '#FECC38',
    '#FC4384',
]


def plot_spectra(config, spec_st, specnoise_st=None, ncols=4,
                 stack_plots=False, plottype='regular'):
    """
    Plot spectra for signal and noise.

    Display to screen and/or save to file.
    """
    # Unload matplotlib modules (which have been loaded by obspy.signal).
    unload_matplotlib()
    # Check config, if we need to plot at all
    if not config.PLOT_SHOW and not config.PLOT_SAVE:
        return
    # Re-import matplotlib
    import matplotlib
    matplotlib.rcParams['pdf.fonttype'] = 42  # to edit text in Illustrator
    # If we do not need to show the plot, we use the 'agg' backend
    if not config.PLOT_SHOW:
        matplotlib.use('agg')
    # Finally we load the pyplot module
    import matplotlib.pyplot as plt

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
    fig = plt.figure(figsize=figsize)
    fig.subplots_adjust(hspace=.025, wspace=.03)

    # Plot!
    axes = []
    plotn = 0
    for station in sorted(set(x.stats.station for x in spec_st.traces)):
        spec_st_sel = spec_st.select(station=station)
        # 'code' is band+instrument code
        for code in set(x.stats.channel[0:2] for x in spec_st_sel):
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
            plt.setp(ax.get_xticklabels(), visible=False)
            plt.setp(ax.get_yticklabels(), visible=False)
            ax.tick_params(width=2)  # FIXME: ticks are below grid lines!

            # ax2 has magnitude units
            if plottype != 'weight':
                if ((stack_plots and plotn == 1) or not stack_plots):
                    ax2 = ax.twinx()
                    ax2.set_ylim(mag_minmax)
                    plt.setp(ax2.get_xticklabels(), visible=False)
                    plt.setp(ax2.get_yticklabels(), visible=True)
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
                if orientation in ['N', '2']:
                    color = 'green'
                if orientation in ['E', '3']:
                    color = 'blue'
                if orientation == 'H':
                    color = 'red'
                if orientation == 'S':
                    if stack_plots:
                        color = synth_colors[(plotn-1) % len(synth_colors)]
                    else:
                        color = 'black'
                if plottype in ['regular', 'noise']:
                    ax.loglog(spec.get_freq(), spec.data, color=color,
                              zorder=20)
                    if orientation == 'S':
                        ax.axvline(spec.stats.par['fc'], color='#999999',
                                   linewidth=2., zorder=1)
                elif plottype == 'weight':
                    ax.semilogx(spec.get_freq(), spec.data, color=color,
                                zorder=20)
                else:
                    raise ValueError('Unknown plot type: %s' % plottype)
                #leg = ax.legend(('N', 'E', 'H'),
                #    'lower right')

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
                    if stack_plots:
                        text_y = 0.05 + (plotn-1) * 0.05
                    else:
                        text_y = 0.1
                        color = 'black'
                    ax_text = '%s %s' % (spec.id[0:-1], spec.stats.instrtype)
                    ax.text(0.05, text_y, ax_text,
                            horizontalalignment='left',
                            verticalalignment='bottom',
                            color=color,
                            #backgroundcolor=(1, 1, 1, 0.7),
                            transform=ax.transAxes,
                            zorder=50)
                    ax_text = True

                if orientation == 'S':
                    if stack_plots:
                        text_y2 = text_y - 0.02
                    else:
                        text_y2 = 0.04
                        color = 'black'
                    fc = spec.stats.par['fc']
                    Mw = spec.stats.par['Mw']
                    Mo = mag_to_moment(Mw)
                    t_star = spec.stats.par['t_star']
                    ax.text(0.05, text_y2,
                            'Mo: %.2g Mw: %.1f fc: %.2fHz t*: %.2fs' %
                            (Mo, Mw, fc, t_star),
                            horizontalalignment='left',
                            verticalalignment='bottom',
                            #FIXME: does not work in interactive plots
                            #backgroundcolor = (1, 1, 1, 0.7),
                            color=color,
                            fontsize=9,
                            transform=ax.transAxes,
                            zorder=50)

    # Show the x-labels only for the last row
    for ax, ax2 in axes[-ncols:]:
        plt.setp(ax.get_xticklabels(), visible=True)
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
        plt.setp(ax.get_yticklabels(), visible=True)
        ax.set_ylabel('Weight')
        if plottype != 'weight':
            ax.set_ylabel('Seismic moment (Nm)')
            if ax2:
                plt.setp(ax2.get_yticklabels(), visible=True)
                ax2.set_ylabel('Magnitude')

    if config.PLOT_SHOW:
        plt.show()
    if config.PLOT_SAVE:
        #TODO: improve this:
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
        fig.savefig(figurefile, bbox_inches='tight')
        logging.info(message + ' plots saved to: ' + figurefile)
    fig.clf()
