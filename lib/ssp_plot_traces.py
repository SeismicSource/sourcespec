# -*- coding: utf-8 -*-
# ssp_plot_traces.py
#
# (c) 2015 Claudio Satriano <satriano@ipgp.fr>
'''
Trace plotting routine.
'''
from __future__ import division
import os
import math
import logging
from ssp_setup import unload_matplotlib


def plot_traces(config, st, ncols=4, stack_plots=False):
    '''
    Plot displacement traces.
    Display to screen and/or save to file.
    '''
    # Unload matplotlib modules (which have been loaded by obspy.signal).
    unload_matplotlib()
    # Check config, if we need to plot at all
    if config.PLOT_SHOW == False and config.PLOT_SAVE == False:
        return
    # Re-import matplotlib
    import matplotlib
    matplotlib.rcParams['pdf.fonttype'] = 42 #to edit text in Illustrator
    # If we do not need to show the plot, we use the 'agg' backend
    if config.PLOT_SHOW == False:
        matplotlib.use('agg')
    # Finally we load the pyplot module
    import matplotlib.pyplot as plt

    # Determine the number of plots and axes min and max:
    nplots=0
    for station in set(x.stats.station for x in st.traces):
        st_sel = st.select(station=station)
        for instrtype in set(x.stats.instrtype for x in st_sel):
            nplots += 1
    nlines = int(math.ceil(nplots/ncols))

    # OK, now we can plot!
    if nlines <= 3 or stack_plots:
        figsize=(16,9)
    else:
        figsize=(16,18)
    fig = plt.figure(figsize=figsize)
    fig.subplots_adjust(hspace = .025, wspace = .03)

    # Plot!
    axes=[]
    plotn = 0
    for station in sorted(set(x.stats.station for x in st.traces)):
        st_sel = st.select(station=station)
        for instrtype in set(x.stats.instrtype for x in st_sel):
            plotn += 1
            ax_text = False

            if plotn == 1:
                if stack_plots:
                    ax = fig.add_subplot(1, 1, 1)
                else:
                    ax = fig.add_subplot(nlines, ncols, plotn)
            else:
                if not stack_plots:
                    ax = fig.add_subplot(nlines, ncols, plotn, sharex=axes[0], sharey=axes[0])
            ax.grid(True, which='both', linestyle='solid', color='#DDDDDD', zorder=0)
            ax.set_axisbelow(True)
            plt.setp(ax.get_xticklabels(), visible=False)
            plt.setp(ax.get_yticklabels(), visible=False)
            ax.tick_params(width=2) #FIXME: ticks are below grid lines!
            ax.ticklabel_format(style='scientific', axis='y', scilimits=(-3,3))
            axes.append(ax)

            for trace in st_sel.traces:
                if trace.stats.instrtype != instrtype:
                    continue
                if trace.stats.channel == 'Synth':
                    orientation = 'Synth'
                else:
                    orientation = trace.stats.channel[-1]
                if orientation == 'Z':
                    color='purple'
                if orientation == 'N':
                    color='green'
                if orientation == 'E':
                    color='blue'
                ax.plot(trace.times(), trace.data, color=color, zorder=20)

                if not ax_text:
                    if stack_plots:
                        text_y = 0.05 + (plotn-1) * 0.05
                    else:
                        text_y = 0.1
                        color = 'black'
                    ax.text(0.05, text_y, '%s %s' % (trace.stats.station, trace.stats.instrtype),
                            horizontalalignment='left',
                            verticalalignment='bottom',
                            color = color,
                            #backgroundcolor = (1, 1, 1, 0.7),
                            transform = ax.transAxes,
                            zorder = 50)
                    ax_text = True

                if orientation == 'Synth':
                    if stack_plots:
                        pass
                        #text_y2 = text_y - 0.02
                    else:
                        #text_y2 = 0.04
                        color = 'black'
                    #ax.text(0.05, text_y2, 'Mo: %.2g Mw: %.1f fc: %.2fHz t*: %.2fs' % (Mo, Mw, fc, t_star),
                    #        horizontalalignment='left',
                    #        verticalalignment='bottom',
                    #        #backgroundcolor = (1, 1, 1, 0.7), #FIXME: does not work in interactive plots
                    #        color = color,
                    #        fontsize = 9,
                    #        transform = ax.transAxes,
                    #        zorder = 50)

    # Show the x-labels only for the last row
    for ax in axes[-ncols:]:
        plt.setp(ax.get_xticklabels(), visible=True)
        ax.set_xlabel('Time (s)')
    # Show the y-labels only for the first column
    for i in range(0, len(axes)+ncols, ncols):
        try:
            ax = axes[i]
        except IndexError:
            continue
        plt.setp(ax.get_yticklabels(), visible=True)
        ax.set_ylabel('Displacement (m)')

    if config.PLOT_SHOW:
        plt.show()
    if config.PLOT_SAVE:
        #TODO: improve this:
        evid = st.traces[0].stats.hypo.evid
        figurefile = os.path.join(config.options.outdir, evid +\
                                  config.PLOT_SAVE_FORMAT)
        fig.savefig(figurefile, bbox_inches='tight')
        logging.info('Trace plots saved to: ' + figurefile)
