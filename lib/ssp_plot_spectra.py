# -*- coding: utf-8 -*-
# ssp_plot_spectra.py
#
# Spectral plotting for source_spec
# (c) 2012 Claudio Satriano <satriano@ipgp.fr>
# (c) 2013 Claudio Satriano <satriano@ipgp.fr>,
#          Emanuela Matrullo <matrullo@geologie.ens.fr>
from __future__ import division
import sys
import math
import logging
#from matplotlib.ticker import MaxNLocator
from ssp_util import spec_minmax, moment_to_mag

synth_colors = [
    '#201F1F', 
    '#94F75B', 
    '#3EC2AA', 
    'FECC38',
    '#FC4384', 
]

def plot_spectra(config, spec_st, specnoise_st=None, ncols=4, stack_plots=False):
    # Unload matplotlib modules (which have been presumably loaded by
    # ObsPy).
    # Source:
    #   http://stackoverflow.com/questions/3285193/how-to-switch-backends-in-matplotlib-python
    modules = []
    for module in sys.modules:
        if module.startswith('matplotlib'):
            modules.append(module)
    for module in modules:
        sys.modules.pop(module)
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

    # OK, now we can plot!
    fig = plt.figure(figsize=(16,9))
    fig.subplots_adjust(hspace = .025, wspace = .03)

    # Determine the number of plots and axes min and max:
    nplots=0
    moment_minmax = None
    mag_minmax = None
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
            mag_minmax, dum =\
                spec_minmax(spec.data_mag, spec.get_freq(),
                            mag_minmax, freq_minmax) 
        for instrtype in set(x.stats.instrtype for x in spec_st_sel):
            nplots += 1
    nlines = int(math.ceil(nplots/ncols))
    moment_minmax[1] *= 10
    mag_minmax = moment_to_mag(moment_minmax)

    # Plot!
    plotn=1
    axes=[]
    for station in sorted(set(x.stats.station for x in spec_st.traces)):
        spec_st_sel = spec_st.select(station=station)
        for instrtype in set(x.stats.instrtype for x in spec_st_sel):
            ax_text = False
            if plotn==1:
                if stack_plots:
                    ax = fig.add_subplot(1, 1, 1)
                else:
                    ax = fig.add_subplot(nlines, ncols, plotn)
            else:
                if not stack_plots:
                    ax = fig.add_subplot(nlines, ncols, plotn, sharex=axes[0][0], sharey=axes[0][0])
            ax.set_xlim(freq_minmax)
            ax.set_ylim(moment_minmax)
            ax.grid(True, which='both')
            plt.setp(ax.get_xticklabels(), visible=False)
            plt.setp(ax.get_yticklabels(), visible=False)
            ax.tick_params(width=2)
            ax2 = ax.twinx()
            ax2.set_ylim(mag_minmax)
            plt.setp(ax2.get_xticklabels(), visible=False)
            plt.setp(ax2.get_yticklabels(), visible=True)
            for tick in ax2.yaxis.get_major_ticks():
                tick.set_pad(-1)
                tick.label2.set_horizontalalignment('right')
            ax2.yaxis.set_tick_params(width=0)
            axes.append((ax, ax2))
            for spec in spec_st_sel.traces:
                amp = None
                if spec.stats.instrtype != instrtype: continue
                if spec.stats.channel == 'Synth':
                    orientation = 'Synth'
                    amp = spec.data[0]
                else: orientation = spec.stats.channel[-1]
                if orientation == 'Z': color='purple'
                if orientation == 'N': color='green'
                if orientation == 'E': color='blue'
                if orientation == 'H': color='red'
                if orientation == 'Synth':
                    if stack_plots:
                        color = synth_colors[(plotn-1)%len(synth_colors)]
                    else:
                        color='black'
                ax.loglog(spec.get_freq(), spec.data, color=color)
                #ax2.semilogx(spec.get_freq(), spec.data_mag, color=color)
                #leg = ax.legend(('N', 'E', 'H'),
                #    'lower right')

                if specnoise_st:
                    if (spec.stats.channel == 'Synth' or
                            spec.stats.channel == 'H'):
                        continue
                    specid = spec.getId()
                    print specid
                    sp_noise = specnoise_st.select(id=specid)[0]
                    orientation = sp_noise.stats.channel[-1]
                    if orientation == 'Z': color='purple'
                    if orientation == 'N': color='green'
                    if orientation == 'E': color='blue'
                    ax.loglog(sp_noise.get_freq(), sp_noise.data, '--', color=color)

                if not ax_text:
                    if stack_plots:
                        text_y = 0.1 + (plotn-1) * 0.05
                    else:
                        text_y = 0.1
                    if not stack_plots:
                        color = 'black'
                    ax.text(0.05, text_y, '%s %s' % (spec.stats.station, spec.stats.instrtype),
                            horizontalalignment='left',
                            verticalalignment='bottom',
                            color = color,
                            transform = ax.transAxes)
                    ax_text = True
                if amp:
                    text_y = 0.05
                    ax.text(0.05, text_y, '%g' % (amp),
                            horizontalalignment='left',
                            verticalalignment='bottom',
                            color = 'black',
                            transform = ax.transAxes)
            plotn+=1

    # Show the x-labels only for the last row
    for ax, ax2 in axes[-ncols:]:
        plt.setp(ax.get_xticklabels(), visible=True)
        ax.set_xlabel('Frequency (Hz)')
    # Show the y-labels only for the first column
    for i in range(0, len(axes)+ncols, ncols):
        try:
            ax, dum = axes[i]
        except IndexError:
            continue
        try:
            dum, ax2 = axes[i-1]
        except IndexError:
            continue
        plt.setp(ax.get_yticklabels(), visible=True)
        plt.setp(ax2.get_yticklabels(), visible=True)
        ax.set_ylabel('Seismic moment (Nm)')
        ax2.set_ylabel('Magnitude')
    
    if config.PLOT_SHOW:
        plt.show()
    if config.PLOT_SAVE:
        #TODO: improve this:
        evid = spec_st.traces[0].stats.hypo.evid
        figurefile = config.options.outdir + '/' + evid + '.ssp.' +\
            config.PLOT_SAVE_FORMAT
        fig.savefig(figurefile, bbox_inches='tight')
        logging.info('Spectral plots saved to: ' + figurefile)

