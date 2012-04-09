# -*- coding: utf-8 -*-
# ssp_plot_spectra.py
#
# Spectral plotting for source_spec
# (c) 2012 Claudio Satriano <satriano@ipgp.fr>
from __future__ import division
import sys
import math
import logging
#from matplotlib.ticker import MaxNLocator
from ssp_util import spec_minmax

def plot_spectra(config, spec_st):
	# Unload matplotlib modules (which have been presumably loaded by
	# ObsPy.
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
	amp_minmax  = None
	freq_minmax = None
	nplots=0
	for station in set(x.stats.station for x in spec_st.traces):
		spec_st_sel = spec_st.select(station=station)
		for spec in spec_st_sel.traces:
			amp_minmax, freq_minmax =\
				spec_minmax(spec.data, spec.get_freq(), amp_minmax, freq_minmax) 
		for instrtype in set(x.stats.instrtype for x in spec_st_sel):
			nplots += 1
	ncols = 4
	nlines = int(math.ceil(nplots/ncols))
	plotn=1
	axes=[]

	for station in sorted(set(x.stats.station for x in spec_st.traces)):
		spec_st_sel = spec_st.select(station=station)
		for instrtype in set(x.stats.instrtype for x in spec_st_sel):
			if plotn==1:
				ax = fig.add_subplot(nlines, ncols, plotn)
			else:
				ax = fig.add_subplot(nlines, ncols, plotn, sharex=axes[0], sharey=axes[0])
			axes.append(ax)
			ax.set_xlim(freq_minmax)
			ax.set_ylim(amp_minmax)
			#ax.yaxis.set_major_locator(MaxNLocator(1))
			ax.grid(True, which='both')
			plt.setp(ax.get_xticklabels(), visible=False)
			plt.setp(ax.get_yticklabels(), visible=False)
			plotn+=1
			for spec in spec_st_sel.traces:
				if spec.stats.instrtype != instrtype: continue
				if spec.stats.channel == 'Synth': orientation = 'Synth'
				else: orientation = spec.stats.channel[-1]
				if orientation == 'N': color='green'
				if orientation == 'E': color='blue'
				if orientation == 'H': color='red'
				if orientation == 'Synth': color='black'
				ax.semilogx(spec.get_freq(), spec.data, color=color)
				ax.text(0.05, 0.1, '%s %s' % (spec.stats.station, spec.stats.instrtype),
				        horizontalalignment='left',
					verticalalignment='bottom',
					transform = ax.transAxes)

	# Show the x-labels only for the last row
	for ax in axes[-ncols:]:
		plt.setp(ax.get_xticklabels(), visible=True)
		ax.set_xlabel('Frequency (Hz)')
	# Show the y-labels only for the first column
	for i in range(0, len(axes)+ncols, ncols):
		try: ax = axes[i]
		except IndexError: continue
		plt.setp(ax.get_yticklabels(), visible=True)
		ax.set_ylabel('Magnitude')
	
	if config.PLOT_SHOW:
		plt.show()
	if config.PLOT_SAVE:
		#TODO: improve this:
		evid = spec_st.traces[0].stats.hypo.evid
		figurefile = config.options.outdir + '/' + evid + '.ssp.' +\
			config.PLOT_SAVE_FORMAT
		fig.savefig(figurefile, bbox_inches='tight')
		logging.info('Spectral plots saved to: ' + figurefile)

