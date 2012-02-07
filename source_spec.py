#!/usr/bin/env python
# -*- coding: utf-8 -*-

# source_spec.py
# Python code to invert S-wave displacement spectra
# Derived from sspec_v1.0.sh by Aldo Zollo and Claudio Satriano
# 
# (c) 2011-2012 Claudio Satriano <satriano@ipgp.fr>
# v 0.2 - 2012-02-06 - Extended and generalized for the CRL application 
# v 0.1 - 2012-01-17 - Initial Python port

# Computes the S-wave displacement spectra from stations recording a single event. 
#
# Reads as input a file, a tgz archive or a directory (which can, in turn, contain
# files and/or tgz archives) with traces in any format supported by ObsPy.
# Optionally, one can specify:
# - a file or a dir path containing station dataless
# - a hypocenter file
# - a phase file with P and S arrivals

# The code computes spectra of the two horizontal components and then modulus as:
#     sqrt(c1(w)^2+c2(w)^2)
# It then inverts spectra for a 3-parameter model (Mw,Fc,t*=T/Qs) using initial
# values for Mw, fc and t*:
#      log S(w)= log(coeff*Mo) + log((1/(1+(w/wc)^2)) + log (exp (- w *t_star/2)) 
# Plots all spectra on a single log-log graph (Mw vs log-frequency)
# Plot obs vs theo spectra for each vectorial component
# Computes average and st.dev of Mw, Mo, fc, source radius and Brune sd
#
# Original csh version by Aldo Zollo <zollo@unina.it>
# Bash version by Claudio Satriano <satriano@na.infn.it>
# v 1.0 - 2008-10-29 - Ported to bash
#                    - Added more handling of exceptions to run safely in
#		       automatic mode
#		       Claudio Satriano <satriano@na.infn.it>
from __future__ import division
import sys
import os
from optparse import OptionParser
from imp import load_source
import math
import numpy as np
from scipy.optimize import curve_fit
from obspy.core import Stream
from ssp_read_traces import read_traces
from ssp_util import *
from ssp_plot_spectra import *
import spectrum

try:
        # ipython >= 0.11
        from IPython.frontend.terminal.embed import InteractiveShellEmbed
        ipshell = InteractiveShellEmbed()
except ImportError:
        # ipython < 0.11
        from IPython.Shell import IPShellEmbed
        ipshell = IPShellEmbed()


def dprint(string):
	if DEBUG:
		sys.stderr.write(string)
		sys.stderr.write('\n')


def main():
	global DEBUG
	usage = "usage: %prog [options] trace_file(s) | trace_dir"

	parser = OptionParser(usage=usage);
	parser.add_option("-c", "--configfile", dest="config_file", action="store", default=None,
                  help="Load configuration from FILE", metavar="DIR | FILE")
	parser.add_option("-d", "--dataless", dest="dataless", action="store", default=None,
                  help="Search for dataless in DIR or in FILE", metavar="DIR | FILE")
	parser.add_option("-H", "--hypocenter", dest="hypo_file", action="store", default=None,
                  help="Get hypocenter information from FILE", metavar="FILE")
	parser.add_option("-p", "--pickfile", dest="pick_file", action="store", default=None,
                  help="Get picks from FILE", metavar="FILE")
	parser.add_option("-o", "--outdir", dest="outdir", action="store", default=None,
                  help="Save output to OUTDIR", metavar="OUTDIR")

	(options, args) = parser.parse_args();
	if len(args) < 1:
		parser.print_usage(file=sys.stderr)
		sys.stderr.write("\tUse '-h' for help\n\n")
		sys.exit(1)

	st = read_traces(args, options)
	#for trace in st.traces:
	#	print trace.getId(), trace.stats.paz
	#sys.exit()

	if options.config_file == None:
		options.config_file = 'config.py'
	try:
		config=load_source("config", options.config_file)
		DEBUG = config.DEBUG
        except:
                sys.stderr.write("Unable to open file: %s\n" % options.config_file)
                sys.exit(1)
	
	# Loop on stations for building spectra and the spec_st object 
	spec_st = Stream()
	for trace in st.traces:
		# Remove instrument response
		if remove_instr_response(trace,
				config.correct_sensitivity_only, config.pre_filt) == None:
			dprint('Undefined instrument response: continue')
			continue

		stats = trace.stats
		comp  = stats.channel
		# skip vertical components
		if comp[-1] == 'Z':
			continue
		station = stats.station
		dprint('%s %s' % (station, comp))

		# compute hypocentral distance hd
		hd = hypo_dist(trace)
		if hd == None: continue
		hd_m = hd*1000

		# S time window
		s_arrival_time = swave_arrival(trace, config.vs)
		t1 = s_arrival_time - config.pre_s_time
		t2 = t1 + config.s_win_length
		trace_cut = trace.slice(t1, t2)
		npts = len(trace_cut.data)
		if npts == 0:
			dprint('No data for the select cut interval: continue')
			continue
		
		# TODO: parameterize
		# coefficient for converting displ spectrum
		# to seismic moment (Aki&Richards,1980)
		vs_m = config.vs*1000
		vs3  = pow(vs_m,3) 
		coeff = 4 * math.pi * vs3 * config.rho /config.rps/2.
		dprint('coeff= %f' % coeff)

		# compute S-wave (displacement) spectra from
		# accelerometers and velocimeters, uncorrected for attenuation,
		# corrected for instrumental constants, normalized by
		# hypocentral distance
		instrtype = stats.instrtype
		if instrtype == 'acc':
			nint = 2 #number of intergrations to perform
			# band-pass frequencies:
			# TODO: calculate from sampling rate?
			bp_freqmin = config.bp_freqmin_acc
			bp_freqmax = config.bp_freqmax_acc
			# cut frequencies:
			freq1 = config.freq1_acc
			freq2 = config.freq2_acc
		elif instrtype == 'vel':
			nint = 1 #number of intergrations to perform
			# band-pass frequencies:
			# TODO: calculate from sampling rate?
			bp_freqmin = config.bp_freqmin_vel
			bp_freqmax = config.bp_freqmax_vel
			# cut frequencies:
			freq1 = config.freq1_vel
			freq2 = config.freq2_vel
		else: continue

		# remove the mean...
		trace_cut.detrend(type='constant')
		# ...and the linear trend...
		trace_cut.detrend(type='linear')
		trace_cut.filter(type='bandpass', freqmin=bp_freqmin, freqmax=bp_freqmax)
		# ...and taper
		cosine_taper(trace_cut.data, width=0.5)

		# normalization for the hypocentral distance
		trace_cut.data *= hd_m

		# calculate fft
		spec = spectrum.do_spectrum(trace_cut)
		spec.stats.instrtype = instrtype
		spec.stats.coords    = stats.coords
		spec.stats.hypo      = stats.hypo
		spec.stats.hypo_dist = stats.hypo_dist

		# Integrate in frequency domain (divide by the pulsation omega)
		for i in range(0,nint):
			spec.data /= (2 * math.pi * spec.get_freq())

		# TODO: konno-omachi
		# smooth the abs of fft
		data_smooth = smooth(spec.data, 6)

		# Uncomment these lines to see the effect of smoothing
		#plt.figure()
		#plt.loglog(spec.get_freq(), spec.data, color='gray')
		#plt.loglog(spec.get_freq(), data_smooth)
		#plt.grid(True)
		#plt.show()
		#sys.exit()

		spec.data = data_smooth
		# convert to seismic moment
		spec.data *= coeff

		# Cut the spectrum between freq1 and freq2
		spec_cut = spec.slice(freq1, freq2)

		spec_st.append(spec_cut)

	#end of loop on stations for building spectra

	# Add to spec_st the "horizontal" component, obtained from the
	# modulus of the N-S and E-W components.
	for station in set(x.stats.station for x in spec_st.traces):
		spec_st_sel = spec_st.select(station=station)
		for instrtype in set(x.stats.instrtype for x in spec_st_sel):
			spec_h = None
			for spec in spec_st_sel.traces:
				if spec.stats.instrtype != instrtype: continue
				if spec_h == None:
					spec_h = spec.copy()
					spec_h.stats.channel = 'H'
				else:
					data_h = spec_h.data
					data   = spec.data
					data_h = np.power(data_h, 2) + np.power(data, 2)
					data_h = np.sqrt(data_h / 2) #divide by the number of horizontal components
					spec_h.data = data_h
			spec_st.append(spec_h)
			
	# convert the spectral amplitudes to moment magnitude
	for spec in spec_st.traces:
		spec.data = (np.log10(spec.data) - 9.1 ) / 1.5 

	# Inversion of displacement spectra
	loge = math.log10(math.e)
	# Spectral weighting:
	#   weight for f<=f_weight
	#   1      for f> f_weight
	f_weight = config.f_weight
	weight   = config.weight
	sourcepar = dict()
	for station in set(x.stats.station for x in spec_st.traces):
		spec_st_sel = spec_st.select(station=station)
		for spec in spec_st_sel.traces:
			if spec.stats.channel != 'H': continue
			dprint(station)

			# spectral amplitude is in Mw units
			amp = spec.data

			# We calculate the initial value for Mw,
			# as an average of the first 5 values of amp
			Mw_0 = amp[0:5].mean()

			# We try to retrieve fc_0 from the configuration...
			try: fc_0 = config.fc_0
			# ...if it is not available, we calculate it
			except: 
				l_m0 = Mw_0 * 1.5 + 9.1
				l_beta = math.log10(vs_m)
				l_bsd = math.log10(config.bsd)
				l_fc = l_bsd - l_m0 + 3*l_beta - 0.935
				l_fc /= 3.
				fc_0 = math.pow(10, l_fc)

			# initial value for t_star
			t_star_0 = config.t_star_0
			
			# print initial values of fc, M0 and t_star
			dprint('INITIAL CORNER FREQUENCY= %f' % fc_0)
			dprint('INITIAL MOMENT MAGNITUDE= %f' % Mw_0)
			dprint('INITIAL T_STAR= %f' % t_star_0)

			hd   = spec.stats.hypo_dist
			hd_m = hd*1000
			az   = 0 #TODO: add azimuth computation
			dprint('%s %s %f %f' % (station, spec.stats.instrtype, hd, az))
			coeff = math.pi*loge*hd_m/vs_m
			dprint('coeff= %f' % coeff)

			def spectral_model(freq, Mw, fc, t_star):
				# log S(w)= log(coeff*Mo) + log((1/(1+(w/wc)^2)) + log (exp (- w *t_star/2)) 
				# attenuation model: exp[-pi t* f] with t*=T /Q
				return Mw - np.log10(1. + np.power((freq/fc),2) ) - loge*(math.pi*t_star*freq) 

			params_name = ('Mw', 'fc', 't_star')
			params_0 = np.array([Mw_0, fc_0, t_star_0])

			xdata = spec.get_freq()
			ydata = amp
			yerr = np.ones(len(ydata))
			# 'curve_fit' interprets 'yerr' as standard deviation vector and calculates
			# weights as 1/yerr^2 . Therefore we build yerr as:
			yerr[xdata<=f_weight] = 1./math.sqrt(weight)
			# Curve fitting using the Levenburg-Marquardt algorithm
			params_opt, params_cov = curve_fit(spectral_model, xdata, ydata, p0=params_0, sigma=yerr)
			#print params_opt, params_cov
			par = dict(zip(params_name, params_opt))
			par['hyp_dist'] = hd
			par['az'] = az
			chanId = '%s.%s' % (station, spec.stats.instrtype)
			sourcepar[chanId] = par

			spec_synth = spec.copy()
			spec_synth.stats.channel = 'Synth'
			spec_synth.data = spectral_model(xdata, *params_opt)
			spec_st.append(spec_synth)


	# Filter stations with negative t_star
	# or with anomalous corner frequencies
	f1 = config.min_corner_freq
	f2 = config.max_corner_freq
	for statId in sourcepar.keys():
		par = sourcepar[statId]
		t_star = par['t_star']
		if t_star < 0:
			print 'Ignoring station: %s t_star: %f' % (statId, t_star)
			sourcepar.pop(statId, None)
		fc = par['fc']
		if fc < f1 or fc > f2:
			print 'Ignoring station: %s fc: %f' % (statId, fc)
			sourcepar.pop(statId, None)

	if len(sourcepar) == 0: sys.exit()
	# Print source parameters
	print '*** Source parameters ***'
	for statId in sorted(sourcepar.keys()):
		par = sourcepar[statId]
		sys.stdout.write('%s ' % statId)
		for key in par:
			sys.stdout.write('  %s %6.3f ' % (key, par[key]))
		sys.stdout.write('\n')


	## Compute average source parameters
	print '*** Average source parameters ***'
	# Mw
	Mw_array = np.array(list(x['Mw'] for x in sourcepar.values()))
	Mw_mean  = Mw_array.mean()
	Mw_std   = Mw_array.std()
	print 'Mw: %.2f +/- %.2f' % (Mw_mean, Mw_std)

	# Mo (N.m)
	Mo_array = np.power(10, Mw_array*1.5 + 9.1)
	Mo_mean  = Mo_array.mean()
	Mo_std   = Mo_array.std()
	print 'Mo: %.3e +/- %.3e N.m' % (Mo_mean, Mo_std)

	# fc , hertz
	fc_array = np.array(list(x['fc'] for x in sourcepar.values()))
	fc_mean  = fc_array.mean()
	fc_std   = fc_array.std()
	print 'fc: %.3f +/- %.3f Hz' % (fc_mean, fc_std)

	# ra, radius (meters)
	ra_array = 0.37 * vs_m / fc_array
	ra_mean  = ra_array.mean()
	ra_std   = ra_array.std()
	print 'Source radius: %.3f +/- %.3f m' % (ra_mean, ra_std)

	# bds, Brune stress drop (MPa)
	bsd_array = 0.4375 * Mo_array / np.power(ra_array, 3) * 1e-6
	bsd_mean  = bsd_array.mean()
	bsd_std   = bsd_array.std()
	print 'Brune stress drop: %.3f +/- %.3f MPa' % (bsd_mean, bsd_std)


	# (Optional) Plotting
	if config.DOPLOTS:
		plotspectra(spec_st)

        # Remove the bytecoded version of the config file,
        # generated by "load_source"
        try: os.remove("%sc" % options.config_file)
        except: pass


if __name__ == '__main__':
	try:
		thismodule=os.path.basename(__file__).replace('.py','')
		this=__import__(thismodule)
		main()
	except KeyboardInterrupt:
		sys.exit(1)
