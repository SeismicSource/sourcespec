# -*- coding: utf-8 -*-
# spectrum.py
# Introduces the class Spectrum as a modification of the ObsPy class Trace()
# Provides the high-level function do_spectrum() and the low-level funciton do_fft()
#
# v 0.1 - 2012-01-17 - Claudio Satriano <satriano@ipgp.fr>
import numpy as np
from obspy.core import Trace
import matplotlib.pyplot as plt

def do_fft(signal, delta):
	'''Computes the complex Fourier transform of a signal'''
	npts=len(signal)
	# if npts is even, we make it odd
	# so that we do not have a negative frequency
	# in the last point 
	# (see http://docs.scipy.org/doc/numpy/reference/generated/numpy.fft.rfft.html)
	if not npts%2: npts -= 1

	fft = np.fft.rfft(signal, n=npts) * delta
	fftfreq = np.fft.fftfreq(len(signal), d=delta)
	fftfreq = fftfreq[0:fft.size]

	return fft, fftfreq

def do_spectrum(trace):
	'''Computes the spectrum of an ObsPy Trace object'''
	signal = trace.data
	delta  = trace.stats.delta
	amp, freq = do_fft(signal, delta)

	tr = Spectrum()
	# remove DC component (freq=0)
	tr.data = abs(amp)[1:]
	tr.stats.delta    = 1. / (len(signal) * delta)
	tr.stats.begin    = tr.stats.delta #the first frrquency is not 0!
	# copy some relevant header field
	tr.stats.station  = trace.stats.station
	tr.stats.network  = trace.stats.network
	tr.stats.location = trace.stats.location
	tr.stats.channel  = trace.stats.channel

	return tr

class Spectrum(Trace):
	def get_freq(self):
 		fdelta = self.stats.delta
		freq = np.arange(0, len(self.data)*fdelta, fdelta) 
		freq = freq[0:len(self.data)]
		freq += self.stats.begin
		return freq
	def plot(self, **kwargs):
		freq = self.get_freq()
		plt.loglog(freq, self.data, **kwargs)
		plt.grid(True)
		plt.xlabel('Frequency (Hz)')
		plt.ylabel('Amplitude')
		plt.show()
	def slice(self, fmin, fmax):
		t = self.stats.starttime
		freq = self.get_freq()
		begin = self.stats.begin
		spec_slice = super(Spectrum, self).slice(t-begin+fmin, t-begin+fmax)
		#find the closest frequency to fmin:
		idx=(np.abs(freq-fmin)).argmin()
		spec_slice.stats.begin = freq[idx]
		return spec_slice
