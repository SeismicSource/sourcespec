# -*- coding: utf-8 -*-
"""
A Spectrum() class defined as a modification of the ObsPy class Trace().

Provides the high-level function do_spectrum()
and the low-level funciton do_fft().

:copyright:
    2012-2021 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement, Version 2.1
    (http://www.cecill.info/index.en.html)
"""
import numpy as np
from copy import copy, deepcopy
from obspy.core import Trace


def do_fft(signal, delta):
    """Compute the complex Fourier transform of a signal."""
    npts = len(signal)
    # if npts is even, we make it odd
    # so that we do not have a negative frequency in the last point
    # (see numpy.fft.rfft doc)
    if not npts % 2:
        npts -= 1

    fft = np.fft.rfft(signal, n=npts) * delta
    fftfreq = np.fft.fftfreq(len(signal), d=delta)
    fftfreq = fftfreq[0:fft.size]
    return fft, fftfreq


def do_spectrum(trace):
    """Compute the spectrum of an ObsPy Trace object."""
    signal = trace.data
    delta = trace.stats.delta
    amp, freq = do_fft(signal, delta)

    tr = Spectrum()
    # remove DC component (freq=0)
    tr.data = abs(amp)[1:]
    tr.stats.delta = 1. / (len(signal) * delta)
    tr.stats.begin = tr.stats.delta  # the first frequency is not 0!
    tr.stats.npts = len(tr.data)
    # copy some relevant header field
    tr.stats.station = trace.stats.station
    tr.stats.network = trace.stats.network
    tr.stats.location = trace.stats.location
    tr.stats.channel = trace.stats.channel
    return tr


class Spectrum(Trace):

    def get_freq(self):
        fdelta = self.stats.delta
        freq = np.arange(0, self.stats.npts*fdelta, fdelta)
        freq = freq[0:self.stats.npts]
        freq += self.stats.begin
        return freq

    def plot(self, **kwargs):
        freq = self.get_freq()
        import matplotlib.pyplot as plt
        plt.loglog(freq, self.data, **kwargs)
        plt.grid(True)
        plt.xlabel('Frequency (Hz)')
        plt.ylabel('Amplitude')
        plt.show()

    def slice(self, fmin, fmax, pad=False, nearest_sample=True,
              fill_value=None):
        t = self.stats.starttime
        freq = self.get_freq()
        begin = self.stats.begin
        spec_slice = copy(self)
        spec_slice.stats = deepcopy(self.stats)
        spec_slice.trim(t-begin+fmin, t-begin+fmax, pad, nearest_sample,
                        fill_value)
        delta_t = t - spec_slice.stats.starttime
        if delta_t > 0:
            spec_slice.stats.begin = begin - delta_t
        else:
            # find the closest frequency to fmin:
            idx = (np.abs(spec_slice.get_freq()-fmin)).argmin()
            spec_slice.stats.begin = freq[idx]
        return spec_slice
