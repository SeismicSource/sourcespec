#!/usr/bin/env python
# -*- coding: utf8 -*- 
from obspy.core import Stream
from copy import deepcopy
from lib.ssp_setup import configure
from lib.ssp_read_traces import read_traces
from lib.ssp_process_traces import process_traces
from lib.ssp_build_spectra import build_spectra
from lib.ssp_spectral_model import spectral_model
from lib.ssp_plot_spectra import *
from lib.spectrum import Spectrum

def main():
    config = configure('source_model')
    fdelta = 0.01
    fmin = config.options.fmin
    fmax = config.options.fmax + fdelta
    fc = config.options.fc
    mag = config.options.mag
    t_star = config.options.t_star

    if len(config.args) > 0:
        st = read_traces(config)
        # Deconvolve, filter, cut traces:
        proc_st = process_traces(config, st)
        # Build spectra (amplitude in magnitude units)
        spec_st = build_spectra(config, proc_st)

        for tspec in spec_st:
            orientation = tspec.stats.channel[-1]
            if orientation != 'H': continue

            spec = Spectrum()
            spec.stats = deepcopy(tspec.stats)
            spec.stats.channel = 'Synth'
    
            freq = spec.get_freq()
            spec.data = spectral_model(freq, mag, fc, t_star)
    
            spec_st.append(spec)
    else:
        spec_st = Stream()
    
        spec = Spectrum()
        spec.stats.begin = fmin
        spec.stats.delta = fdelta
        spec.stats.npts = int((fmax-fmin)/fdelta)
        spec.stats.instrtype = 'Synth'
        spec.stats.channel = 'Synth'
    
        freq = spec.get_freq()
        spec.data = spectral_model(freq, mag, fc, t_star)
    
        spec_st.append(spec)
    
    class C(): pass
    config = C()
    config.PLOT_SHOW = True
    config.PLOT_SAVE = False
    
    plot_spectra(config, spec_st, ncols=1)

if __name__ == '__main__':
    main()
