#!/usr/bin/env python
# -*- coding: utf8 -*- 
from copy import deepcopy
from obspy.core import Stream
from lib.ssp_setup import configure
from lib.ssp_read_traces import read_traces
from lib.ssp_process_traces import process_traces
from lib.ssp_build_spectra import build_spectra
from lib.ssp_spectral_model import spectral_model
from lib.ssp_plot_spectra import plot_spectra
from lib.spectrum import Spectrum

def make_synth(config, spec_st, trace_spec=None):
    fdelta = 0.01
    fmin = config.options.fmin
    fmax = config.options.fmax + fdelta
    for n, fc in enumerate(config.options.fc):
        mag = config.options.mag[n]
        t_star = config.options.t_star[n]

        spec = Spectrum()
        if trace_spec:
            spec.stats = deepcopy(trace_spec.stats)
        else:
            spec.stats.begin = fmin
            spec.stats.delta = fdelta
            spec.stats.npts = int((fmax-fmin)/fdelta)

        spec.stats.station = 'mag: %.1f fc: %.2fHz t*: %.2fs' % (mag, fc, t_star)
        spec.stats.instrtype = 'Synth'
        spec.stats.channel = 'Synth'
    
        freq = spec.get_freq()
        spec.data = spectral_model(freq, mag, fc, t_star)
        spec_st.append(spec)

def main():
    config = configure('source_model')

    if len(config.args) > 0:
        st = read_traces(config)
        # Deconvolve, filter, cut traces:
        proc_st = process_traces(config, st)
        # Build spectra (amplitude in magnitude units)
        spec_st = build_spectra(config, proc_st)
        # We keep just horizontal component:
        spec_st = Stream([ x for x in spec_st.traces if (x.stats.channel[-1] == 'H') ])

        for trace_spec in spec_st:
            orientation = trace_spec.stats.channel[-1]
            if orientation != 'H': continue
            make_synth(config, spec_st, trace_spec)
    else:
        spec_st = Stream()
        make_synth(config, spec_st)
    
    class C(): pass
    config = C()
    config.PLOT_SHOW = True
    config.PLOT_SAVE = False
    
    plot_spectra(config, spec_st, ncols=1, stack_plots=True)

if __name__ == '__main__':
    main()
