# -*- coding: utf8 -*-
"""
Direct spectral modelling.

:copyright:
    2013-2014 Claudio Satriano <satriano@ipgp.fr>,
              Agnes Chounet <chounet@ipgp.fr>

    2015-2021 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement, Version 2.1
    (http://www.cecill.info/index.en.html)
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import math
import numpy as np
# import cPickle as pickle
from copy import deepcopy
from obspy.core import Stream
from sourcespec.spectrum import Spectrum
from sourcespec.ssp_setup import configure, ssp_exit
from sourcespec.ssp_read_traces import read_traces
from sourcespec.ssp_process_traces import process_traces
from sourcespec.ssp_build_spectra import build_spectra
from sourcespec.ssp_spectral_model import spectral_model, objective_func
from sourcespec.ssp_util import mag_to_moment, moment_to_mag
from sourcespec.ssp_plot_spectra import plot_spectra
from sourcespec.ssp_plot_traces import plot_traces


def make_synth(config, spec_st, trace_spec=None):
    fdelta = 0.01
    fmin = config.options.fmin
    fmax = config.options.fmax + fdelta

    residuals = list()
    for fc, mag, Mo, t_star, alpha in zip(
            config.options.fc, config.options.mag, config.options.Mo,
            config.options.t_star, config.options.alpha):
        spec = Spectrum()
        if trace_spec:
            spec.stats = deepcopy(trace_spec.stats)
        else:
            spec.stats.begin = fmin
            spec.stats.delta = fdelta
            spec.stats.npts = int((fmax-fmin)/fdelta)

        if math.isnan(Mo):
            Mo = mag_to_moment(mag)
        else:
            mag = moment_to_mag(Mo)

        spec.stats.station = 'Mw: %.1f fc: %.2fHz t*: %.2fs alpha: %.2f' %\
            (mag, fc, t_star, alpha)
        spec.stats.instrtype = 'Synth'
        spec.stats.channel = spec.stats.channel[:-1] + 'S'
        spec.stats.par = {'Mw': mag, 'fc': fc, 't_star': t_star,
                          'alpha': alpha}

        freq = spec.get_freq()
        spec.data_mag = spectral_model(freq, mag, fc, t_star, alpha)
        spec.data = mag_to_moment(spec.data_mag)
        spec_st.append(spec)

        if trace_spec:
            objective_func2 = objective_func(freq, trace_spec.data_mag,
                                             np.ones_like(trace_spec.data_mag))
            print(Mo, mag, fc, t_star, alpha,
                  objective_func2((mag, fc, t_star, alpha)))
            residuals.append([Mo, mag, fc, t_star,
                             objective_func2((mag, fc, t_star, alpha))])

    # figurefile = config.options.station+'-'+config.options.evid+'-res.pickle'
    # fp = open(figurefile,'wb')
    # pickle.dump(residuals,fp)
    # fp.close()
    # print 'residuals caculated. exit code'


def main():
    config = configure('source_model')

    if len(config.options.trace_path) > 0:
        st = read_traces(config)
        # Deconvolve, filter, cut traces:
        proc_st = process_traces(config, st)
        # Build spectra (amplitude in magnitude units)
        spec_st = build_spectra(config, proc_st)
        if len(spec_st) == 0:
            ssp_exit()
        # We keep just horizontal component:
        spec_st = Stream([x for x in spec_st.traces
                          if (x.stats.channel[-1] == 'H')])

        for trace_spec in spec_st:
            orientation = trace_spec.stats.channel[-1]
            if orientation != 'H':
                continue
            make_synth(config, spec_st, trace_spec)
    else:
        spec_st = Stream()
        make_synth(config, spec_st)

    if config.options.plot:
        config.PLOT_SHOW = True
    else:
        config.PLOT_SHOW = False
    config.PLOT_SAVE = False

    plot_traces(config, proc_st, ncols=2, block=False)
    plot_spectra(config, spec_st, ncols=1, stack_plots=True)


if __name__ == '__main__':
    main()
