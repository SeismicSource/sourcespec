#!/usr/bin/env python
# -*- coding: utf-8 -*-
# source_spec.py
#
# Main function for source_spec
# (c) 2012 Claudio Satriano <satriano@ipgp.fr>
# (c) 2013 Claudio Satriano <satriano@ipgp.fr>,
#          Emanuela Matrullo <matrullo@geologie.ens.fr>,
#          Agnes Chounet <chounet@ipgp.fr>
from __future__ import division
import multiprocessing
from lib.ssp_setup import configure, setup_logging, ssp_exit
from lib.ssp_read_traces import read_traces
from lib.ssp_process_traces import process_traces
from lib.ssp_build_spectra import build_spectra
from lib.ssp_local_magnitude import local_magnitude
from lib.ssp_inversion import spectral_inversion
from lib.ssp_output import write_output
from lib.ssp_residuals import spectral_residuals
from lib.ssp_plot_spectra import plot_spectra

def main():
    # Setup stage
    config = configure()
    setup_logging(config)
    st = read_traces(config)

    # Now that we (hopefully) have the evid
    # we rename the logfile to use the evid
    #TODO: improve this:
    evid = st.traces[0].stats.hypo.evid
    setup_logging(config, evid)

    # Deconvolve, filter, cut traces:
    proc_st, noise_st = process_traces(config, st)
    # Build spectra (amplitude in magnitude units)
    spec_st, specnoise_st, weight_st = build_spectra(config, proc_st, noise_st)

    #Ml = local_magnitude(config, proc_st)
    Ml = local_magnitude(config, st, deconvolve=True)

    # Spectral inversion
    sourcepar = spectral_inversion(config, spec_st, weight_st, Ml)

    # Save output
    sourcepar_mean = write_output(config, evid, sourcepar) 

    # Save residuals
    spectral_residuals(config, spec_st, evid, sourcepar_mean)

    # Plotting
    pool = multiprocessing.Pool()
    args = [
            (config, spec_st, specnoise_st, 'regular'),
            (config, specnoise_st, None, 'noise'),
            (config, weight_st, None, 'weight')
           ]
    pool.map(__call_plot_spectra__, args)

    #plot_spectra(config, spec_st, specnoise_st=specnoise_st)
    #plot_spectra(config, specnoise_st, plottype='noise')
    #plot_spectra(config, weight_st, plottype='weight')

    ssp_exit()

def __call_plot_spectra__(args):
    config, spec_st, specnoise_st, plottype = args
    plot_spectra(config, spec_st,
            specnoise_st=specnoise_st,
            plottype=plottype)

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        ssp_exit(1)
