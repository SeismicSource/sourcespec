# -*- coding: utf-8 -*-
"""
Earthquake source parameters from inversion of S-wave spectra.

:copyright:
    2012 Claudio Satriano <satriano@ipgp.fr>

    2013-2014 Claudio Satriano <satriano@ipgp.fr>,
              Emanuela Matrullo <matrullo@geologie.ens.fr>,
              Agnes Chounet <chounet@ipgp.fr>

    2015-2021 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement, Version 2.1
    (http://www.cecill.info/index.en.html)
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from sourcespec.ssp_setup import (configure, setup_logging, save_config,
                                  init_plotting, ssp_exit)
from sourcespec.ssp_read_traces import read_traces
from sourcespec.ssp_process_traces import process_traces
from sourcespec.ssp_build_spectra import build_spectra
from sourcespec.ssp_local_magnitude import local_magnitude
from sourcespec.ssp_inversion import spectral_inversion
from sourcespec.ssp_output import write_output
from sourcespec.ssp_residuals import spectral_residuals
from sourcespec.ssp_plot_spectra import plot_spectra
from sourcespec.ssp_plot_traces import plot_traces
from sourcespec.ssp_html_report import html_report


def main():
    # Setup stage
    config = configure()
    setup_logging(config)

    st = read_traces(config)

    setup_logging(config, config.hypo.evid)

    # Save config to out dir
    save_config(config)

    # Deconvolve, filter, cut traces:
    proc_st = process_traces(config, st)

    # Build spectra (amplitude in magnitude units)
    spec_st, specnoise_st, weight_st =\
        build_spectra(config, proc_st, noise_weight=True)

    plotter = init_plotting(config)
    ntr = len(set(t.id[:-1] for t in proc_st))
    ncols = 4 if ntr > 6 else 3
    plot_traces(config, proc_st, ncols=ncols, async_plotter=plotter)

    # Spectral inversion
    sourcepar, sourcepar_err =\
        spectral_inversion(config, spec_st, weight_st)

    # Local magnitude
    if config.compute_local_magnitude:
        local_magnitude(config, st, proc_st, sourcepar, sourcepar_err)

    # Save output
    sourcepar_mean = write_output(config, sourcepar, sourcepar_err)

    # Save residuals
    spectral_residuals(config, spec_st, sourcepar_mean)

    # Plotting
    nspec = len(set(s.id[:-1] for s in spec_st))
    ncols = 4 if nspec > 6 else 3
    plot_spectra(config, spec_st, specnoise_st, plottype='regular',
                 ncols=ncols, async_plotter=plotter)
    plot_spectra(config, weight_st, plottype='weight',
                 ncols=ncols, async_plotter=plotter)
    if config.plot_station_map:
        # We import here, cause we check for Cartopy at runtime
        from sourcespec.ssp_plot_stations import plot_stations
        plot_stations(config, sourcepar)

    if config.html_report:
        html_report(config, sourcepar, sourcepar_err)

    ssp_exit()


if __name__ == '__main__':
    main()
