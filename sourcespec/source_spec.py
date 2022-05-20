# -*- coding: utf-8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Earthquake source parameters from inversion of S-wave spectra.

:copyright:
    2012 Claudio Satriano <satriano@ipgp.fr>

    2013-2014 Claudio Satriano <satriano@ipgp.fr>,
              Emanuela Matrullo <matrullo@geologie.ens.fr>,
              Agnes Chounet <chounet@ipgp.fr>

    2015-2022 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)


def main():
    # Lazy-import modules for speed
    from sourcespec.ssp_parse_arguments import parse_args
    options = parse_args(progname='source_spec')

    # Setup stage
    from sourcespec.ssp_setup import (
        configure, move_outdir, remove_old_outdir, setup_logging,
        save_config, init_plotting, ssp_exit)
    config = configure(options, progname='source_spec')
    setup_logging(config)

    from sourcespec.ssp_read_traces import read_traces
    st = read_traces(config)

    # Now that we have an evid, we can rename the outdir and the log file
    move_outdir(config)
    setup_logging(config, config.hypo.evid)
    remove_old_outdir(config)

    # Save config to out dir
    save_config(config)

    # Deconvolve, filter, cut traces:
    from sourcespec.ssp_process_traces import process_traces
    proc_st = process_traces(config, st)

    # Build spectra (amplitude in magnitude units)
    from sourcespec.ssp_build_spectra import build_spectra
    spec_st, specnoise_st, weight_st = build_spectra(config, proc_st)

    plotter = init_plotting(config)
    ntr = len(set(t.id[:-1] for t in proc_st))
    ncols = 4 if ntr > 6 else 3
    from sourcespec.ssp_plot_traces import plot_traces
    plot_traces(config, proc_st, ncols=ncols, async_plotter=plotter)

    # Spectral inversion
    from sourcespec.ssp_inversion import spectral_inversion
    sourcepar = spectral_inversion(config, spec_st, weight_st)

    # Local magnitude
    from sourcespec.ssp_local_magnitude import local_magnitude
    if config.compute_local_magnitude:
        local_magnitude(config, st, proc_st, sourcepar)

    # Compute averages, find outliers
    from sourcespec.ssp_averages import compute_averages
    compute_averages(config, sourcepar)

    # Save output
    from sourcespec.ssp_output import write_output
    write_output(config, sourcepar)

    # Save residuals
    from sourcespec.ssp_residuals import spectral_residuals
    spectral_residuals(config, spec_st, sourcepar)

    # Plotting
    from sourcespec.ssp_plot_spectra import plot_spectra
    nspec = len(set(s.id[:-1] for s in spec_st))
    ncols = 4 if nspec > 6 else 3
    plot_spectra(config, spec_st, specnoise_st, plot_type='regular',
                 ncols=ncols, async_plotter=plotter)
    plot_spectra(config, weight_st, plot_type='weight',
                 ncols=ncols, async_plotter=plotter)
    if config.plot_station_map:
        # We import here, cause we check for Cartopy at runtime
        from sourcespec.ssp_plot_stations import plot_stations
        plot_stations(config, sourcepar)

    if config.html_report:
        from sourcespec.ssp_html_report import html_report
        html_report(config, sourcepar)

    ssp_exit()


if __name__ == '__main__':
    main()
