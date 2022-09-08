# -*- coding: utf-8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Earthquake source parameters from inversion of P- or S-wave spectra.

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


def main():
    """Main routine for source_spec."""
    # Lazy-import modules for speed
    from sourcespec.ssp_parse_arguments import parse_args
    options = parse_args(progname='source_spec')

    # Setup stage
    from sourcespec.ssp_setup import (
        configure, move_outdir, remove_old_outdir, setup_logging,
        save_config, ssp_exit)
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

    from sourcespec.ssp_plot_traces import plot_traces
    plot_traces(config, proc_st)

    # Spectral inversion
    from sourcespec.ssp_inversion import spectral_inversion
    sspec_output = spectral_inversion(config, spec_st, weight_st)

    # Radiated energy
    from sourcespec.ssp_radiated_energy import radiated_energy
    radiated_energy(config, spec_st, specnoise_st, sspec_output)

    # Local magnitude
    if config.compute_local_magnitude:
        from sourcespec.ssp_local_magnitude import local_magnitude
        local_magnitude(config, st, proc_st, sspec_output)

    # Compute summary statistics from station spectral parameters
    from sourcespec.ssp_summary_statistics import compute_summary_statistics
    compute_summary_statistics(config, sspec_output)

    # Save output
    from sourcespec.ssp_output import write_output
    write_output(config, sspec_output)

    # Save residuals
    from sourcespec.ssp_residuals import spectral_residuals
    spectral_residuals(config, spec_st, sspec_output)

    # Plotting
    from sourcespec.ssp_plot_spectra import plot_spectra
    plot_spectra(config, spec_st, specnoise_st, plot_type='regular')
    plot_spectra(config, weight_st, plot_type='weight')
    from sourcespec.ssp_plot_params_stats import box_plots
    box_plots(config, sspec_output)
    if config.plot_station_map:
        from sourcespec.ssp_plot_stations import plot_stations
        plot_stations(config, sspec_output)

    if config.html_report:
        from sourcespec.ssp_html_report import html_report
        html_report(config, sspec_output)

    ssp_exit()


if __name__ == '__main__':
    main()
