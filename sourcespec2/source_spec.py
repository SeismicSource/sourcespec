# -*- coding: utf-8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Earthquake source parameters from inversion of P- or S-wave spectra.

:copyright:
    2012 Claudio Satriano <satriano@ipgp.fr>

    2013-2014 Claudio Satriano <satriano@ipgp.fr>,
              Emanuela Matrullo <matrullo@geologie.ens.fr>,
              Agnes Chounet <chounet@ipgp.fr>

    2015-2026 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""


def main():
    """Main routine for source_spec."""
    # pylint: disable=import-outside-toplevel
    # Lazy-import modules for speed
    from .ssp_parse_arguments import parse_args
    options = parse_args(progname='source_spec')

    # Setup stage
    from .config import config, configure_cli
    configure_cli(options, progname='source_spec')
    from .ssp_setup import (
        move_outdir, remove_old_outdir, setup_logging,
        save_config, ssp_exit)
    setup_logging()

    from .ssp_read_traces import read_traces
    st = read_traces()

    # Now that we have an evid, we can rename the outdir and the log file
    move_outdir()
    setup_logging(config.event.event_id)
    remove_old_outdir()

    # Save config to out dir
    save_config()

    # Deconvolve, filter, cut traces:
    from .ssp_process_traces import process_traces
    proc_st = process_traces(st)

    # Build spectra (amplitude in magnitude units)
    from .ssp_build_spectra import build_spectra
    spec_st, specnoise_st, weight_st = build_spectra(proc_st)

    from .ssp_plot_traces import plot_traces
    plot_traces(st, suffix='raw')
    plot_traces(proc_st)

    # Spectral inversion
    from .ssp_inversion import spectral_inversion
    sspec_output = spectral_inversion(spec_st, weight_st)

    # Radiated energy and apparent stress
    from .ssp_radiated_energy import radiated_energy_and_apparent_stress
    radiated_energy_and_apparent_stress(spec_st, specnoise_st, sspec_output)

    # Local magnitude
    if config.compute_local_magnitude:
        from .ssp_local_magnitude import local_magnitude
        local_magnitude(st, proc_st, sspec_output)

    # Compute summary statistics from station spectral parameters
    from .ssp_summary_statistics import compute_summary_statistics
    compute_summary_statistics(sspec_output, spec_st, weight_st)

    # Save output
    from .ssp_output import write_output, save_spectra
    write_output(sspec_output)
    save_spectra(spec_st)

    # Save residuals
    from .ssp_residuals import spectral_residuals
    spectral_residuals(spec_st, sspec_output)

    # Plotting
    from .ssp_plot_spectra import plot_spectra
    plot_spectra(spec_st, specnoise_st, plot_type='regular')
    plot_spectra(weight_st, plot_type='weight')
    from .ssp_plot_stacked_spectra import plot_stacked_spectra
    plot_stacked_spectra(spec_st, weight_st, sspec_output)
    from .ssp_plot_params_stats import box_plots
    box_plots(sspec_output)
    if config.plot_station_map:
        from .ssp_plot_stations import plot_stations
        plot_stations(sspec_output)

    if config.html_report:
        from .ssp_html_report import html_report
        html_report(config, sspec_output)

    ssp_exit()


if __name__ == '__main__':
    main()
