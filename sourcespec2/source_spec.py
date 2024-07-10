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


def ssp_run(st, inventory, ssp_event, picks, allow_exit=False):
    """
    Run source_spec as function with collected traces, station inventory,
    event and picks

    :param st: Traces to be processed
    :type st: :class:`obspy.core.stream.Stream`
    :param inventory: Station metadata
    :type inventory: :class:`obspy.core.inventory.Inventory`
    :param ssp_event: Event information
    :type ssp_event: :class:`sourcespec.ssp_event.SSPEvent`
    :param picks: List of picks
    :type picks: list of :class:`sourcespec.ssp_event.Pick`
    :param allow_exit: whether to allow hard exit (without returning)
    :type allow_exit: bool

    :return: (proc_st, spec_st, specnoise_st, weight_st, sspec_output)
    :rtype: tuple of
        :class:`obspy.core.stream.Stream`,
        :class:`sourcespec.spectrum.SpectrumStream`,
        :class:`sourcespec.spectrum.SpectrumStream`,
        :class:`sourcespec.spectrum.SpectrumStream`,
        :class:`sourcespec.ssp_data_types.SourceSpecOutput`
    """
    # pylint: disable=import-outside-toplevel
    # Lazy-import modules for speed
    import os
    from .config import config

    # Avoid hard exit when ssp_exit is called somewhere
    if not allow_exit:
        from . import ssp_setup
        ssp_setup.SSP_EXIT_CALLED = True

    # Create output folder if required, save config and setup logging
    from .ssp_setup import (get_outdir_path, save_config, setup_logging)
    if config.options.outdir:
        outdir = get_outdir_path()
        if not os.path.exists(outdir):
            os.makedirs(outdir)
            # Setup logging
            # (it is assumed this is already done if outdir exists)
            setup_logging(config.event.event_id)
        # Save config to out dir
        save_config()

    # Preprocessing
    from .ssp_read_traces import (augment_event, augment_traces,
                                  select_components)
    augment_event(ssp_event)
    st = select_components(st)
    augment_traces(st, inventory, ssp_event, picks)

    # Deconvolve, filter, cut traces:
    from .ssp_process_traces import process_traces
    proc_st = process_traces(st)

    # Build spectra (amplitude in magnitude units)
    from .ssp_build_spectra import build_spectra
    spec_st, specnoise_st, weight_st = build_spectra(proc_st)

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

    return (proc_st, spec_st, specnoise_st, weight_st, sspec_output)


def ssp_output(st, proc_st, spec_st, specnoise_st, weight_st, sspec_output):
    """
    Function writing all output to disk.

    If no output directory is specified, no output is written.

    :param st: Original traces
    :type st: :class:`obspy.core.stream.Stream`
    :param proc_st: Processed traces
    :type proc_st: :class:`obspy.core.stream.Stream`
    :param spec_st: Signal spectra
    :type spec_st: :class:`sourcespec.spectrum.SpectrumStream`
    :param specnoise_st: Noise spectra
    :type specnoise_st: :class:`sourcespec.spectrum.SpectrumStream`
    :param weight_st: Spectral weights
    :type weight_st: :class:`sourcespec.spectrum.SpectrumStream`
    :param sspec_output: SourceSpec output
    :type sspec_output: :class:`sourcespec.ssp_data_types.SourceSpecOutput`
    """
    # pylint: disable=import-outside-toplevel
    # Lazy-import modules for speed
    from .config import config

    # TODO: check (again) if outdir exists
    if not config.options.outdir:
        return

    # Save output
    from .ssp_output import write_output, save_spectra
    write_output(sspec_output)
    save_spectra(spec_st)

    # Save residuals
    from .ssp_residuals import spectral_residuals
    spectral_residuals(spec_st, sspec_output)

    # Plotting
    from .ssp_plot_traces import plot_traces
    plot_traces(st, suffix='raw')
    plot_traces(proc_st)
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
        html_report(sspec_output)


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
        ssp_exit)
    setup_logging()

    # Read all required information from disk
    from .ssp_read_traces import (read_traces, read_station_inventory,
                                  read_event_and_picks)
    st = read_traces()
    trace1 = st[0] if len(st) else None
    st.sort()
    inventory = read_station_inventory()
    ssp_event, picks = read_event_and_picks(trace1)

    # Now that we have an evid, we can rename the outdir and the log file
    move_outdir()
    setup_logging(config.event.event_id)
    remove_old_outdir()

    # Run sourcespec function
    result = ssp_run(st, inventory, ssp_event, picks, allow_exit=True)
    (proc_st, spec_st, specnoise_st, weight_st, sspec_output) = result

    # Generate output
    ssp_output(st, proc_st, spec_st, specnoise_st, weight_st, sspec_output)

    ssp_exit()


if __name__ == '__main__':
    main()
