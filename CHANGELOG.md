# SourceSpec Changelog

Earthquake source parameters from P- or S-wave displacement spectra

Copyright (c) 2011-2025 Claudio Satriano <satriano@ipgp.fr>

## unreleased

### Input/output

- Fix `np.float64` being printed in logs and YAML output when using NumPy>=2
- HTML report improvements:
  - Link to supplementary files and plots in the HTML report
  - Display configuration, log and output files in a modal window with syntax
    highlighting (only when the HTML report is served by a web server)

### Processing

- New option `r_power_n_segmented` for the `geom_spread_model` config parameter
  to use a segmented geometrical spreading model with different powers for
  different distance ranges
- New config parameter `refine_theoretical_arrivals` to refine the
  theoretical P and S arrival times using a simple autopicker based on the
  smoothed envelope of the trace
- New config parameter `clipping_min_amplitude_ratio` to set a threshold for
  trace amplitude below which the trace is not checked for clipping

### Plotting

- New plot: raw traces
- Stacked spectra: color spectral curves according to the weighting function
- Spectral plots: show information on the reason why a fit failed
- `plot_sourcepars`: possibility of selecting the latest runid for each event
- `plot_sourcepars`: new plot type: static stress drop vs. depth
- `plot_sourcepars`: option to color data points in scatter plot by another
  parameter
- Support for Matplotlib 3.9
- Maps: use PlateCarree projection for small maps (map diagonal below 100 km),
  which speeds up the plotting
- Maps: possibility to use a GeoTIFF file as a basemap
- Use threading to speed up saving plots to disk

### Config file

- New option `r_power_n_segmented` for the `geom_spread_model` config parameter
- New config parameters: `geom_spread_n_exponents`, `geom_spread_n_distances`
- New config parameters: `refine_theoretical_arrivals`, `autopick_freqmin`,
  `autopick_debug_plot`
- New config parameter `clipping_min_amplitude_ratio`
- New option `geotiff` for `plot_map_style` to use a GeoTIFF file as a basemap
- New config parameters: `plot_map_geotiff_filepath`,
  `plot_map_geotiff_grayscale`, `plot_map_geotiff_attribution`
- Improved documentation for the `sn_min` and `spectral_sn_min` parameters

### Bugfixes

- Fix for rejected spectra still being plotted in the stacked spectra plot
- Fix for corner case where all the inversion errors are zero
- Fix I/O error when reading PAZ files
- Fix for event ids in SourceSpec event file being only numbers: they are now
  correctly interpreted as strings
- Fix for ignored picks for certain kind of QuakeML files (like the ones from
  USGS) where the phase name is not in the `<phase_hint>` attribute of the
  `<pick>` element but in the `<phase>` attribute of the corresponding
  `<arrival>` element
- Stabilize the computation of derivatives in the Okal geometrical spreading
  model, thus avoiding numerical instabilities.
- Speed up map plotting when using coastlines
- Gently skip building spectra for station with very short data windows (less
  than 10 samples)

## v1.8 - 2024-04-07

This long overdue release brings many improvements, new features and bugfixes
gradually introduced during the last year.

Release highlights:

- New file formats for events and spectra
- New configuration options to better specifify velocity and density models
- Better support for P-wave inversion and teleseismic events
- Radiation pattern correction from focal mechanism in radiated energy
  computation
- Option for travel time-based signal window length
- More options to control the calculation of source radius and stress drop
- Improved estimation of radiated energy
- New source parameter: apparent stress
- Fix for map tiles not plotted anymore

This release requires at least Python 3.7.

Warning: the SQLite database used by this version is not compatible with
previous versions. You will need to upgrade your old database manually or using
`source_spec -u DATABASE_FILE_NAME`.

Make sure to read the detailed Changelog below 👇

### v1.8: Input/output

- Introducing a new file format for providing event information
  (hypocentral location, magnitude, focal mechanism, moment tensor):
  the [SourceSpec Event File].
- New HDF5 and TEXT file formats to store spectra
- Station residuals are now saved in an HDF5 spectrum file, instead of a
  pickle file
- New config file option `save_spectra` to save the spectra to an HDF5 file
  in the output directory
- Changes in the YAML output file:
  - `bsd` (Brune stress drop) parameter renamed to `ssd` (static stress drop)
  - Store in the `event_info` section the values of vp, vs and rho close to
    the hypocenter
  - Store in the `inversion_info` section the type of wave used for the
    inversion (P, S, SV or SH)
- Changes in the SQLite database (warning: these changes break compatibility
  with previous database versions):
  - `bsd` (Brune stress drop) parameter renamed to `ssd` (static stress drop)
  - Store the `Stations` table information on whether each parameter is an
    outlier (see [#38])
  - Make place in the `Stations` table for station-level errors on radiated
    energy (even if they are currently not computed)
  - Store in the `Events` table the number of observations used for computing
    each summary parameter
  - Store in the `Events` table weighted means for radiated energy and local
    magnitude, even if those means are currently the same as the simple means,
    since those parameters do not have station-level errors defined
  - New columns for apparent stress in both `Events` and `Stations` tables
  - Store in the `Events` table the values of vp, vs and rho close to the
    hypocenter
  - Store in the `Events` table the type of wave used for the inversion
    (P, S, SV or SH)
- New command line option (`-u` or `--updatedb`) to update an existing database
  from a previous version
- Input files are now linked symbolically in the `input_files` subdirectory
  of the output directory (not implemented for Windows)
- New command line option (`-R` or `--run_id_subdir`) to use `run_id`
  (if defined) as a subdirectory of the event directory
- Print event info to console and to log file
- HTML report improvements:
  - Event name in the summary table, if available
  - Author, agency and run completion date in the summary table
  - SourceSpec version in the inversion information table
  - Link to input files
  - Information on the type of wave used for the inversion (P, S, SV or SH)

### v1.8: Processing

- Use all the available components to compute P-wave spectra (previously,
  only the vertical component was used)
- Possibility of specifying a free surface amplification factor different
  from 2
- Possibility of specifying a layered velocity and density model for the
  source
- Possibility of specifying a different density for the source and for the
  stations
- If density is not provided (i.e., it is `None`), use the density from the
  global velocity model "iasp91"
- Teleseismic geometrical spreading model (Okal, 1992)
- New weighting option based on inverse frequency, so that lower frequencies
  have larger weight in the inversion. If traces contain noise, weights will
  be set to zero where SNR < 3 (see [#37])
- For weights computed from spectral S/N ratio (noise weighting), set to zero
  all the weights below 20% of the maximum weight, so that these weakly
  constrained parts of the spectrum are ignored in the inversion
- Possibility of using variable signal window lengths for each station
  as a function of the travel time of the P or S wave (see [#48])

### v1.8: Inversion

- Possibility of using the magnitude (or scalar moment) provided in the event
  file as initial Mw value for the inversion
- Reintroduced the possibility of providing the variability around the initial
  Mw value
- By combining the previous options, it is now possible to fix the Mw value
  during the inversion to the value provided in the event file

### v1.8: Post-Inversion

- Possibility of choosing the "k" coefficient to compute source radius from
  corner frequency (Kaneko and Shearer, 2014)
- Better control on the frequency range used for computing radiated energy
  (see [#49])
- Use station-specific radiation pattern (when available) for computing
  radiated energy
- Take into account for energy partition when computing radiated energy
  (Boatwright and Choy, 1986). This affects mostly the radiated energy
  computed from P waves
- New source parameter: apparent stress
- For parameters with no station-level uncertainty defined (currently,
  radiated energy and local magnitude), use simple mean when computing summary
  weighted averages (the previous behavior was to not compute weighted averages
  for these parameters)

### v1.8: Plotting

- Show the station radiated energy (Er) value on the station spectra plots
- Show the summary radiated energy (Er) value on the stacked spectra plot
- Station maps improvements:
  - Possibility of choosing a basemap style or no basemap
  - Possibility of not plotting the coastlines
  - Exclude outliers when computing colorbar limits
  - Improved computation of bounding box for regional and teleseismic events
  - Use a global orthographic projection when using stations at large
    teleseismic epicentral distances (more than 3000 km)
- Changes to `plot_sourcepars`:
  - Read vp, vs and rho from the SQLite database (previously: vs was hardcoded
    to 3.5 km/s, rho to 2700 kg/m3 and vp was not used)
  - Read the source radius "k" coefficient from the SQLite database
    (previously: "k" was hardcoded to 0.3724, value for the Brune model)
  - New command line option `--wave_type` to select the wave type (P, S, SV
    or SH) for plots involving the corner frequency
  - Possibility of plotting histogram of apparent stress
  - Option to filter events by apparent stress

### v1.8: Config file

- New config parameter `epi_dist_ranges` to select stations within one or
  more ranges of epicentral distances. It replaces the old parameter
  `max_epi_dist`.
- New config parameter `free_surface_amplification` to specify the free surface
  amplification factor (default: 2)
- New config parameter `layer_top_depths` to specify the depth of the top of
  the layers in a layered velocity and density model
- The config parameters `vp_source`, `vs_source` and `rho_source` can now be
  lists of values, to specify a layered velocity and density model for the
  source
- Config parameter `rho` renamed to `rho_source`
- New config parameter `rho_stations`
- New config parameter `geom_spread_min_teleseismic_distance` to set the
  minimum epicentral distance for using the teleseismic geometrical
  spreading model
- New config parameters `kp` and `ks` to set the "k" coefficient for
  computing source radius from corner frequency
- Config parameter `pi_bsd_min_max` renamed to `pi_ssd_min_max`
- New option `inv_frequency` for the config parameter `weighting` (see [#37])
- Config parameter `max_freq_Er` replaced by `Er_freq_range` (see [#49])
- New parameters, `qml_event_description` and `qml_event_description_regex`,
  to obtain the event name from the QuakeML event "description" tag
- New parameter `Mw_0_from_event_file` to use the magnitude (or scalar moment)
  provided in the event file as initial Mw value for the inversion
- Reintroduced the parameter `Mw_0_variability` to set the variability around
  the initial Mw value
- New parameter `plot_save_asap` to save plots as soon as they are ready.
  This uses less memory but slows down the code.
- New parameter `plot_map_style` to choose the map style
- New parameter `plot_map_api_key` to provide a Stadia Maps
  api key for Stamen Terrain basemap
- New option for the parameter `plot_coastline_resolution`: `no_coastline`
- New config parameter `variable_win_length_factor` to specify window
  length as a fraction of the travel time of the P/S wave (see [#48])

### v1.8: Bugfixes

- Fix source radius computation when using P waves (use P-wave velocity instead
  of S-wave velocity)
- Do not ignore picks labeled with lowercase "p" or "s"
- Fixed: config parameter `p_arrival_tolerance` was used also for S waves,
  instead of `s_arrival_tolerance` (see [#35])
- Fix Boatwright spreading model (log10 instead of natural log)
- Fix bug where signal and noise windows were plotted with the wrong length,
  under certain circumstances (see [#35])
- Fixes related to records with short signal windows (see [#39])
- Fix for beachball not plotted anymore with recent versions of Matplotlib.
- Fix bug where traces ignored because of low spectral S/N ratio, where still
  plotted as if they were valid traces
- Fix bug when specifying an absolute path for output directory: the path
  was treated as relative (see [#40])
- Fix bug where paths starting with tilde (~) were not parsed correctly
  (see [#43] and [#44])
- Fix bug where local magnitude was not written to the HYPO71 output file,
  when using weighted mean as reference statistics
- Fix for Stamen Terrain basemap now requiring an API key from Stadia Maps

### v1.8: Requirements

- Python minimum version raised to 3.7
- Matplotlib minimum version raised to 3.2
- Cartopy minimum version raised to 0.21

## v1.7 - 2023-03-31

This release improves trace processing through the use of modern routines for
instrument correction, optional baseline removal, a new clipping detection
algorithm and a better definition of signal and noise time windows.

A new plot is introduced, "stacked spectra" which allows to compare all the
spectra at once (and easily detect problematic stations 😉).
Also, a new command line tool, `plot_sourcepars`, allows making aggregate
plots of source parameters for many events (starting from the SQLite database).

New config file parameters have been added, while some have been removed.
Please run `source_spec -U CONFIG_FILE_NAME` to update your old config file.

As always, many bugfixes and improvements have been made in this release.
Thanks to all the users who took time to write and ask questions (by mail
or using the official SourceSpec [Discussions]).

Big kudos to Kris Vanneste [@krisvanneste] who helped all along the development
with code review and testing and submitted pull requests on noise windows and
clipping detection.

Below is the detailed Changelog 👇

### v1.7: Input/output

- Possibility of using a single PAZ file as a "generic" PAZ file for all the
  stations
- Command line option `--station_metadata` (or `-w`) for overriding the config
  file parameter with the same name (see pull request [#16])
- Removed command line option `--no-response` for avoiding removing instrument
  response (use the config option `correct_instrumental_response` instead)
- New output file in YAML format. The old `.out` file is still available but
  deprecated.
- Information on the inversion procedure in YAML and HTML output
- Option to add an agency logo to the HTML page
- Possibility of generating HTML report without figures (see [#30])

### v1.7: Processing

- Use modern ObsPy `trace.remove_response()` routine for instrument correction
  (see [#27])
- Option to remove the trace baseline after instrument correction and before
  filtering (`remove_baseline` config parameter) (see [#25])
- New algorithms for clipping detection based on kernel density estimation of
  the trace amplitude values (see [#23], [#24], [#25])
  - Two methods are available:
    - `clipping_score`: compute a trace clipping score based on the shape of
      the kernel density estimation.
    - `clipping_peaks`: check if trace is clipped, based on the number of peaks
      in the kernel density estimation;
  - Use `clipping_detection_algorithm` in the config file to choose the
    algorithm and the other `clipping_*` parameters to adjust the results.
  - The algorithms can also be called from the command line, e.g. for debug
    purposes, using the shell command `clipping_detection`.
- Relax noise window requirements if noise weighting is not used.
  This is useful for older triggered records with noise windows
  that are short or even missing entirely (see pull request [#18])
- Some small improvements were made in the window definitions
  (see pull request [#18]):
  - Generate error if signal window is incomplete (P- or S-arrival
    before the start time of the trace)
  - Warn if noise window overlaps with P-window (instead of S-window, as in
    previous versions)
  - Constrain `signal_pre_time` for S-phase to half the S-P interval, if this
    interval is shorter than `signal_pre_time` (i.e., for short-distance
    records with short S-P interval)
- Extract source and station P and S velocities from global 'iasp91' velocity
  model, if both `v(p,s)_source` and `v(p,s)_stations` are set to `None`
  (see [#20])
- Magnitude limits for inversion are now autoset between 90% of the minimum
  of the spectral plateau and 110% of its maximum (see [#22])

### v1.7: Post-Inversion

- Possibility of choosing the reference summary statistics that will be used
  for map plots, QuakeML and HYPO output, as well as for the "Event Summary"
  section in HTML report and for computing station spectral residuals.
  Available summary statistics are:
  - mean
  - weighted_mean
  - percentiles (new!)
- Possibility of defining the number of sigmas for uncertainties on event means
  and weighted means

### v1.7: Plotting

- New plot: stacked spectra
- Do not zero-pad traces to common length when plotting, so that missing data
  at beginning or at the end can be easily detected (see [#21])
- Plot noise and signal windows separately for each component (see [#21])
- Show on the trace plot the reason why a trace has been ignored
- Logscale for boxplots, if parameters span a large interval
  (see pull request [#15])
- Support for SVG format for plot files (can be used in HTML output as
  alternative to PNG)
- Improved trace plot quality for vector formats (PDF, SVG)
- New command line tool: `plot_sourcepars` to make 1D or 2D plot of source
  parameters from a sqlite parameter file.

### v1.7: Config file

- Removed `sensitivity_only` option from `correct_instrumental_response`
- Removed config parameter: `Mw_0_variability`
- Removed config parameter: `clip_max_percent`
- New config parameter: `remove_baseline`
- New config parameters for clipping detection:
  - `clipping_detection_algorithm`
  - `clipping_debug_plot`
  - `clipping_score_threshold`
  - `clipping_peaks_sensitivity`
  - `clipping_peaks_percentile`
- Config file section `AVERAGES PARAMETERS` renamed to
  `SUMMARY STATISTICS PARAMETERS`
- New config parameter: `reference_statistics`
- New config parameter: `n_sigma`
- New config parameters for percentiles calculation: `lower_percentage`,
  `mid_percentage` and `upper_percentage`
- New config parameters for filtering and spectral windowing of displacement
  signals:
  - `bp_freqmin_disp`, `bp_freqmax_disp`
  - `freq1_disp`, `freq2_disp`
- Default values for `t_star_min_max` (instead of `None`)

### v1.7: Code improvements

- Large refactoring of the whole codebase, to make the code more modern and
  easier to maintain (see [#28])

### v1.7: Bugfixes

- Properly ignore vertical components when `ignore_vertical` is `True`
- Fix a bug preventing reading phase picks from HYPOINVERSE-2000 files
- Fix for noise window not showing up in PNG trace plots in some cases
- Fix reading velocities from NLL model (see [#20])
- HTML report: better scrollbars for station table across all the browsers
- Fix for cropped map for very large station-to-event distances (greater
  than 500 km)
- Fix a bug in generating evid form origin time when reading origin time from
  SAC header and the number of seconds was 59
- Fix a crash when no map tiles were available at the selected zoom level
- Fix for a corner case where the three components of the same instrument
  have different trace length (see [#31])
- Fix `source_residuals`, which didn't work anymore

## v1.6 - 2022-08-02

This release introduces several modifications to the config file.
You will need to upgrade your old config files manually or using
`source_spec -U CONFIG_FILE_NAME`.

This release requires at least Python 3.6.

A lot of effort has been devoted to improve the documentation.
Please check it out on
[https://sourcespec.readthedocs.io/](https://sourcespec.readthedocs.io/)

### v1.6: Input/output

- QuakeML output (when using QuakeML input)
- Command line option `--run-id` to provide a string identifying the
  current run (see pull request [#6])
- Write SourceSpec version to parfile
- Write SourceSpec version and run complete time to SQLite file
- Write author and agency info (if specified) to output files, figures and
  HTML report
- HTML page for misfit plots (when using grid search or importance sampling)
- Station table in HTML report is now sortable (and its header remains fixed)!
- Reduce PNG figures file size, while improving their resolution 😃
- Removed option to read event information and traces from a pickle file
  (rarely used)

### v1.6: Processing

- Support for P-wave spectral inversion (see pull request [#9])
- It is now possible to provide different vp and vs velocities, close to the
  source and close to the stations (see the new config options above and
  issue [#5])
- Possibility to choose a geometrical spreading model between
  (see issue [#8]):
  - rⁿ (default: n=1 – body waves)
  - Boatwright et al. (2002): "r" below a cutoff distance, frequency-dependent
    above the cutoff distance
- Use travel time to compute quality factor from t* (and viceversa)
  (see issue [#5])
- Compute travel time from pick and origin time, when possible
  (see issue [#10])
- Warn if noise window ends after P or S arrival

### v1.6: Post-Inversion

- Subtract the integral of noise spectrum from the integral of signal spectrum
  when computing radiated energy, under the hypothesis that energy is additive
  and noise is stationary

### v1.6: Config file

- Config parameter `paz` has been removed and merged into `station_metadata`
- Config parameters `vp` and `vs` have been renamed to `vp_source`
  and `vs_source` (see issue [#5])
- New, optional, config parameter `vp_stations` and `vs_stations`
- Config parameter `pre_p_time` and `pre_s_time` have been renamed to
  `noise_pre_time` and `signal_pre_time`, respectively (see pull request [#9])
- Config parameter `rps_from_focal_mechanism` renamed to
  `rp_from_focal_mechanism` (see pull request [#9])
- New config parameter: `geom_spread_model` (see issue [#8])
- Config parameters `PLOT_SHOW`, `PLOT_SAVE` and `PLOT_SAVE_FORMAT` are now
  lowercase (`plot_show`, `plot_save` and `plot_save_format`)
- New, optional, general config parameters for specifying author and agency
  information. This information is written to output files and figures,
  if specified
- New config parameter, `event_url`, to link the event page from the HTML
  report
- Removed `DEBUG` config parameter
- Parameters from `GENERAL PARAMETERS` section reorganized into a new
  section called `TRACE AND METADATA PARAMETERS`
- Some parameters from `INVERSION PARAMETERS` moved into a new section
  called `SPECTRAL MODEL PARAMETERS`

### v1.6: Bugfixes

- Fix for not working `weighting` options: `frequency` and `no_weight`
- Fix for negative weights occasionally generated by interpolation
- Fix bug when event coordinates are written into sqlite as binary blobs

## v1.5 - 2022-05-22

This is a pretty big release coming after several months of work on a large
dataset of more than 5000 events in Mayotte.

Please read through the changelog to discover all the improvements and new
features.

You will need to update your old config files via
`source_spec -U CONFIG_FILE_NAME`

Note that v1.5 is no more compatible with Python 2!

### v1.5: Input/output

- Write output files into a subdirectory of OUTDIR, whose name is
  the event id
- Support for HYPOINVERSE-2000 output files
- Removed autodetection of hypo71 file paths (specific to CRL case)
- SQLite output: added radiated energy, weighted averages,
  errors on parameters, number of observations, hypocentral location
  and origin time
- Removed `-C` argument to apply station correction to spectra. Now spectra
  are automatically corrected if `residuals_filepath` is specified in the
  configuration file
- Save an additional event parameter to output files: average quality factor
- Save additional station parameters to output files:
  source radius, Brune stress drop, source radius, quality factor
- Mark outliers in `.out` file and in html report
- Colored console output for log messages! (Not supported on Windows)

### v1.5: Processing

- New parameter for setting the width of the spectral smoothing window
  in terms of frequency decades: `spectral_smooth_width_decades`
  (see issue [#2])
- Compute spectral weights after spectral correction (when a station
  residuals file is specified via `residuals_filepath`)
- Removed configuration parameter `trace_format`
- New configuration parameter `sensitivity` to provide a constant sensor
  sensitivity (flat response curve), which overrides any response curve
  provided in metadata
- New parameter for manually specifying trace units: `trace_units`
  (defaults to `auto`)
- New approach for trace clipping detection (requires just one configuration
  parameter, named `clip_max_percent`)
  - Check for trace clipping only in the processing window
  - Use histogram of samples to detect clipping
- Fix for wrong component used for 'SV' spectra (see issue [#3])

### v1.5: Inversion

- New config option: `Mw_0_variability`. Allowed variability around `Mw_0`
  during the main inversion. Previously hardcoded to 0.1
- New inversion methods for grid sampling:
  - grid search (very slow!)
  - importance sampling of the misfit grid using k-d tree
    (faster, but less accurate)
- Fix for Basin-hopping algorithm not running

### v1.5: Post-Inversion

- New set of post-inversion parameters to reject certain inversion results,
  per-station: `pi_fc_min_max`, `pi_t_star_min_max`, `pi_bsd_min_max`,
  `pi_misfit_max`
- Reject inversion results when inverted `fc` is within 10% of `fc_min` or
  `fc_max`
- Fix: use logarithmic error width for weighted logarithmic averages
  - previous way of computing weighted logarithmic averages was not correct!
- Option to reject outliers using the IQR (interquartile range) method:
  parameter `nIQR`
- Support for non symmetric error on station spectral parameters
- Compute additional, per-station parameters: source radius, Brune stress
  drop and quality factor
- Compute errors for all station parameters
- Compute weighted averages for all event parameters (except radiated energy)
- Compute spectral residuals using weighted average spectral parameters

### v1.5: Plotting

- Source parameter box plots to evaluate parameter dispersion across stations
  and visually detect outliers
- Misfit plot (2D and 1D) when using grid sampling
- `cartopy` removed as installation dependency, since it is not easily
  installable via `pip`
- Use GSHHS database to draw coastlines.
  - New config option: `plot_coastline_resolution`
- Correctly show circles on maps with diagonal smaller than 1 km
- Fix plotting map colorbar on Matplotlib 3.5
- Make average and errorbar lines more visible on map colorbar
- Fix for error on plotting fc map, when only one station is available
- Fix trace plot scaling for traces with larger signal outside the plot
  window
- Do not plot 'H' spectrum if there is only one instrument component
  (since it will coincide with the only component)
- Plot uncorrected spectrum when station correction is used

## v1.4 - 2021-10-13

- New config option `rps_from_focal_mechanism` to compute station-specific
  S-wave radiation pattern from focal mechanism, if a focal mechanism is
  available in the QuakeML file
- Plot the focal mechanism on maps, if it is available
- Change default inversion algorithm to TNC (truncated
  Newton algorithm)
- Config option `dataless` renamed to `station_metadata`
- Config option `traceids` renamed to `traceid_mapping_file`
- Config options `ignore_stations` and `use_stations` renamed to
  `ignore_traceids` and `use_traceids`, respectively
- Support for 2D NonLinLoc grids (via `nllgrid >= 1.4.1`)
- Possibility of using a generic `DEFAULT` NonLinLoc time grid
- Added `cartopy` as an installation dependency
- Fixed: `nllgrid` was always requested at runtime
- Fixed: gracefully handle the case when there is no internet connection and
  map tiles cannot be downloaded
- Fixed (Windows): suppress colored terminal output, which is not supported
- Fixed (Windows): it is now possible to relaunch the same run, without
  having to delete the output directory first
- Fixed (Windows): use same timezone names than Linux and macOS

## v1.3.1 - 2021-09-13

- Fix for HTML report not showing trace and spectral plots
- HTML report: add Corner Frequency in Event Summary

## v1.3 - 2021-08-20

- HTML reports
- Option to provide station-specific spectral windowing

## v1.2 - 2021-05-20

- Use `python-versioneer` to manage version numbers

## v1.1 - 2021-04-07

- Bug fixes:
  - Accept band code `C` for broadband seismometers sampled at >=250 Hz
  - Require `cartopy>=0.18` for compatibility with `matplotlib>=3.3`

## v1.0 - 2021-03-03

- Simplification of time window parameters:
  - an unique window length, `win_length`, is used for time-domain S/N ratio
    calculation and spectral analysis
  - the old parameters, `noise_win_length` and `s_win_length`, are no more
    supported and generate an error
- Reorganized config file and improved inline documentation within the file
  - removed unused option: `fc_0`
  - removed `pre_filt` option
  - post-processing check on `fc` bounds has been removed
- Option to specify non standard instrument codes (e.g., `L` for
  acceleration)
- Option to specify filter frequencies for local magnitude computation
- Plotting improvements:
  - Show S/N ratio on trace plots
  - Show spectral S/N ratio on spectrum plots
  - Option to show/hide ignored traces/spectra (low S/N)
- Bug fixes:
  - Pay attention to location code when building average spectra
  - Plotting: avoid overlapping traces with different location code
  - Plotting: avoid overlapping spectra with different location code

## v0.9 - 2020-04-24

- Support for QuakeML input and StationXML
- Support for Python 3.5
- Only compatible with ObsPy >= 1.1.0
- Project reorganization:
  - Project renamed to SourceSpec
  - `ssp_residuals` renamed to `source_residuals`
  - New installable package (e.g., via `pip`)
- Spectra are smoothed in log-freq (no more Konno-Ohmachi)
- Inversion is performed in a log-freq space
- Option to invert for `t_star_0` on the plateau level
- Traces are filtered before computing S/N ratio
- Trace clipping detection
- Traces are always plotted, even if no inversion is performed
- Use by default a global model for theoretical travel time calculation
- Possibility of using NonLinLoc travel time grids (requires `nllgrid`)
- New options for P and S arrival time tolerance
- New option for maximum epicentral distance for trace processing
- Possibility of using a NonLinLoc model grid for obtaining vs at the source
  and at the receiver (requires `nllgrid`)
- Use `log10` of weighting function, since the inversion is done in magnitude
  units
- Use `json` format for `traceid` correction file
- Save config file to output dir (only for `source_spec`)
- Save run completion time into output file
- Logarithmic average and logarithmic (asymmetric) error bars for Mo, fc and
  source radius
- Computation of radiated energy
- Station-specific filters
- New parameters: `gap_max`, `overlap_max`
- Add legend to spectral plots
- Add event information, code version and run completion time to plots
- Multifigure plotting for traces and spectra (for large number of stations)
- New option to plot a station map, color-coded by station magnitude
  (requires Cartopy)
- Refactoring of local magnitude computation code:
  - Wood-Anderson amplitude is computed from the whole trace converted to W-A
    and not converting only the min and max peaks (which is the default in
    ObsPy)
  - Trace windowing is computed on HF envelopes
  - New option for custom Ml coefficients
  - Remove unnecessary Ml options
- Code cleaning and optimization:
  - Switch from `optparse` (deprecated) to `argparse`
  - Code style fixes and refactoring
- New option to update a config file from previous versions
- BUGFIX: Fix a major bug in reading hypo pick file
- BUGFIX: Fix for pick file not read in `source_model`

## v0.8 - 2014-07-11

- Trace plot showing S and noise time windows
- Improved handling of paz files
- Per-station Mo on output file
- Code cleaning and optimization

## v0.7 - 2014-04-07

- Code reorganization:
  - inversion code split to its own functions
- Option to use bounded inversion
- Station residuals, through `ssp_residuals.py`
- `source_model` can now output tables of trade-off between parameters
- Fix in the way noise trace is computed and processed
- Documentation!

## v0.6 - 2013-06-05

- Signal to noise weighting
- Improved local magnitude computation
- New options:
  - time domain integration
  - vertical component
- Largely improved plotting function
- More data formats supported (Antilles, IPOC)
- `source_model`: a code for plotting theoretical spectra
- Code refactoring

## v0.5 - 2013-02-10

- Azimuth computation
- Construction of an overall database
- Local magnitude computation
- konnoOhmachiSmoothing

## v0.4 - 2012-04-10

- Logging infrastructure
- Code reorganization

## v0.3 - 2012-02-10

- Output is no more printed at screen, but to file
- The plots can be saved to a file as well
- We differentiate between short periods and broad bands

## v0.2 - 2012-02-06

Extended and generalized for the CRL application.

## v0.1 - 2012-01-17

Initial Python port.

[Discussions]: https://github.com/SeismicSource/sourcespec/discussions
[@krisvanneste]: https://github.com/krisvanneste
[SourceSpec Event File]: https://sourcespec.readthedocs.io/en/latest/source_spec_event_file.html

[#2]: https://github.com/SeismicSource/sourcespec/issues/2
[#3]: https://github.com/SeismicSource/sourcespec/issues/3
[#5]: https://github.com/SeismicSource/sourcespec/issues/5
[#6]: https://github.com/SeismicSource/sourcespec/issues/6
[#8]: https://github.com/SeismicSource/sourcespec/issues/8
[#9]: https://github.com/SeismicSource/sourcespec/issues/9
[#10]: https://github.com/SeismicSource/sourcespec/issues/10
[#15]: https://github.com/SeismicSource/sourcespec/issues/15
[#16]: https://github.com/SeismicSource/sourcespec/issues/16
[#18]: https://github.com/SeismicSource/sourcespec/issues/18
[#20]: https://github.com/SeismicSource/sourcespec/issues/20
[#21]: https://github.com/SeismicSource/sourcespec/issues/21
[#22]: https://github.com/SeismicSource/sourcespec/issues/22
[#23]: https://github.com/SeismicSource/sourcespec/issues/23
[#24]: https://github.com/SeismicSource/sourcespec/issues/24
[#25]: https://github.com/SeismicSource/sourcespec/issues/25
[#27]: https://github.com/SeismicSource/sourcespec/issues/27
[#28]: https://github.com/SeismicSource/sourcespec/issues/28
[#30]: https://github.com/SeismicSource/sourcespec/issues/30
[#31]: https://github.com/SeismicSource/sourcespec/issues/31
[#35]: https://github.com/SeismicSource/sourcespec/issues/35
[#37]: https://github.com/SeismicSource/sourcespec/issues/37
[#38]: https://github.com/SeismicSource/sourcespec/issues/38
[#39]: https://github.com/SeismicSource/sourcespec/issues/39
[#40]: https://github.com/SeismicSource/sourcespec/issues/40
[#43]: https://github.com/SeismicSource/sourcespec/issues/43
[#44]: https://github.com/SeismicSource/sourcespec/issues/44
[#48]: https://github.com/SeismicSource/sourcespec/issues/48
[#49]: https://github.com/SeismicSource/sourcespec/issues/49
