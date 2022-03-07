# SourceSpec
Earthquake source parameters from S-wave displacement spectra

(c) 2011-2022 Claudio Satriano <satriano@ipgp.fr>


## unreleased

- Input/output:
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
  - Colored console output for log messages! (Not supported on Windows)
- Processing:
  - New parameter for setting the width of the spectral smoothing window
    in terms of frequency decades: `spectral_smooth_width_decades`
  - Compute spectral weights after spectral correction (when a station
    residuals file is specified via `residuals_filepath`)
  - New approach for trace clipping detection (requires just one configuration
    parameter, named `clip_max_percent`)
    - Check for trace clipping only in the processing window
    - Use histogram of samples to detect clipping
  - Fix for wrong component used for 'SV' spectra (#3)
- Inversion:
  - New config option: `Mw_0_variability`. Allowed variability around `Mw_0`
    during the main inversion. Previously hardcoded to 0.1
  - New inversion methods for grid sampling:
    - grid search (very slow!)
    - importance sampling of the misfit grid using k-d tree
      (faster, but less accurate)
  - Fix for Basin-hopping algorithm not running
- Post-Inversion:
  - New set of post-inversion parameters to reject certain inversion results,
    per-station: `pi_fc_min_max`, `pi_t_star_min_max`, `pi_bsd_min_max`,
    `pi_misfit_max`
  - Reject inversion results when inverted `fc` is within 10% of `fc_min` or
    `fc_max`
  - Fix: use logarithmic error width for weighted logarithmic averages
    - previous way of computing weighted logarithmic averages was not correct!
  - Remove outliers (>2sigma) when computing event averages
  - Support for non symmetric error on station spectral parameters
  - Compute additional, per-station parameters: source radius, Brune stress
    drop and quality factor
  - Compute errors for all station parameters
  - Compute weighted averages for all event parameters (except radiated energy)
- Plotting:
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
