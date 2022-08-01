<img src="logo/SourceSpec_logo.svg" width="600">

# SourceSpec

Earthquake source parameters from P- or S-wave displacement spectra

[![PyPI-badge]][PyPI-link]
[![license-badge]][license-link]
[![docs-badge]][docs-link]
[![DOI-badge]][DOI-link]

(c) 2011-2022 Claudio Satriano <satriano@ipgp.fr>

## Description

SourceSpec is a collection of command line tools to compute earthquake source
parameteres (seismic moment, corner frequency, radiated energy, source size,
stress drop) from the inversion of P-wave and S-wave displacement spectra
recorded at one or more seismic stations.
SourceSpec also computes attenuation parameters (t-star, quality factor) and,
as a bonus, local magnitude.

See [Madariaga (2010)][Madariaga2010] for a primer on earthquake source
parameters and scaling laws.

Go to section [Theoretical background](#theoretical-background) below to get
more information on how the code works. More details are available on the
official SourceSpec [documentation].

SourceSpec is written in Python and requires a working Python environment to
run (see [Installation](#installation) below). However, since SourceSpec is
based on command line, you don't have to know how to code in Python to use it.

The SourceSpec package is made of three command line tools:

- `source_spec`: Compute earthquake source parameters from the inversion
  of P- or S-wave spectra.
- `source_model`: Direct modelling of P- or S-wave spectra, based on
  user-defined earthquake source parameters.
- `source_residuals`: Compute station residuals from `source_spec` output.

## Getting started

### For the impatient

If you have seismic recordings in [miniSEED] format (e.g., `traces.mseed`),
metadata in [StationXML] format (e.g., `station.xml`) and event information in
[QuakeML] format (e.g., `event.xml`), then:

1. Generate a config file via `source_spec -S`;
2. Edit the config file variable `station_metadata` to point to `station.xml` file;
3. Run `source_spec -t traces.mseed -q event.xml`.

### Command line arguments

After successfully installed SourceSpec (see [Installation](#installation)
below), you can get help on the command line arguments used by each code by
typing from your terminal:

    source_spec -h

(or `source_model -h`, or `source_residuals -h`).

`source_spec` and `source_model` require you to provide the path to seismic
traces via the `--trace_path` command line argument (see [Supported
file formats](#supported-file-formats) below).

Information on the seismic event can be stored in the trace header ([SAC]
format), or provided through a [QuakeML] file (`--qmlfile`) or a [HYPO71] or
[HYPOINVERSE-2000] file (`--hypocenter`). See
[Supported file formats](#supported-file-formats) below for more information on
the supported file formats.

### Configuration file

`source_spec` and `source_model` require a configuration file. The default file
name is `source_spec.conf`, other file names can be specified via the
`--configfile` command line argument.

You can generate a sample configuration file through:

    source_spec -S

Take your time to go through the genereated configuration file (named
`source_spec.conf`): the comments within the file will guide you on how to set
up the different parameters.

## Supported file formats

### Trace formats

SourceSpec can read all the
[trace formats supported by ObsPy][obspy_trace_formats].

Two very common choices are:

- [miniSEED]
- [SAC]

The SAC format can carry additional information in its header, like event
location and origin time, phase picks, instrument sensitivity.

### Event formats

SourceSpec can read event information (event ID, location, origin time) in the
following formats:

- [QuakeML]: SourceSpec will also read phase picks and
  focal mechanism, if available
- [HYPO71]
- [HYPOINVERSE-2000]: SourceSpec will also read phase picks, if available

Event information can also be stored in the SAC file headers (header fields:
`EVLA`, `EVLO`, `EVDP`, `O`, `KEVNM`).

### Phase pick formats

Phase picks for P and S waves can be read from one of the following formats:

- [QuakeML]
- [HYPO71]
- [HYPOINVERSE-2000]

Phase picks can also be stored in the SAC file headers (header fields: `A` and
`T0`).

### Station metadata formats

Station metadata (coordinates, instrumental response) can be provided in one of
the following formats:

- [StationXML]
- [Dataless SEED]
- [SEED RESP]
- [SAC polezero (PAZ)]

Note that SEED RESP and PAZ formats do not contain station coordinates, which
should therefore be in the trace header (traces in SAC format).

The station metadata file name or file directory is provided in the
configuration file through the parameter `station_metadata`.

Alternatively, instrument sensitivity can be provided in the SAC header or as a
constant in the configuration file. In both cases, use the configuration
parameter `sensitivity`.

### Output files

The SourceSpec main code, `source_spec` will produce the following output files
(`EVID` is replaced by the actual event ID):

- `EVID.ssp.out`: text file containing the estimated earthquake source
  parameters (per station and average)
- `EVID.ssp.log`: log file in text format (including the command line arguments,
  for [reproducibility])
- `EVID.ssp.conf`: the input config file (for [reproducibility])
- `EVID-residuals.pickle`: station residuals in
  [Python pickle format][pickle]
- `EVID.xml`: updated StationXML file with the results of the SourceSpec
  inversion (only if an input StationXML file is provided)

The following plots will be created, in png or pdf format:

- `EVID.traces.png[.pdf]`: trace plots
- `EVID.ssp.png[.pdf]`: spectral plots
- `EVID.sspweight.png[.pdf]`: spectral weight plots
- `EVID.boxplot.png[.pdf]`: [box plots][box_plot] for the earthquake source
  parameters retrieved at each station
- Misfit plots, when using "grid search" or "importance sampling" for the
  spectral inversion

As an option, station maps can be created (requires [Cartopy]):

- `EVID.map_mag.png[.pdf]`: station map with symbols colored by estimated
  moment magnitude
- `EVID.map_fc.png[.pdf]`: station map with symbols colored by estimated
  corner frequency

As an option, the retrieved source parameters (per station and average) can be
appended to a [SQLite] database, whose path is defined in the configuration
file.

Finally, always as an option, `source_spec` can generate a report in HTML
format.

## Theoretical background

The code computes spectra of the two horizontal components (and optionally of
the vertical component, as well), and then modulus as:

    sqrt(c1(w)^2+c2(w)^2)

It then inverts spectra for a 3-parameter model (Mw, Fc, t﹡ = T/Qs) using
initial values for Mw, fc and t﹡:

    log S(w)= log(coeff*Mo) + log((1/(1+(w/wc)^2)) + log(exp(-w*t_star/2))

It plots observed and inverted spectra on a single log-log graph (Mo vs
log-frequency).
Computes average and standard deviation of Mw, Mo, fc, t*, source radius and
Brune stress drop.

To get help:

    source_spec -h
    source_model -h
    source_residuals -h

## Installation

SourceSpec requires at least Python 3.6. All the required dependencies will be
downloaded and installed during the setup process.

### Using pip and PyPI (preferred method)

The latest release of SourceSpec is available on the
[Python Package Index](https://pypi.org/project/sourcespec/).

You can install it easily through `pip`:

    pip install sourcespec

### From SourceSpec GitHub releases

Download the latest release from the
[releases page](https://github.com/SeismicSource/sourcespec/releases),
in `zip` or `tar.gz` format, then:

    pip install sourcespec-X.Y.zip

or

    pip install sourcespec-X.Y.tar.gz

Where, `X.Y` is the version number (e.g., `1.2`).
You don't need to uncompress the release files yourself.

### From SourceSpec GitHub repository

If you need a recent feature that is not in the latest release (see the
`unreleased` section in [CHANGELOG](CHANGELOG.md)), you want to use the source
code from the
[SourceSpec GitHub repository](https://github.com/SeismicSource/sourcespec).

For that, clone the project:

    git clone https://github.com/SeismicSource/sourcespec.git

(avoid using the "Download ZIP" option from the green "Code" button, since
version number is lost), then install the code from within the `sourcespec`
main directory by running:

    pip install .

## Documentation

The offical SourceSpec documentation can be find at
[sourcespec.readthedocs.io][documentation].

## Sample runs

Several sample runs are available in the
[sourcespec_testruns](https://github.com/SeismicSource/sourcespec_testruns)
repository.

## How to cite

If you used SourceSpec for a scientific paper, please cite it as:

> Satriano, C. (2022). SourceSpec – Earthquake source parameters from
> P- or S-wave displacement spectra (X.Y).
> [doi: 10.5281/ZENODO.3688587]

Please replace `X.Y` with the SourceSpec version number you used.

## References

Madariaga, R. (2010). Earthquake Scaling Laws. In "Extreme Environmental
Events", pp. 364–383, [doi: 10.1007/978-1-4419-7695-6_22]. Available on
[ResearchGate][Madariaga2010].

<!-- Document links -->
[PyPI-badge]: http://img.shields.io/pypi/v/sourcespec.svg
[PyPI-link]: https://pypi.python.org/pypi/sourcespec
[license-badge]: https://img.shields.io/badge/license-CeCILL--2.1-green.svg
[license-link]: http://www.cecill.info/licences.en.html
[docs-badge]: https://readthedocs.org/projects/sourcespec/badge/?version=latest
[docs-link]: https://sourcespec.readthedocs.io/en/latest/?badge=latest
[DOI-badge]: https://zenodo.org/badge/DOI/10.5281/zenodo.3688587.svg
[DOI-link]: https://doi.org/10.5281/zenodo.3688587
[documentation]: https://sourcespec.readthedocs.io

[obspy_trace_formats]: https://docs.obspy.org/packages/autogen/obspy.core.stream.read.html
[miniSEED]: http://ds.iris.edu/ds/nodes/dmc/data/formats/miniseed/
[SAC]: https://ds.iris.edu/ds/support/faq/17/sac-file-format/
[QuakeML]: https://quake.ethz.ch/quakeml/
[HYPO71]: https://pubs.er.usgs.gov/publication/ofr72224
[HYPOINVERSE-2000]: https://pubs.er.usgs.gov/publication/ofr02171
[StationXML]: http://docs.fdsn.org/projects/stationxml/en/latest/
[Dataless SEED]: https://ds.iris.edu/ds/nodes/dmc/data/formats/dataless-seed/
[SEED resp]: https://ds.iris.edu/ds/nodes/dmc/data/formats/resp/
[SAC polezero (PAZ)]: https://www.jakewalter.net/sacresponse.html
[pickle]: https://docs.python.org/3/library/pickle.html
[box_plot]: https://en.wikipedia.org/wiki/Box_plot
[Cartopy]: https://scitools.org.uk/cartopy/docs/latest
[reproducibility]: https://en.wikipedia.org/wiki/Reproducibility
[SQLite]: https://www.sqlite.org

[doi: 10.5281/ZENODO.3688587]: https://doi.org/10.5281/ZENODO.3688587
[doi: 10.1007/978-1-4419-7695-6_22]: https://doi.org/10.1007/978-1-4419-7695-6_22
[Madariaga2010]: https://www.researchgate.net/publication/226065848_Earthquake_Scaling_Laws
