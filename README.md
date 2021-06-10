<img src="logo/SourceSpec_logo.png" width="600">

# SourceSpec

Earthquake source parameters from S-wave displacement spectra

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3688587.svg)](
  https://doi.org/10.5281/zenodo.3688587)

(c) 2011-2012 Claudio Satriano <satriano@ipgp.fr>

(c) 2013-2014 Claudio Satriano <satriano@ipgp.fr>,
              Emanuela Matrullo <matrullo@geologie.ens.fr>,
              Agnès Chounet <chounet@ipgp.fr>

(c) 2015-2021 Claudio Satriano <satriano@ipgp.fr>

## Description

The SourceSpec package is made of three codes:

- source_spec
- source_model
- source_residuals

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
`unreleased` section in [ChangeLog](ChangeLog.txt)), you want to use the source
code from the SourceSpec GitHub repository.

For that, clone the project:

    git clone https://github.com/SeismicSource/sourcespec.git

(avoid using the "Download ZIP" option from the green "Code" button, since
version number is lost), then install the code from within the `sourcespec`
main directory by running:

    pip install .

## Documentation

A very incomplete documentation can be found at
[sourcespec.readthedocs.io](https://sourcespec.readthedocs.io).

## Sample runs

Several sample runs are available in the
[sourcespec_testruns](https://github.com/SeismicSource/sourcespec_testruns)
repository.
