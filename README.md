<img src="logo/SourceSpec_logo.png" width="600">

# SourceSpec
**Earthquake source parameters from S-wave displacement spectra**

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3688587.svg)](
  https://doi.org/10.5281/zenodo.3688587)

(c) 2011-2012 Claudio Satriano <satriano@ipgp.fr>

(c) 2013-2014 Claudio Satriano <satriano@ipgp.fr>,
              Emanuela Matrullo <matrullo@geologie.ens.fr>,
              Agnès Chounet <chounet@ipgp.fr>

(c) 2015-2020 Claudio Satriano <satriano@ipgp.fr>


## Description
The SourceSpec package is made of three codes:

 - source_spec
 - source_model
 - source_residuals

The code computes spectra of the two horizontal components (and optionally of
the vertical component, as well), and then modulus as:

    sqrt(c1(w)^2+c2(w)^2)

It then inverts spectra for a 3-parameter model (Mw, Fc, t* = T/Qs) using
initial values for Mw, fc and t*:

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
Simply uncompress the archive and run the codes from `bin` directory.
Optionally, you can install the codes by running:

    pip install .

(use `pip install -e .` to install in developer mode), or:

    pip install sourcespec-x.x.tar.gz

(where `x.x` is the version number).


## Sample runs
Several sample runs are available in the
[sourcespec_testruns](https://github.com/claudiodsf/sourcespec_testruns)
package.
