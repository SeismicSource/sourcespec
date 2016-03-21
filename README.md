# source_spec
**Earthquake source parameters from S-wave displacement spectra**

(c) 2011-2012 Claudio Satriano <satriano@ipgp.fr>

(c) 2013-2014 Claudio Satriano <satriano@ipgp.fr>,
              Emanuela Matrullo <matrullo@geologie.ens.fr>,
              Agnès Chounet <chounet@ipgp.fr>

(c) 2015-2016 Claudio Satriano <satriano@ipgp.fr>


## Description
The source_spec package is made of three codes:

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


## Sample runs
Several sample runs are available in the `source_spec_testruns` package.
