## source\_spec: a python code for modelling S-wave displacement spectra and inverting source parameters
(c) 2011-2012 Claudio Satriano <satriano@ipgp.fr>
(c) 2013-2014 Claudio Satriano <satriano@ipgp.fr>,
              Emanuela Matrullo <matrullo@geologie.ens.fr>,
              Agnès Chounet <chounet@ipgp.fr>

Derived from sspec_v1.0.sh by Aldo Zollo and Claudio Satriano

### Description:
source_spec computes the S-wave displacement spectra from stations recording a single event.

It reads as input a file, a tgz archive or a directory (which can, in turn, contain
files and/or tgz archives) with traces in any format supported by ObsPy.

Optionally, one can specify:
 - a file or a dir path containing station dataless
 - a hypocenter file
 - a phase file with P and S arrivals

The code computes spectra of the two horizontal components (and optionally of the vertical
component, as well), and then modulus as:

     sqrt(c1(w)^2+c2(w)^2)

It then inverts spectra for a 3-parameter model (Mw,Fc,t*=T/Qs) using initial
values for Mw, fc and t*:

     log S(w)= log(coeff*Mo) + log((1/(1+(w/wc)^2)) + log (exp (- w *t_star/2))

It plots observed and inverted spectra on a single log-log graph (Mo vs log-frequency).
Computes average and standard deviation of Mw, Mo, fc, t*, source radius and Brune stress drop.

### Sample runs:
To run the CRL test:

     ./source_spec.py testdata/CRL/2010.01.20-08.10.27 -c testconfig/config_CRL.conf\
         -H testdata/CRL/2010.01.20-08.10.27.phs.h -p testdata/CRL/2010.01.20-08.10.27.phs

To run the ISNet test:

     ./source_spec.py testdata/ISNet/14641r.full.sac.tgz -c testconfig/config_ISNet.conf

To run the IPOC test:

     ./source_spec.py testdata/IPOC/324_0051-PB05-03077_tr14/ -c testconfig/config_IPOC.conf

To get help:

     ./source_spec.py -h




### source\_model:
source\_model plots theoretical spectra (and optionally observed ones), given one or more
values for Mw, fc and t*.

To run source\_model test on CRL:

     ./source_model.py -c testconfig/config_CRL.conf testdata/CRL/2010.01.20-08.10.27/2010.01.20-08.10.27.TRIZ.HHE.SAC\
          -H testdata/CRL/2010.01.20-08.10.27.phs.h --mag=3.0 --fc=6,10 --tstar=0.05,0.07 -p


### ssp\_residuals:
ssp\_residuals computes station residuals from the output of source\_spec.
It takes multiple pickle files in the form:

     EVID-residual.pickle

containing single-event station residual, computes average station residuals, and store them into
a pickle file called:

     residual_mean.pickle
