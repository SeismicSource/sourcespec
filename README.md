## source\_spec: a python code for modelling S-wave displacement spectra and inverting source parameters 
(c) 2011-2012 Claudio Satriano <satriano@ipgp.fr>  
(c) 2013-     Claudio Satriano <satriano@ipgp.fr>,
              Emanuela Matrullo <matrullo@geologie.ens.fr>  
Derived from sspec_v1.0.sh by Aldo Zollo and Claudio Satriano  

### Description:
source_spec computes the S-wave displacement spectra from stations recording a single event. 

It reads as input a file, a tgz archive or a directory (which can, in turn, contain
files and/or tgz archives) with traces in any format supported by ObsPy.

Optionally, one can specify:  
 - a file or a dir path containing station dataless  
 - a hypocenter file  
 - a phase file with P and S arrivals  

The code computes spectra of the two horizontal components and then modulus as:  

     sqrt(c1(w)^2+c2(w)^2)

It then inverts spectra for a 3-parameter model (Mw,Fc,t*=T/Qs) using initial
values for Mw, fc and t*:  

     log S(w)= log(coeff*Mo) + log((1/(1+(w/wc)^2)) + log (exp (- w *t_star/2)) 

It plots all spectra on a single log-log graph (Mw vs log-frequency).  
Plot obs vs theo spectra for each vectorial component.  
Computes average and st.dev of Mw, Mo, fc, source radius and Brune sd.  

### Sample runs:
To run the CRL test:

     ./source_spec.py testdata/CRL/2010.01.20-08.10.27 -c testconfig/config_CRL.py\
         -H testdata/CRL/2010.01.20-08.10.27.phs.h -p testdata/CRL/2010.01.20-08.10.27.phs 

To run the ISNet test:

     ./source_spec.py testdata/ISNet/14641r.full.sac.tgz -c testconfig/sample_config.py

To run the IPOC test:

     ./source_spec.py testdata/IPOC/324_0051-PB05-03077_tr14/ -c testconfig/config_IPOC.py

To get help:

     ./source_spec.py -h

To run source\_model test on CRL:

     ./source_model.py -c testconfig/config_CRL.py testdata/CRL/2010.01.20-08.10.27/2010.01.20-08.10.27.TRIZ.HHE.SAC\
          -H testdata/CRL/2010.01.20-08.10.27.phs.h --mag=2.6 --fc=8,12 --tstar=0.02,0.03   
