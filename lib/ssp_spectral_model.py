# -*- coding: utf-8 -*-
# ssp_spectral_model.py
#
# Spectral model for source_spec
# (c) 2012 Claudio Satriano <satriano@ipgp.fr>
# (c) 2013 Claudio Satriano <satriano@ipgp.fr>,
#          Emanuela Matrullo <matrullo@geologie.ens.fr>
import math
import numpy as np

def spectral_model(freq, Mw, fc, t_star):
    # log S(w)= log(coeff*Mo) + log((1/(1+(w/wc)^2)) + log (exp (- w *t_star/2)) 
    # attenuation model: exp[-pi t* f] with t*=T /Q
    loge = math.log10(math.e)
    return Mw - np.log10(1. + np.power((freq/fc),2) ) - loge*(math.pi*t_star*freq) 
