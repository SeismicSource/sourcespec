# -*- coding: utf-8 -*-
# ssp_spectral_model.py
#
# (c) 2012 Claudio Satriano <satriano@ipgp.fr>
# (c) 2013 Claudio Satriano <satriano@ipgp.fr>,
#          Emanuela Matrullo <matrullo@geologie.ens.fr>,
#          Agnes Chounet <chounet@ipgp.fr>
'''
Spectral model and objective function.
'''
from __future__ import division
import math
import numpy as np

def spectral_model(freq, Mw, fc, t_star, alpha=1.):
    '''
    .. math::

       Y_{data} = M_w + \\frac{2}{3} \\left[ - \\log_{10} \\left(
                        1+\\left(\\frac{f}{f_c}\\right)^2 \\right) -
                        \\pi \\, f t^* \\log_{10} e \\right]

    see :ref:`source_spec <source_spec>` for a detailed derivation of this model.
    '''

    # log S(w)= log(coeff*Mo) + log((1/(1+(w/wc)^2)) + log (exp (- w *t_star/2))
    # attenuation model: exp[-pi t* f] with t*=T /Q
    loge = math.log10(math.e)
    return Mw -\
           (2./3.)*np.log10(1. + np.power((freq/fc), 2)) -\
           (2./3.)*loge * (math.pi * np.power(freq, alpha) * t_star)

def objective_func(xdata, ydata, weight):
    '''
    Objective function generator for bounded inversion.
    '''
    errsum = np.sum(weight)
    def _objective_func(params):
        #params components should be np.float
        if len(params)==4:
            model = spectral_model(xdata, params[0], params[1], params[2], params[3])
        else:
            model = spectral_model(xdata, params[0], params[1], params[2])
        res = np.array(ydata) - np.array(model)
        res2 = np.power(res, 2)
        wres = np.array(weight) * np.array(res2)
        return np.sqrt(np.sum(wres)/errsum)
    return _objective_func

def callback(x):
    pass
    #print 'parameters:', x
