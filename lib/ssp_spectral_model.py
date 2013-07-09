# -*- coding: utf-8 -*-
# ssp_spectral_model.py
#
# Spectral model for source_spec
# (c) 2012 Claudio Satriano <satriano@ipgp.fr>
# (c) 2013 Claudio Satriano <satriano@ipgp.fr>,
#          Emanuela Matrullo <matrullo@geologie.ens.fr>
#          Agnes Chounet <chounet@ipgp.fr>
import math
import numpy as np

def spectral_model(freq, Mw, fc, t_star):
    # log S(w)= log(coeff*Mo) + log((1/(1+(w/wc)^2)) + log (exp (- w *t_star/2)) 
    # attenuation model: exp[-pi t* f] with t*=T /Q
    loge = math.log10(math.e)
    return Mw - np.log10(1. + np.power((freq/fc), 2)) - loge * (math.pi * freq * t_star)
    
def objective_func(xdata, ydata, weight):
###A### here ydata, xdata should be np.array ; and params_0 components should be np.float.
###A### and the 'func' argument should be the 'spectral_model' function, or a function 
###A### that takes exactly 4 parameters. (objective_fun is designed for this application only)
    errsum = np.sum(weight)
    def objective_func2(params):
        model = spectral_model(xdata, params[0], params[1], params[2])
        res = np.array(ydata) - np.array(model)
        res2 = np.power(res, 2)
        wres = np.array(weight) * np.array(res2)
        return np.sqrt(np.sum(wres)/errsum)
    return objective_func2
    
def callback(x):
    pass
    #print 'parameters:', x
