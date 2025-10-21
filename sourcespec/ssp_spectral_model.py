# -*- coding: utf-8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Spectral model and objective function.

:copyright:
    2012 Claudio Satriano <satriano@ipgp.fr>

    2013-2014 Claudio Satriano <satriano@ipgp.fr>,
              Emanuela Matrullo <matrullo@geologie.ens.fr>,
              Agnes Chounet <chounet@ipgp.fr>

    2015-2025 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
import math
import numpy as np


def spectral_model(freq, Mw, fc, t_star, alpha=1.):
    r"""
    Spectral model for seismic source spectrum.

    This function computes the theoretical spectral model based on the Brune
    source model with attenuation, expressed in terms of moment magnitude.

    .. math::

        Y_{data} = M_w + \frac{2}{3} \left[ - \log_{10} \left(
                              1+\left(\frac{f}{f_c}\right)^2 \right) -
                              \pi \, f^\alpha t^*(f) \log_{10} e \right]

    See :ref:`Theoretical Background <theoretical_background>`

    Parameters
    ----------
    freq : float or array-like
         Frequency or array of frequencies (in Hz) at which to compute the
         spectral model.
    Mw : float
         Moment magnitude of the seismic event.
    fc : float
         Corner frequency (in Hz) of the source spectrum.
    t_star : float or callable
         Attenuation parameter t* (in seconds), equal to travel time divided
         by quality factor (tt/Q). Can be a constant value or a function of
         frequency.
    alpha : float, optional
         Frequency exponent for the attenuation term. Default is 1.0 for
         standard attenuation model.

    Returns
    -------
    float or array-like
         Logarithmic spectral amplitude(s) at the specified
         frequency/frequencies.
    """
    # log S(w)= log(coeff*Mo) + log((1/(1+(w/wc)^2)) + \
    #           log (exp (- w *t_star/2))
    # attenuation model: exp[-pi t* f] with t*=T /Q
    loge = math.log10(math.e)
    # Check if t_star is callable (function) or a constant
    t_star_val = t_star(freq) if callable(t_star) else t_star
    return (
        Mw -
        (2. / 3.) * np.log10(1. + np.power((freq / fc), 2)) -
        (2. / 3.) * loge * (math.pi * np.power(freq, alpha) * t_star_val)
    )


def objective_func(xdata, ydata, weight):
    """
    Objective function generator for bounded inversion.

    Parameters
    ----------
    xdata : array-like
        x data (frequencies)
    ydata : array-like
        y data (log spectral amplitudes)
    weight : array-like
        weights for each data point

    Returns
    -------
    callable
        objective function
    """
    errsum = np.sum(weight)

    def _objective_func(params):
        """
        Objective function for bounded inversion.

        Parameters
        ----------
        params : array-like
            Model parameters (Mw, fc, t_star, [alpha])

        Returns
        -------
        float
            RMS misfit
        """
        # params components should be np.float
        if len(params) == 4:
            model = spectral_model(xdata, params[0], params[1],
                                   params[2], params[3])
        else:
            model = spectral_model(xdata, params[0], params[1], params[2])
        res = np.array(ydata) - np.array(model)
        res2 = np.power(res, 2)
        wres = np.array(weight) * np.array(res2)
        return np.sqrt(np.sum(wres) / errsum)
    return _objective_func


def callback(_):
    """Empty callback function for bounded inversion."""
