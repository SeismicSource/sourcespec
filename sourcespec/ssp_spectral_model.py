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


def objective_func(freq, ampl, weight, t_star=None):
    """
    Objective function generator for bounded inversion.

    Parameters
    ----------
    freq : array-like
        Frequencies (in Hz)
    ampl : array-like
        Log spectral amplitudes
    weight : array-like
        Weights for each data point
    t_star : float or callable, optional
        Set a fixed t_star value or function of frequencies;
        if None, t_star is a parameter to invert for

    Returns
    -------
    callable
        Objective function
    """
    errsum = np.sum(weight)

    def _objective_func(params):
        """
        Objective function for bounded inversion.

        Parameters
        ----------
        params : array-like
            Model parameters. The layout depends on whether t_star is fixed:

            - If t_star is fixed (passed to outer function):
              [Mw, fc] or [Mw, fc, alpha]
            - If t_star is free (None in outer function):
              [Mw, fc, t_star] or [Mw, fc, t_star, alpha]

        Returns
        -------
        float
            Root mean square (RMS) misfit between observed and modeled
            spectral amplitudes.
        """
        Mw, fc = params[0], params[1]
        # Determine parameter layout based on whether t_star is fixed
        if t_star is None:
            # params: [Mw, fc, t_star, (optional) alpha]
            t_star_param = params[2]
            alpha = params[3] if len(params) == 4 else None
        else:
            # params: [Mw, fc, (optional) alpha]
            t_star_param = t_star
            alpha = params[2] if len(params) == 3 else None
        kwargs = {
            'freq': freq,
            'Mw': Mw,
            'fc': fc,
            't_star': t_star_param
        }
        if alpha is not None:
            kwargs['alpha'] = alpha
        model = spectral_model(**kwargs)
        res = np.array(ampl) - np.array(model)
        res2 = np.power(res, 2)
        wres = np.array(weight) * np.array(res2)
        return np.sqrt(np.sum(wres) / errsum)
    return _objective_func


def callback(_):
    """Empty callback function for bounded inversion."""
