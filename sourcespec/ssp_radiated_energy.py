# -*- coding: utf8 -*-
"""
Compute radiated energy from spectral integration.

:copyright:
    2012 Claudio Satriano <satriano@ipgp.fr>

    2013-2014 Claudio Satriano <satriano@ipgp.fr>,
              Emanuela Matrullo <matrullo@geologie.ens.fr>,
              Agnes Chounet <chounet@ipgp.fr>

    2015-2021 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement, Version 2.1
    (http://www.cecill.info/index.en.html)
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np


def radiated_energy(config, spec_st, sourcepar):
    """Compute radiated energy, using eq. (3) in Lancieri et al. (2012)."""
    # Loop on '??H' spectra
    for spec in [spec for spec in spec_st if (spec.stats.channel[-1] == 'H')]:
        statId = '%s %s' % (spec.id, spec.stats.instrtype)
        try:
            par = sourcepar[statId]
        except KeyError:
            continue
        t_star = par['t_star']
        deltaf = spec.stats.delta
        freq = spec.get_freq()
        rho = config.rho
        vs = config.hypo.vs * 1000.
        rps = config.rps
        # Data is in moment units. Let's put it back to displacement units,
        # and derive it to velocity through multiplication by 2*pi*freq:
        # (2. is the free-surface amplification factor)
        coeff = 4 * np.pi * vs**3 * rho / (rps * 2.)
        data = (spec.data/coeff) * (2*np.pi*freq)
        # Correct data for attenuation:
        data *= np.exp(np.pi * t_star * freq)
        # Compute the energy integral, up to fmax:
        fmax = config.max_freq_Er
        if fmax is not None:
            data[freq > fmax] = 0.
        Er = np.sum((data**2) * deltaf)
        # Apply corrections:
        # From Lancieri et al. (2012), eq. (3), the correction term is:
        #       8 * pi * r**2 * C**2 * rho * vs
        # We do not multiply by r**2, since data is already distance-corrected.
        # From Boatwright et al. (2002), eq. (2), C = <Fs>/(Fs * S),
        # where <Fs> is the average radiation pattern in the focal sphere
        # and Fs is the radiation pattern for the given angle between source
        # and receiver. Here we put <Fs>/Fs = 1, meaning that we rely on the
        # averaging between measurements at different stations, instead of
        # precise measurements at a single station.
        # S is the free-surface amplification factor, whihch we put = 2
        Er *= 8 * np.pi * (1./2.)**2 * rho * vs

        # Let's now compute finite bandwidth correction
        # expressed as the ratio R between the estimated energy
        # and the true energy (Di Bona & Rovelli 1988):
        if fmax is None:
            fmax = freq[-1]
        fc = par['fc']
        R = 2./np.pi * (np.arctan2(fmax, fc) - (fmax/fc)/(1+(fmax/fc)**2.))

        # And apply the correction
        Er /= R

        # Store in the parameter dictionary
        par['Er'] = Er
