# -*- coding: utf8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Compute radiated energy from spectral integration.

:copyright:
    2012 Claudio Satriano <satriano@ipgp.fr>

    2013-2014 Claudio Satriano <satriano@ipgp.fr>,
              Emanuela Matrullo <matrullo@geologie.ens.fr>,
              Agnes Chounet <chounet@ipgp.fr>

    2015-2022 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
import logging
import numpy as np
logger = logging.getLogger(__name__.split('.')[-1])


def _spectral_integral(spec, t_star, fmax):
    """Compute spectral integral in eq. (3) from Lancieri et al. (2012)."""
    # Note: eq. (3) from Lancieri et al. (2012) is the same as
    # eq. (1) in Boatwright et al. (2002), but expressed in frequency,
    # instead of angular frequency (2pi factor).
    deltaf = spec.stats.delta
    freq = spec.get_freq()
    # Data is in moment units. Let's put it back to displacement units,
    # and derive it to velocity through multiplication by 2*pi*freq:
    # (2. is the free-surface amplification factor)
    data = (spec.data/spec.coeff) * (2*np.pi*freq)
    # Correct data for attenuation:
    data *= np.exp(np.pi * t_star * freq)
    # Compute the energy integral, up to fmax:
    if fmax is not None:
        data[freq > fmax] = 0.
    integral = np.sum((data**2) * deltaf)
    return integral


def _radiated_energy_coefficient(rho, vel):
    """Compute coefficient in eq. (3) from Lancieri et al. (2012)."""
    # Note: eq. (3) from Lancieri et al. (2012) is the same as
    # eq. (1) in Boatwright et al. (2002), but expressed in frequency,
    # instead of angular frequency (2pi factor).
    # In the original eq. (3) from Lancieri et al. (2012), eq. (3),
    # the correction term is:
    #       8 * pi * r**2 * C**2 * rho * vs
    # We do not multiply by r**2, since data is already distance-corrected.
    # From Boatwright et al. (2002), eq. (2), C = <Fs>/(Fs * S),
    # where <Fs> is the average radiation pattern in the focal sphere
    # and Fs is the radiation pattern for the given angle between source
    # and receiver. Here we put <Fs>/Fs = 1, meaning that we rely on the
    # averaging between measurements at different stations, instead of
    # precise measurements at a single station.
    # S is the free-surface amplification factor, whihch we put = 2
    coeff = 8 * np.pi * (1./2.)**2 * rho * vel
    return coeff


def _finite_bandwidth_correction(spec, fc, fmax):
    """
    Compute finite bandwidth correction.

    Expressed as the ratio R between the estimated energy
    and the true energy (Di Bona & Rovelli 1988)
    """
    if fmax is None:
        fmax = spec.get_freq()[-1]
    R = 2./np.pi * (np.arctan2(fmax, fc) - (fmax/fc)/(1+(fmax/fc)**2.))
    return R


def radiated_energy(config, spec_st, specnoise_st, sourcepar):
    """Compute radiated energy, using eq. (3) in Lancieri et al. (2012)."""
    logger.info('Computing radiated energy...')
    if config.wave_type == 'P':
        logger.warning(
            'Warning: computing radiated energy from P waves might lead to '
            'an underestimation')
    # Select specids with channel code: '??H'
    spec_ids = [spec.id for spec in spec_st if spec.id[-1] == 'H']
    for spec_id in spec_ids:
        spec = spec_st.select(id=spec_id)[0]
        specnoise = specnoise_st.select(id=spec_id)[0]

        statId = '{} {}'.format(spec_id, spec.stats.instrtype)
        try:
            par = sourcepar.station_parameters[statId]
        except KeyError:
            continue

        t_star = par['t_star']
        fc = par['fc']
        fmax = config.max_freq_Er
        rho = config.rho
        if config.wave_type == 'P':
            vel = config.hypo.vp * 1000.
        elif config.wave_type in ['S', 'SV', 'SH']:
            vel = config.hypo.vs * 1000.

        # Compute signal and noise integrals and subtract noise from signal,
        # under the hypothesis that energy is additive and noise is stationary
        signal_integral = _spectral_integral(spec, t_star, fmax)
        noise_integral = _spectral_integral(specnoise, t_star, fmax)
        coeff = _radiated_energy_coefficient(rho, vel)
        Er = coeff * (signal_integral - noise_integral)
        if Er < 0:
            msg = '{}: noise energy is larger than signal energy: '.format(
                statId)
            msg += 'skipping spectrum.'
            logger.warning(msg)
            par['Er'] = np.nan
            continue

        R = _finite_bandwidth_correction(spec, fc, fmax)
        Er /= R

        # Store in the parameter dictionary
        par['Er'] = Er
    logger.info('Computing radiated energy: done')
