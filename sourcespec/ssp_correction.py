# -*- coding: utf8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Spectral station correction calculated from ssp_residuals.

:copyright:
    2013-2014 Claudio Satriano <satriano@ipgp.fr>,
              Agnes Chounet <chounet@ipgp.fr>

    2015-2025 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
import logging
import numpy as np
from scipy.interpolate import interp1d
from sourcespec.spectrum import read_spectra
from sourcespec.ssp_util import mag_to_moment
from sourcespec.ssp_setup import ssp_exit
logger = logging.getLogger(__name__.rsplit('.', maxsplit=1)[-1])


def _interpolate_spectrum_to_new_freq_range(spec, new_freq):
    """
    Interpolate spectrum to a new frequency range.

    Parameters
    ----------
    spec : Spectrum
        Spectrum to be interpolated.
    new_freq : array_like
        New frequency range.

    Returns
    -------
    Spectrum
        Interpolated spectrum.
    """
    spec_interp = spec.copy()
    spec_interp.freq = new_freq
    # spec_interp.data must exist, so we create it with zeros
    spec_interp.data = np.zeros_like(new_freq)
    f = interp1d(
        spec.freq, spec.data_mag, fill_value='extrapolate')
    spec_interp.data_mag = f(new_freq)
    return spec_interp


def station_correction(spec_st, config):
    """
    Correct spectra using station-average residuals.

    Residuals are obtained from a previous run.

    Parameters
    ----------
    spec_st : SpectrumStream or list of Spectrum
        List of spectra to be corrected. Corrected spectra are appended
        to the list.
    config : Config
        Configuration object containing the residuals file path.
    """
    res_filepath = config.residuals_filepath
    if res_filepath is None:
        return
    try:
        residual = read_spectra(res_filepath)
    except Exception as msg:
        logger.error(msg)
        ssp_exit(1)

    H_specs = [spec for spec in spec_st if spec.stats.channel[-1] == 'H']
    for spec in H_specs:
        try:
            corr = residual.select(id=spec.id)[0]
        except IndexError:
            continue
        # Define a common frequency range for the correction and the spectrum
        freq_min = max(spec.freq.min(), corr.freq.min())
        freq_max = min(spec.freq.max(), corr.freq.max())
        delta_f = spec.stats.delta
        freq = np.arange(freq_min, freq_max, delta_f)
        # Interpolate residual to the common frequency range
        corr_interp = _interpolate_spectrum_to_new_freq_range(corr, freq)
        # Interpolate original spectrum before correction
        spec_corr = _interpolate_spectrum_to_new_freq_range(spec, freq)
        # Uncorrected spectrum will have component name 'h',
        # while corrected spectrum will have component name 'H'
        spec.stats.channel = f'{spec.stats.channel[:-1]}h'
        # Apply correction
        try:
            spec_corr.data_mag -= corr_interp.data_mag
        except ValueError as msg:
            logger.error(f'Cannot correct spectrum {spec.id}: {msg}')
            continue
        # Interpolate the corrected data_mag to logspaced frequencies
        # Note that spec_corr.freq_logspaced must not change, otherwise
        # it will be out of sync with the weight used in the inversion.
        # Instead, we set it to NaN outside the valid frequency range,
        # through the parameter bounds_error=False.
        f = interp1d(spec_corr.freq, spec_corr.data_mag, bounds_error=False)
        spec_corr.data_mag_logspaced = f(spec_corr.freq_logspaced)
        # Convert mag to moment
        spec_corr.data = mag_to_moment(spec_corr.data_mag)
        spec_corr.data_logspaced = mag_to_moment(spec_corr.data_mag_logspaced)
        spec_st.append(spec_corr)
        fmin = freq.min()
        fmax = freq.max()
        logger.info(
            f'{spec_corr.id}: corrected, frequency range is: '
            f'{fmin:.2f} {fmax:.2f} Hz')
