# -*- coding: utf8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Spectral station correction calculated from ssp_residuals.

:copyright:
    2013-2014 Claudio Satriano <satriano@ipgp.fr>,
              Agnes Chounet <chounet@ipgp.fr>

    2015-2026 Claudio Satriano <satriano@ipgp.fr>
              Kris Vanneste <kris.vanneste@oma.be>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
import logging
import numpy as np
from .spectrum import read_spectra
from .ssp_util import mag_to_moment
from .setup import config, ssp_exit
logger = logging.getLogger(__name__.rsplit('.', maxsplit=1)[-1])


def station_correction(spec_st):
    """
    Correct spectra using station-average residuals.

    Residuals are obtained from a previous run.

    Parameters
    ----------
    spec_st : SpectrumStream or list of Spectrum
        List of spectra to be corrected. Corrected spectra are appended
        to the list (component code of uncorrected spectra is renamed to 'h').
    """
    res_filepath = config.residuals_filepath
    if res_filepath is None:
        return
    try:
        residual = read_spectra(res_filepath)
    except Exception as msg:
        ssp_exit(msg)

    H_specs = [
        spec for spec in spec_st
        if spec.stats.channel[-1] == 'H' and not spec.stats.ignore
    ]
    for spec in H_specs:
        try:
            corr = residual.select(id=spec.id)[0]
        except IndexError:
            continue
        # Define common frequency range for the correction and the spectrum
        freq_min = max(spec.freq.min(), corr.freq.min())
        freq_max = min(spec.freq.max(), corr.freq.max())
        # Note that frequency range of corrected spectrum must not change,
        # otherwise it will be out of sync with the noise spectrum
        # and with the weight used in the inversion
        # Instead, we use NaN values outside the common frequency range
        corr_interp = corr.copy()
        corr_interp.interp_data_to_new_freq(spec.freq)
        corr_interp.data_mag[corr_interp.freq < freq_min] = np.nan
        corr_interp.data_mag[corr_interp.freq > freq_max] = np.nan
        # Copy spectrum before correction
        spec_corr = spec.copy()
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
        spec_corr.make_logspaced_from_linear(which='data_mag')
        # Convert mag to moment
        spec_corr.data = mag_to_moment(spec_corr.data_mag)
        spec_corr.data_logspaced = mag_to_moment(spec_corr.data_mag_logspaced)
        spec_st.append(spec_corr)
        # Find frequency range of the corrected spectrum and log the result
        nan_idxs = np.isnan(spec_corr.data_mag)
        fmin = spec_corr.freq[~nan_idxs].min()
        fmax = spec_corr.freq[~nan_idxs].max()
        logger.info(
            f'{spec_corr.id}: corrected, frequency range is: '
            f'{fmin:.2f} {fmax:.2f} Hz')
