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


def station_correction(spec_st, config):
    """
    Correct spectra using station-average residuals.

    Residuals are obtained from a previous run.
    """
    res_filepath = config.residuals_filepath
    if res_filepath is None:
        return spec_st
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
        freq = spec.freq
        corr_interp = corr.copy()
        corr_interp.freq = freq
        # corr_interp.data must exist, so we create it with zeros
        corr_interp.data = np.zeros_like(freq)
        # interpolate the correction to the same frequencies as the spectrum
        f = interp1d(
            corr.freq, corr.data_mag, fill_value='extrapolate')
        corr_interp.data_mag = f(freq)
        spec_corr = spec.copy()
        # uncorrected spectrum will have component name 'h'
        spec.stats.channel = f'{spec.stats.channel[:-1]}h'
        try:
            spec_corr.data_mag -= corr_interp.data_mag
        except ValueError as msg:
            logger.error(f'Cannot correct spectrum {spec.id}: {msg}')
            continue
        # interpolate the corrected data_mag to logspaced frequencies
        f = interp1d(freq, spec_corr.data_mag, fill_value='extrapolate')
        spec_corr.data_mag_logspaced = f(spec_corr.freq_logspaced)
        # convert mag to moment
        spec_corr.data = mag_to_moment(spec_corr.data_mag)
        spec_corr.data_logspaced = mag_to_moment(spec_corr.data_mag_logspaced)
        spec_st.append(spec_corr)
        fmin = freq.min()
        fmax = freq.max()
        logger.info(
            f'{spec_corr.id}: corrected, frequency range is: '
            f'{fmin:.2f} {fmax:.2f} Hz')
    return spec_st
