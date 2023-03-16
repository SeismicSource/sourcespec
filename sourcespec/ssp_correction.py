# -*- coding: utf8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Spectral station correction calculated from ssp_residuals.

:copyright:
    2013-2014 Claudio Satriano <satriano@ipgp.fr>,
              Agnes Chounet <chounet@ipgp.fr>

    2015-2023 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
import pickle
import logging
from sourcespec.ssp_util import moment_to_mag, mag_to_moment
from sourcespec.ssp_setup import ssp_exit
from scipy.interpolate import interp1d
logger = logging.getLogger(__name__.split('.')[-1])


def station_correction(spec_st, config):
    """
    Correct spectra using station-average residuals.

    Residuals are obtained from a previous run.
    """
    res_filepath = config.residuals_filepath
    if res_filepath is None:
        return spec_st
    try:
        with open(res_filepath, 'rb') as fp:
            residual = pickle.load(fp)
    except Exception as msg:
        logger.error(msg)
        ssp_exit(1)

    H_specs = [spec for spec in spec_st if (spec.stats.channel[-1] == 'H')]
    for spec in H_specs:
        try:
            corr = residual.select(id=spec.id)[0]
        except IndexError:
            continue
        freq = spec.get_freq()
        fmin = freq.min()
        fmax = freq.max()
        corr = corr.slice(fmin, fmax)
        corr.data_mag = moment_to_mag(corr.data)
        spec_corr = spec.copy()
        # uncorrected spectrum will have component name 'h'
        spec.stats.channel = f'{spec.stats.channel[:-1]}h'
        spec_corr.data_mag -= corr.data_mag
        # interpolate the corrected data_mag to log frequencies
        f = interp1d(freq, spec_corr.data_mag, fill_value='extrapolate')
        spec_corr.data_log_mag = f(spec_corr.freq_log)
        # convert mag to moment
        spec_corr.data = mag_to_moment(spec_corr.data_mag)
        spec_corr.data_log = mag_to_moment(spec_corr.data_log_mag)
        spec_st.append(spec_corr)
        logger.info(
            f'{spec_corr.id} corrected, frequency range is: '
            f'{fmin:.2f} {fmax:.2f} Hz')
    return spec_st
