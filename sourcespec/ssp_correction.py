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
import yaml
import math
import logging
from obspy import Stream
from sourcespec.ssp_setup import ssp_exit
from sourcespec.ssp_util import moment_to_mag, mag_to_moment
from scipy.interpolate import interp1d
logger = logging.getLogger(__name__.split('.')[-1])


def _correct_spectrum(spec, corr):
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
    return spec_corr


def _get_corrections(config, spec_st):
    """
    Get station corrections from residuals or kappa0 files.

    If residuals file is given, it is used to calculate the station
    corrections. Otherwise, kappa0 file is used to calculate the
    station corrections.

    :param config: Config object
    :param spec_st: Stream of spectra
    :return: Stream of station corrections
    """
    res_filepath = config.residuals_filepath
    if res_filepath is not None:
        logger.info(f'Loading residuals from {res_filepath}')
        with open(res_filepath, 'rb') as fp:
            return pickle.load(fp)
    kappa0_filepath = config.station_kappa0_filepath
    if kappa0_filepath is not None:
        logger.info(f'Loading kappa0 from {kappa0_filepath}')
        with open(kappa0_filepath, 'rb') as fk0:
            kappa0_data = yaml.safe_load(fk0)
        corr_stream = Stream()
        for spec in spec_st:
            spec_id = spec.id
            try:
                # reading kappa0 value for the given station
                sta_id = '.'.join(spec_id.split('.')[:3])
                kappa0 = kappa0_data[sta_id]
            except KeyError:
                continue
            # calculating spectral correction
            corr = spec.copy()
            # kappa0 parametric correction in magnitude units
            corr.data_mag = -(
                2/3 * math.pi * corr.get_freq() *
                kappa0 * math.log10(math.e))
            corr.data = mag_to_moment(corr.data_mag)
            corr_stream.append(corr)
        return corr_stream
    return Stream()


def station_correction(spec_st, config):
    """
    Correct spectra using station-average residuals or kappa0.

    Residuals are obtained from a previous run.
    Kappa0 is obtained from a file.
    """
    H_specs = [spec for spec in spec_st if (spec.stats.channel[-1] == 'H')]
    try:
        corrections = _get_corrections(config, H_specs)
    except Exception as e:
        logger.error(str(e))
        ssp_exit(1)
    for spec in H_specs:
        spec_id = spec.id
        try:
            corr = corrections.select(id=spec_id)[0]
        except IndexError:
            logger.warning(f'No correction found for {spec_id}')
            continue
        spec_corr = _correct_spectrum(spec, corr)
        spec_st.append(spec_corr)
    return spec_st