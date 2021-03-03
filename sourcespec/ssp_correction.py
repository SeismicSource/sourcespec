# -*- coding: utf8 -*-
"""
Spectral station correction calculated from ssp_residuals.

:copyright:
    2013-2014 Claudio Satriano <satriano@ipgp.fr>,
              Agnes Chounet <chounet@ipgp.fr>

    2015-2021 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement, Version 2.1
    (http://www.cecill.info/index.en.html)
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

try:
    import cPickle as pickle
except ImportError:
    import pickle
import logging
from sourcespec.ssp_util import moment_to_mag, mag_to_moment
from scipy.interpolate import interp1d
logger = logging.getLogger(__name__.split('.')[-1])


def station_correction(spec_st, config):
    """
    Correct spectra using station-average residuals.

    Residuals are obtained from a previous run.
    """
    res_filepath = config.residuals_filepath
    if res_filepath is None:
        msg = "'-C' option set, but 'residuals_filepath' not specified "
        msg += "in config file: ignoring station correction"
        logger.warning(msg)
        return spec_st
    with open(res_filepath, 'rb') as fp:
        residual = pickle.load(fp)

    for spec in [spec for spec in spec_st if (spec.stats.channel[-1] == 'H')]:
        station = spec.stats.station
        if station in set(x.stats.station for x in residual):
            # apply correction
            corr = residual.select(station=station)[0]
            freq = spec.get_freq()
            fmin = freq.min()
            fmax = freq.max()
            corr = corr.slice(fmin, fmax)
            corr.data_mag = moment_to_mag(corr.data)
            spec.data_mag -= corr.data_mag
            # interpolate the corrected data_mag to log frequencies
            f = interp1d(freq, spec.data_mag, fill_value='extrapolate')
            spec.data_log_mag = f(spec.freq_log)
            # convert mag to moment
            spec.data = mag_to_moment(spec.data_mag)
            spec.data_log = mag_to_moment(spec.data_log_mag)

            logger.info('%s corrected, frequency range is: %f %f'
                        % (spec.id, fmin, fmax))
    return spec_st
