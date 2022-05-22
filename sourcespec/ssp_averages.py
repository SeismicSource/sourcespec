# -*- coding: utf-8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Post processing of station source parameters.

:copyright:
    2012-2022 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
import logging
import numpy as np
from sourcespec.ssp_setup import ssp_exit
from sourcespec.ssp_util import mag_to_moment
logger = logging.getLogger(__name__.split('.')[-1])


def _avg_and_std(values, errors=None, logarithmic=False):
    """
    Return the average and standard deviation.

    Optionally:
    - errors can be specfied for weighted statistics
    - logarithmic average and standard deviation
    """
    average = std = np.nan
    if len(values) == 0:
        return average, std
    if errors is not None and len(errors) == 0:
        return average, std
    if np.all(np.isnan(values)):
        return average, std
    if errors is None:
        weights = None
    else:
        if logarithmic:
            # compute the width of the error bar in log10 units
            values_minus = values - errors[:, 0]
            # replace negative left values with 1/10 of the central value
            values_minus[values_minus <= 0] = values[values_minus <= 0]/10
            values_log_minus = np.log10(values_minus)
            values_plus = values + errors[:, 1]
            values_log_plus = np.log10(values_plus)
            errors_width = values_log_plus - values_log_minus
        else:
            # compute the width of the error bar
            errors_width = errors[:, 0] + errors[:, 1]
        # fix for infinite weight (zero error width)
        errors_width[errors_width == 0] =\
            np.nanmin(errors_width[errors_width > 0])
        weights = 1./(errors_width**2.)
    if logarithmic:
        values = np.log10(values)
    notnan = ~np.isnan(values)
    if weights is not None:
        notnan = np.logical_and(notnan, ~np.isnan(weights))
        weights = weights[notnan]
    values = values[notnan]
    average = np.average(values, weights=weights)
    variance = np.average((values-average)**2, weights=weights)
    std = np.sqrt(variance)
    if logarithmic:
        log_average = 10.**average
        minus = log_average - 10.**(average-std)
        plus = 10.**(average+std) - log_average
        return log_average, (minus, plus)
    else:
        return average, std


def _M0_avg_and_std(Mw_mean, Mw_error):
    """Compute average and standard deviation for Mo starting from Mw."""
    Mo_mean = mag_to_moment(Mw_mean)
    Mo_min = mag_to_moment(Mw_mean - Mw_error)
    Mo_max = mag_to_moment(Mw_mean + Mw_error)
    return Mo_mean, (Mo_mean - Mo_min, Mo_max - Mo_mean)


def compute_averages(config, sourcepar):
    """Compute average source parameters, find outliers"""
    if len(sourcepar.station_parameters) == 0:
        logger.info('No source parameter calculated')
        ssp_exit()

    means = dict()
    errors = dict()
    means_weight = dict()
    errors_weight = dict()

    nIQR = config.nIQR

    # Mw
    sourcepar.find_outliers('Mw', n=nIQR)
    Mw_values = sourcepar.value_array('Mw', filter_outliers=True)
    Mw_err = sourcepar.error_array('Mw', filter_outliers=True)
    means['Mw'], errors['Mw'] = _avg_and_std(Mw_values)
    means_weight['Mw'], errors_weight['Mw'] = _avg_and_std(Mw_values, Mw_err)

    # Mo (N.m)
    means['Mo'], errors['Mo'] = _M0_avg_and_std(means['Mw'], errors['Mw'])
    means_weight['Mo'], errors_weight['Mo'] = \
        _M0_avg_and_std(means_weight['Mw'], errors_weight['Mw'])

    # fc (Hz)
    sourcepar.find_outliers('fc', n=nIQR)
    fc_values = sourcepar.value_array('fc', filter_outliers=True)
    fc_err = sourcepar.error_array('fc', filter_outliers=True)
    means['fc'], errors['fc'] = _avg_and_std(fc_values, logarithmic=True)
    means_weight['fc'], errors_weight['fc'] =\
        _avg_and_std(fc_values, fc_err, logarithmic=True)

    # t_star (s)
    sourcepar.find_outliers('t_star', n=nIQR)
    t_star_values = sourcepar.value_array('t_star', filter_outliers=True)
    t_star_err = sourcepar.error_array('t_star', filter_outliers=True)
    means['t_star'], errors['t_star'] = _avg_and_std(t_star_values)
    means_weight['t_star'], errors_weight['t_star'] =\
        _avg_and_std(t_star_values, t_star_err)

    # ra, radius (meters)
    sourcepar.find_outliers('ra', n=nIQR)
    ra_values = sourcepar.value_array('ra', filter_outliers=True)
    ra_err = sourcepar.error_array('ra', filter_outliers=True)
    means['ra'], errors['ra'] = _avg_and_std(ra_values, logarithmic=True)
    means_weight['ra'], errors_weight['ra'] =\
        _avg_and_std(ra_values, ra_err, logarithmic=True)

    # bsd, Brune stress drop (MPa)
    sourcepar.find_outliers('bsd', n=nIQR)
    bsd_values = sourcepar.value_array('bsd', filter_outliers=True)
    bsd_err = sourcepar.error_array('bsd', filter_outliers=True)
    means['bsd'], errors['bsd'] = _avg_and_std(bsd_values, logarithmic=True)
    means_weight['bsd'], errors_weight['bsd'] =\
        _avg_and_std(bsd_values, bsd_err, logarithmic=True)

    # Quality factor
    sourcepar.find_outliers('Qo', n=nIQR)
    Qo_values = sourcepar.value_array('Qo', filter_outliers=True)
    Qo_err = sourcepar.error_array('Qo', filter_outliers=True)
    Qo_values[np.isinf(Qo_values)] = np.nan
    means['Qo'], errors['Qo'] = _avg_and_std(Qo_values)
    cond_err = np.logical_or(np.isinf(Qo_err[:, 0]), np.isinf(Qo_err[:, 1]))
    Qo_values[cond_err] = np.nan
    Qo_err[cond_err] = np.nan
    means_weight['Qo'], errors_weight['Qo'] = _avg_and_std(Qo_values, Qo_err)

    # Ml
    sourcepar.find_outliers('Ml', n=nIQR)
    Ml_values = sourcepar.value_array('Ml', filter_outliers=True)
    means['Ml'], errors['Ml'] = _avg_and_std(Ml_values)

    # Er (N.m)
    sourcepar.find_outliers('Er', n=nIQR)
    Er_values = sourcepar.value_array('Er', filter_outliers=True)
    means['Er'], errors['Er'] = _avg_and_std(Er_values, logarithmic=True)

    sourcepar.means = means
    sourcepar.errors = errors
    sourcepar.means_weight = means_weight
    sourcepar.errors_weight = errors_weight

    params_name = ('Mw', 'fc', 't_star')
    sourcepar_mean = dict(
        zip(params_name, [means['Mw'], means['fc'], means['t_star']]))
    logger.info('params_mean: {}'.format(sourcepar_mean))
    sourcepar_mean_weight = dict(
        zip(params_name,
            [means_weight['Mw'], means_weight['fc'], means_weight['t_star']]))
    logger.info('params_mean_weighted: {}'.format(sourcepar_mean_weight))
