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
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import logging
import numpy as np
from sourcespec.ssp_setup import ssp_exit
from sourcespec.ssp_util import mag_to_moment
logger = logging.getLogger(__name__.split('.')[-1])


def _remove_outliers(values, weights, nstd):
    """Remove extreme values that are larger than nstd*std."""
    _values = values.copy()
    nval = len(_values)
    outliers = np.zeros(nval).astype(bool)
    for _ in range(nval):
        # find the extreme value
        vmean = np.nanmean(_values)
        argextreme = np.nanargmax(np.abs(_values-vmean))
        extreme = values[argextreme]
        # trim the extreme value and compute statistics
        # for the remaining values
        values_trim = np.delete(_values, argextreme)
        vmean = np.nanmean(values_trim)
        vstd = np.nanstd(values_trim)
        if np.abs(extreme - vmean) > nstd*vstd:
            _values[argextreme] = np.nan
            outliers[argextreme] = True
        else:
            # if no outlier is found, break
            break
        if np.sum(~outliers) == 3:
            # if only 3 non-outliers remain, break
            break
    values = values[~outliers]
    if weights is not None:
        weights = weights[~outliers]
    return values, weights, outliers


def _avg_and_std(values, errors=None, logarithmic=False, nstd=None):
    """
    Return the average and standard deviation.

    Optionally:
    - errors can be specfied for weighted statistics
    - logarithmic average and standard deviation
    - cutoff outlier values that are larger than nstd*std
    """
    average = std = np.nan
    outliers = np.zeros_like(values).astype(bool)
    if len(values) == 0:
        return average, std, outliers
    if errors is not None and len(errors) == 0:
        return average, std, outliers
    if np.all(np.isnan(values)):
        return average, std, outliers
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
    if nstd is not None:
        values, weights, outliers = _remove_outliers(values, weights, nstd)
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
        return log_average, (minus, plus), outliers
    else:
        return average, std, outliers


def _M0_avg_and_std(Mw_mean, Mw_error):
    """Compute average and standard deviation for Mo starting from Mw."""
    Mo_mean = mag_to_moment(Mw_mean)
    Mo_min = mag_to_moment(Mw_mean - Mw_error)
    Mo_max = mag_to_moment(Mw_mean + Mw_error)
    return Mo_mean, (Mo_mean - Mo_min, Mo_max - Mo_mean)


def compute_averages(sourcepar):
    """Compute average source parameters, find outliers"""
    if len(sourcepar) == 0:
        logger.info('No source parameter calculated')
        ssp_exit()

    means = dict()
    errors = dict()
    means_weight = dict()
    errors_weight = dict()

    # Mw
    Mw_values = sourcepar.value_array('Mw')
    Mw_err = sourcepar.error_array('Mw')
    means['Mw'], errors['Mw'], _ = _avg_and_std(Mw_values, nstd=3)
    means_weight['Mw'], errors_weight['Mw'], outliers = \
        _avg_and_std(Mw_values, errors=Mw_err, nstd=3)
    sourcepar.set_outliers('Mw', outliers)

    # Mo (N.m)
    means['Mo'], errors['Mo'] = _M0_avg_and_std(means['Mw'], errors['Mw'])
    means_weight['Mo'], errors_weight['Mo'] = \
        _M0_avg_and_std(means_weight['Mw'], errors_weight['Mw'])

    # fc (Hz)
    fc_values = sourcepar.value_array('fc')
    fc_err = sourcepar.error_array('fc')
    means['fc'], errors['fc'], _ =\
        _avg_and_std(fc_values, logarithmic=True, nstd=2)
    means_weight['fc'], errors_weight['fc'], outliers = \
        _avg_and_std(fc_values, errors=fc_err, logarithmic=True, nstd=2)
    sourcepar.set_outliers('fc', outliers)

    # t_star (s)
    t_star_values = sourcepar.value_array('t_star')
    t_star_err = sourcepar.error_array('t_star')
    means['t_star'], errors['t_star'], _ = _avg_and_std(t_star_values, nstd=2)
    means_weight['t_star'], errors_weight['t_star'], outliers = \
        _avg_and_std(t_star_values, errors=t_star_err, nstd=2)
    sourcepar.set_outliers('t_star', outliers)

    # ra, radius (meters)
    ra_values = sourcepar.value_array('ra')
    ra_err = sourcepar.error_array('ra')
    means['ra'], errors['ra'], _ =\
        _avg_and_std(ra_values, logarithmic=True, nstd=2)
    means_weight['ra'], errors_weight['ra'], outliers = \
        _avg_and_std(ra_values, errors=ra_err, logarithmic=True, nstd=2)
    sourcepar.set_outliers('ra', outliers)

    # bsd, Brune stress drop (MPa)
    bsd_values = sourcepar.value_array('bsd')
    bsd_err = sourcepar.error_array('bsd')
    means['bsd'], errors['bsd'], _ =\
        _avg_and_std(bsd_values, logarithmic=True, nstd=2)
    means_weight['bsd'], errors_weight['bsd'], outliers = \
        _avg_and_std(bsd_values, errors=bsd_err, logarithmic=True, nstd=2)
    sourcepar.set_outliers('bsd', outliers)

    # Quality factor
    Qo_values = sourcepar.value_array('Qo')
    Qo_err = sourcepar.error_array('Qo')
    Qo_values[np.isinf(Qo_values)] = np.nan
    means['Qo'], errors['Qo'], _ = _avg_and_std(Qo_values, nstd=2)
    cond_err = np.logical_or(np.isinf(Qo_err[:, 0]), np.isinf(Qo_err[:, 1]))
    Qo_values[cond_err] = np.nan
    Qo_err[cond_err] = np.nan
    means_weight['Qo'], errors_weight['Qo'], outliers = \
        _avg_and_std(Qo_values, errors=Qo_err, nstd=2)
    sourcepar.set_outliers('Qo', outliers)

    # Ml
    Ml_values = sourcepar.value_array('Ml')
    means['Ml'], errors['Ml'], outliers = _avg_and_std(Ml_values, nstd=3)
    sourcepar.set_outliers('Ml', outliers)

    # Er (N.m)
    Er_values = sourcepar.value_array('Er')
    means['Er'], errors['Er'], outliers =\
        _avg_and_std(Er_values, logarithmic=True, nstd=2)
    sourcepar.set_outliers('Er', outliers)

    sourcepar['means'] = means
    sourcepar['errors'] = errors
    sourcepar['means_weight'] = means_weight
    sourcepar['errors_weight'] = errors_weight

    params_name = ('Mw', 'fc', 't_star')
    sourcepar_mean = dict(
        zip(params_name, [means['Mw'], means['fc'], means['t_star']]))
    logger.info('params_mean: {}'.format(sourcepar_mean))
    sourcepar_mean_weight = dict(
        zip(params_name,
            [means_weight['Mw'], means_weight['fc'], means_weight['t_star']]))
    logger.info('params_mean_weighted: {}'.format(sourcepar_mean_weight))

    return sourcepar_mean
