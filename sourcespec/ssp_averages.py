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
from sourcespec.ssp_data_types import (
    SummarySpectralParameter, SummaryStatistics)
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


def compute_averages(config, sspec_output):
    """Compute average source parameters, find outliers"""
    if len(sspec_output.station_parameters) == 0:
        logger.info('No source parameter calculated')
        ssp_exit()

    nIQR = config.nIQR

    # Mw
    summary_Mw = SummarySpectralParameter(
        id='Mw', name='moment magnitude', format='{:.2f}')
    sspec_output.find_outliers('Mw', n=nIQR)
    Mw_values = sspec_output.value_array('Mw', filter_outliers=True)
    Mw_err = sspec_output.error_array('Mw', filter_outliers=True)
    nobs = len(Mw_values)
    value, error = _avg_and_std(Mw_values)
    summary_Mw.mean = SummaryStatistics(
        type='mean', value=value, uncertainty=error,
        confidence_level=68.2, nobs=nobs)
    value, error = _avg_and_std(Mw_values, Mw_err)
    summary_Mw.weighted_mean = SummaryStatistics(
        type='weighted_mean', value=value, uncertainty=error,
        confidence_level=68.2, nobs=nobs)
    sspec_output.summary_spectral_parameters.Mw = summary_Mw

    # Mo (N.m)
    summary_Mo = SummarySpectralParameter(
        id='Mo', name='seismic moment', units='N.m', format='{:.3e}')
    value, error = _M0_avg_and_std(
        summary_Mw.mean.value, summary_Mw.mean.uncertainty)
    summary_Mo.mean = SummaryStatistics(
        type='mean', value=value,
        lower_uncertainty=error[0], upper_uncertainty=error[1],
        confidence_level=68.2, nobs=nobs)
    value, error = _M0_avg_and_std(
        summary_Mw.weighted_mean.value, summary_Mw.weighted_mean.uncertainty)
    summary_Mo.weighted_mean = SummaryStatistics(
        type='weighted_mean', value=value,
        lower_uncertainty=error[0], upper_uncertainty=error[1],
        confidence_level=68.2, nobs=nobs)
    sspec_output.summary_spectral_parameters.Mo = summary_Mo

    # fc (Hz)
    summary_fc = SummarySpectralParameter(
        id='fc', name='corner frequency', units='Hz', format='{:.3f}')
    sspec_output.find_outliers('fc', n=nIQR)
    fc_values = sspec_output.value_array('fc', filter_outliers=True)
    fc_err = sspec_output.error_array('fc', filter_outliers=True)
    nobs = len(fc_values)
    value, error = _avg_and_std(fc_values, logarithmic=True)
    summary_fc.mean = SummaryStatistics(
        type='mean', value=value,
        lower_uncertainty=error[0], upper_uncertainty=error[1],
        confidence_level=68.2, nobs=nobs,
        message='logarithmic mean')
    value, error = _avg_and_std(fc_values, fc_err, logarithmic=True)
    summary_fc.weighted_mean = SummaryStatistics(
        type='weighted_mean', value=value,
        lower_uncertainty=error[0], upper_uncertainty=error[1],
        confidence_level=68.2, nobs=nobs,
        message='logarithmic mean')
    sspec_output.summary_spectral_parameters.fc = summary_fc

    # t_star (s)
    summary_t_star = SummarySpectralParameter(
        id='t_star', name='t*', units='s', format='{:.3f}')
    sspec_output.find_outliers('t_star', n=nIQR)
    t_star_values = sspec_output.value_array('t_star', filter_outliers=True)
    t_star_err = sspec_output.error_array('t_star', filter_outliers=True)
    nobs = len(t_star_values)
    value, error = _avg_and_std(t_star_values)
    summary_t_star.mean = SummaryStatistics(
        type='mean', value=value, uncertainty=error,
        confidence_level=68.2, nobs=nobs)
    value, error = _avg_and_std(t_star_values, t_star_err)
    summary_t_star.weighted_mean = SummaryStatistics(
        type='weighted_mean', value=value, uncertainty=error,
        confidence_level=68.2, nobs=nobs)
    sspec_output.summary_spectral_parameters.t_star = summary_t_star

    # radius (meters)
    summary_radius = SummarySpectralParameter(
        id='radius', name='source radius', units='m', format='{:.3f}')
    sspec_output.find_outliers('radius', n=nIQR)
    radius_values = sspec_output.value_array('radius', filter_outliers=True)
    radius_err = sspec_output.error_array('radius', filter_outliers=True)
    nobs = len(radius_values)
    value, error = _avg_and_std(radius_values, logarithmic=True)
    summary_radius.mean = SummaryStatistics(
        type='mean', value=value,
        lower_uncertainty=error[0], upper_uncertainty=error[1],
        confidence_level=68.2, nobs=nobs,
        message='logarithmic mean')
    value, error = _avg_and_std(radius_values, radius_err, logarithmic=True)
    summary_radius.weighted_mean = SummaryStatistics(
        type='weighted_mean', value=value,
        lower_uncertainty=error[0], upper_uncertainty=error[1],
        confidence_level=68.2, nobs=nobs,
        message='logarithmic mean')
    sspec_output.summary_spectral_parameters.radius = summary_radius

    # bsd, Brune stress drop (MPa)
    summary_bsd = SummarySpectralParameter(
        id='bsd', name='Brune stress drop', units='MPa', format='{:.3e}')
    sspec_output.find_outliers('bsd', n=nIQR)
    bsd_values = sspec_output.value_array('bsd', filter_outliers=True)
    bsd_err = sspec_output.error_array('bsd', filter_outliers=True)
    nobs = len(bsd_values)
    value, error = _avg_and_std(bsd_values, logarithmic=True)
    summary_bsd.mean = SummaryStatistics(
        type='mean', value=value,
        lower_uncertainty=error[0], upper_uncertainty=error[1],
        confidence_level=68.2, nobs=nobs,
        message='logarithmic mean')
    value, error = _avg_and_std(bsd_values, bsd_err, logarithmic=True)
    summary_bsd.weighted_mean = SummaryStatistics(
        type='weighted_mean', value=value,
        lower_uncertainty=error[0], upper_uncertainty=error[1],
        confidence_level=68.2, nobs=nobs,
        message='logarithmic mean')
    sspec_output.summary_spectral_parameters.bsd = summary_bsd

    # Quality factor
    summary_Qo = SummarySpectralParameter(
        id='Qo', name='quality factor', format='{:.1f}')
    sspec_output.find_outliers('Qo', n=nIQR)
    Qo_values = sspec_output.value_array('Qo', filter_outliers=True)
    Qo_err = sspec_output.error_array('Qo', filter_outliers=True)
    nobs = len(Qo_values)
    Qo_values[np.isinf(Qo_values)] = np.nan
    value, error = _avg_and_std(Qo_values)
    summary_Qo.mean = SummaryStatistics(
        type='mean', value=value, uncertainty=error,
        confidence_level=68.2, nobs=nobs)
    cond_err = np.logical_or(np.isinf(Qo_err[:, 0]), np.isinf(Qo_err[:, 1]))
    Qo_values[cond_err] = np.nan
    Qo_err[cond_err] = np.nan
    value, error = _avg_and_std(Qo_values, Qo_err)
    summary_Qo.weighted_mean = SummaryStatistics(
        type='weighted_mean', value=value, uncertainty=error,
        confidence_level=68.2, nobs=nobs)
    sspec_output.summary_spectral_parameters.Qo = summary_Qo

    # Er (N.m)
    summary_Er = SummarySpectralParameter(
        id='Er', name='radiated energy', units='N.m', format='{:.3e}')
    sspec_output.find_outliers('Er', n=nIQR)
    Er_values = sspec_output.value_array('Er', filter_outliers=True)
    nobs = len(Er_values)
    value, error = _avg_and_std(Er_values, logarithmic=True)
    summary_Er.mean = SummaryStatistics(
        type='mean', value=value,
        lower_uncertainty=error[0], upper_uncertainty=error[1],
        confidence_level=68.2, nobs=nobs,
        message='logarithmic mean')
    sspec_output.summary_spectral_parameters.Er = summary_Er

    # Ml
    if config.compute_local_magnitude:
        summary_Ml = SummarySpectralParameter(
            id='Ml', name='local magnitude', format='{:.2f}')
        sspec_output.find_outliers('Ml', n=nIQR)
        Ml_values = sspec_output.value_array('Ml', filter_outliers=True)
        nobs = len(Ml_values)
        value, error = _avg_and_std(Ml_values)
        summary_Ml.mean = SummaryStatistics(
            type='mean', value=value, uncertainty=error,
            confidence_level=68.2, nobs=nobs)
        sspec_output.summary_spectral_parameters.Ml = summary_Ml

    params_name = ('Mw', 'fc', 't_star')
    means = sspec_output.mean_values()
    sourcepar_mean = {par: means[par] for par in params_name}
    logger.info('params_mean: {}'.format(sourcepar_mean))
    means_weight = sspec_output.weighted_mean_values()
    sourcepar_mean_weight = {par: means_weight[par] for par in params_name}
    logger.info('params_mean_weighted: {}'.format(sourcepar_mean_weight))
