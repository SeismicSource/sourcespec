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
from scipy.stats import norm
from scipy.integrate import quad
from sourcespec.ssp_setup import ssp_exit
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
        return log_average, np.array((minus, plus))
    else:
        return average, np.array((std, std))


def _normal_confidence_level(n_sigma):
    """
    Compute the confidence level of a normal (Gaussian) distribution
    between -n_sigma and +n_sigma.
    """
    def gauss(x):
        return norm.pdf(x, 0, 1)
    confidence, _ = quad(gauss, -n_sigma, n_sigma)
    return np.round(confidence*100, 2)


def _percentiles(
        values, low_percentage=25, mid_percentage=50, up_percentage=75):
    """Compute lower, mid and upper percentiles."""
    low_percentile, mid_percentile, up_percentile =\
        np.nanpercentile(
            values, (low_percentage, mid_percentage, up_percentage))
    return low_percentile, mid_percentile, up_percentile


def _param_summary_statistics(
        config, sspec_output, id, name, format, units=None, logarithmic=False):
    """Compute summary statistics for one spectral parameter."""
    nIQR = config.nIQR
    summary = SummarySpectralParameter(
        id=id, name=name, format=format, units=units)
    sspec_output.find_outliers(id, n=nIQR)
    values = sspec_output.value_array(id, filter_outliers=True)
    errors = sspec_output.error_array(id, filter_outliers=True)
    # put to NaN infinite values and values whose errors are infinite
    values[np.isinf(values)] = np.nan
    _cond_err = np.logical_or(np.isinf(errors[:, 0]), np.isinf(errors[:, 1]))
    values[_cond_err] = np.nan
    errors[_cond_err] = np.nan
    # only count non-NaN values
    nobs = len(values[~np.isnan(values)])
    # mean
    mean_value, mean_error = _avg_and_std(values, logarithmic=logarithmic)
    mean_error *= config.n_sigma
    conf_level = _normal_confidence_level(config.n_sigma)
    summary.mean = SummaryStatistics(
        type='mean', value=mean_value,
        lower_uncertainty=mean_error[0],
        upper_uncertainty=mean_error[1],
        confidence_level=conf_level, nobs=nobs)
    # weighted mean (only if errors are defined)
    if not np.all(np.isnan(errors)):
        wmean_value, wmean_error = _avg_and_std(
            values, errors, logarithmic=logarithmic)
        wmean_error *= config.n_sigma
        summary.weighted_mean = SummaryStatistics(
            type='weighted_mean', value=wmean_value,
            lower_uncertainty=wmean_error[0],
            upper_uncertainty=wmean_error[1],
            confidence_level=conf_level, nobs=nobs)
    # percentiles
    low_pctl, mid_pctl, up_pctl = _percentiles(
        values, config.lower_percentage, config.mid_percentage,
        config.upper_percentage)
    conf_level = round(
        (config.upper_percentage - config.lower_percentage), 2)
    summary.percentiles = SummaryStatistics(
        type='percentiles', value=mid_pctl,
        lower_uncertainty=mid_pctl-low_pctl,
        upper_uncertainty=up_pctl-mid_pctl,
        confidence_level=conf_level,
        lower_percentage=config.lower_percentage,
        mid_percentage=config.mid_percentage,
        upper_percentage=config.upper_percentage,
        nobs=nobs)
    return summary


def compute_summary_statistics(config, sspec_output):
    """Compute summary statistics from station spectral parameters."""
    if len(sspec_output.station_parameters) == 0:
        logger.info('No source parameter calculated')
        ssp_exit()

    sspec_output.summary_spectral_parameters.reference_statistics =\
        config.reference_statistics

    # Mw
    sspec_output.summary_spectral_parameters.Mw =\
        _param_summary_statistics(
            config, sspec_output,
            id='Mw', name='moment magnitude', format='{:.2f}',
            logarithmic=False
        )

    # Mo (N.m)
    sspec_output.summary_spectral_parameters.Mo =\
        _param_summary_statistics(
            config, sspec_output,
            id='Mo', name='seismic moment', units='N.m', format='{:.3e}',
            logarithmic=True
        )

    # fc (Hz)
    sspec_output.summary_spectral_parameters.fc =\
        _param_summary_statistics(
            config, sspec_output,
            id='fc', name='corner frequency', units='Hz', format='{:.3f}',
            logarithmic=True
        )

    # t_star (s)
    sspec_output.summary_spectral_parameters.t_star =\
        _param_summary_statistics(
            config, sspec_output,
            id='t_star', name='t-star', units='s', format='{:.3f}',
            logarithmic=False
        )

    # radius (meters)
    sspec_output.summary_spectral_parameters.radius =\
        _param_summary_statistics(
            config, sspec_output,
            id='radius', name='source radius', units='m', format='{:.3f}',
            logarithmic=True
        )

    # bsd, Brune stress drop (MPa)
    sspec_output.summary_spectral_parameters.bsd =\
        _param_summary_statistics(
            config, sspec_output,
            id='bsd', name='Brune stress drop', units='MPa', format='{:.3e}',
            logarithmic=True
        )

    # Quality factor
    sspec_output.summary_spectral_parameters.Qo =\
        _param_summary_statistics(
            config, sspec_output,
            id='Qo', name='quality factor', format='{:.1f}',
            logarithmic=False
        )

    # Er (N.m)
    sspec_output.summary_spectral_parameters.Er =\
        _param_summary_statistics(
            config, sspec_output,
            id='Er', name='radiated energy', units='N.m', format='{:.3e}',
            logarithmic=True
        )

    # Ml
    if config.compute_local_magnitude:
        sspec_output.summary_spectral_parameters.Ml =\
            _param_summary_statistics(
                config, sspec_output,
                id='Ml', name='local magnitude', format='{:.2f}',
                logarithmic=False
            )

    params_name = ('Mw', 'fc', 't_star')
    means = sspec_output.mean_values()
    sourcepar_mean = {par: means[par] for par in params_name}
    logger.info('params_mean: {}'.format(sourcepar_mean))
    means_weight = sspec_output.weighted_mean_values()
    sourcepar_mean_weight = {par: means_weight[par] for par in params_name}
    logger.info('params_mean_weighted: {}'.format(sourcepar_mean_weight))
    percentiles = sspec_output.percentiles_values()
    sourcepar_percentiles = {par: percentiles[par] for par in params_name}
    logger.info('params_percentiles: {}'.format(sourcepar_percentiles))
