# -*- coding: utf-8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Post processing of station source parameters.

:copyright:
    2012-2025 Claudio Satriano <satriano@ipgp.fr>
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
logger = logging.getLogger(__name__.rsplit('.', maxsplit=1)[-1])


def _avg_and_std(values, errors=None, logarithmic=False):
    """
    Return the average and standard deviation.

    Optionally:
    - errors can be specified for weighted statistics
    - logarithmic average and standard deviation
    """
    average = std = np.nan
    if len(values) == 0:
        return average, np.array((std, std))
    if errors is not None and len(errors) == 0:
        return average, np.array((std, std))
    if np.all(np.isnan(values)):
        return average, np.array((std, std))
    weights = _weights(values, errors, logarithmic)
    if logarithmic:
        values = np.log10(values)
    notnan = ~np.isnan(values)
    if weights is not None:
        notnan = np.logical_and(notnan, ~np.isnan(weights))
        weights = weights[notnan]
    values = values[notnan]
    average = np.average(values, weights=weights)
    variance = np.average((values - average)**2, weights=weights)
    std = np.sqrt(variance)
    if not logarithmic:
        return average, np.array((std, std))
    log_average = 10.**average
    minus = log_average - 10.**(average - std)
    plus = 10.**(average + std) - log_average
    return log_average, np.array((minus, plus))


def _weights(values, errors=None, logarithmic=False):
    """Compute weights for weighted statistics."""
    if errors is None:
        return None
    # negative errors should not happen
    errors[errors < 0] = 0
    values_minus = values - errors[:, 0]
    values_plus = values + errors[:, 1]
    if logarithmic:
        # compute the width of the error bar in log10 units
        # replace negative left values with 1/10 of the central value
        values_minus[values_minus <= 0] = values[values_minus <= 0] / 10
        values_log_minus = np.log10(values_minus)
        values_log_plus = np.log10(values_plus)
        errors_width = values_log_plus - values_log_minus
    else:
        # compute the width of the error bar in linear units
        errors_width = values_plus - values_minus
    try:
        # fix for infinite weight (zero error width)
        errors_width[errors_width == 0] =\
            np.nanmin(errors_width[errors_width > 0])
    except ValueError:
        # if all errors are zero, return ones
        return np.ones_like(errors_width)
    return 1. / (errors_width**2.)


def _normal_confidence_level(n_sigma):
    """
    Compute the confidence level of a normal (Gaussian) distribution
    between -n_sigma and +n_sigma.
    """
    def gauss(x):
        return norm.pdf(x, 0, 1)
    confidence, _ = quad(gauss, -n_sigma, n_sigma)
    return np.round(confidence * 100, 2)


def _percentiles(
        values, low_percentage=25, mid_percentage=50, up_percentage=75):
    """Compute lower, mid and upper percentiles."""
    if len(values) == 0:
        return np.nan, np.nan, np.nan
    low_percentile, mid_percentile, up_percentile =\
        np.nanpercentile(
            values, (low_percentage, mid_percentage, up_percentage))
    return low_percentile, mid_percentile, up_percentile


def _param_summary_statistics(
        config, sspec_output, param_id, name, format_spec, units=None,
        logarithmic=False):
    """Compute summary statistics for one spectral parameter."""
    nIQR = config.nIQR
    summary = SummarySpectralParameter(
        param_id=param_id, name=name, format_spec=format_spec, units=units)
    sspec_output.find_outliers(param_id, n=nIQR)
    values = sspec_output.value_array(param_id, filter_outliers=True)
    errors = sspec_output.error_array(param_id, filter_outliers=True)
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
        stat_type='mean', value=mean_value,
        lower_uncertainty=mean_error[0],
        upper_uncertainty=mean_error[1],
        confidence_level=conf_level, nobs=nobs)
    if not np.all(np.isnan(errors)):
        # weighted mean (only if errors are defined)
        wmean_value, wmean_error = _avg_and_std(
            values, errors, logarithmic=logarithmic)
    else:
        # else use an unweighted mean
        logger.info(
            f'{param_id} values have no uncertainty: weighted mean cannot be '
            'computed. Using unweighted mean instead.')
        wmean_value, wmean_error = mean_value, mean_error
    summary.weighted_mean = SummaryStatistics(
        stat_type='weighted_mean', value=wmean_value,
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
        stat_type='percentiles', value=mid_pctl,
        lower_uncertainty=mid_pctl - low_pctl,
        upper_uncertainty=up_pctl - mid_pctl,
        confidence_level=conf_level,
        lower_percentage=config.lower_percentage,
        mid_percentage=config.mid_percentage,
        upper_percentage=config.upper_percentage,
        nobs=nobs)
    return summary


def compute_summary_statistics(config, sspec_output):
    """Compute summary statistics from station spectral parameters."""
    logger.info('Computing summary statistics...')
    if len(sspec_output.station_parameters) == 0:
        logger.info('No source parameter calculated')
        ssp_exit()

    sspec_output.summary_spectral_parameters.reference_statistics =\
        config.reference_statistics

    # Mw
    sspec_output.summary_spectral_parameters.Mw =\
        _param_summary_statistics(
            config, sspec_output,
            param_id='Mw', name='moment magnitude', format_spec='{:.2f}',
            logarithmic=False
        )

    # Mo (N·m)
    sspec_output.summary_spectral_parameters.Mo =\
        _param_summary_statistics(
            config, sspec_output,
            param_id='Mo', name='seismic moment', units='N·m',
            format_spec='{:.3e}', logarithmic=True
        )

    # fc (Hz)
    sspec_output.summary_spectral_parameters.fc =\
        _param_summary_statistics(
            config, sspec_output,
            param_id='fc', name='corner frequency', units='Hz',
            format_spec='{:.3f}', logarithmic=True
        )

    # t_star (s)
    sspec_output.summary_spectral_parameters.t_star =\
        _param_summary_statistics(
            config, sspec_output,
            param_id='t_star', name='t-star', units='s', format_spec='{:.3f}',
            logarithmic=False
        )

    # radius (meters)
    sspec_output.summary_spectral_parameters.radius =\
        _param_summary_statistics(
            config, sspec_output,
            param_id='radius', name='source radius', units='m',
            format_spec='{:.3f}', logarithmic=True
        )

    # static stress drop (MPa)
    sspec_output.summary_spectral_parameters.ssd =\
        _param_summary_statistics(
            config, sspec_output,
            param_id='ssd', name='static stress drop',
            units='MPa', format_spec='{:.3e}',
            logarithmic=True
        )

    # Quality factor
    sspec_output.summary_spectral_parameters.Qo =\
        _param_summary_statistics(
            config, sspec_output,
            param_id='Qo', name='quality factor', format_spec='{:.1f}',
            logarithmic=False
        )

    # Er (N·m)
    sspec_output.summary_spectral_parameters.Er =\
        _param_summary_statistics(
            config, sspec_output,
            param_id='Er', name='radiated energy', units='N·m',
            format_spec='{:.3e}', logarithmic=True
        )

    # Apparent stress (MPa)
    sspec_output.summary_spectral_parameters.sigma_a =\
        _param_summary_statistics(
            config, sspec_output,
            param_id='sigma_a', name='apparent stress', units='MPa',
            format_spec='{:.3e}', logarithmic=True
        )

    # Ml
    if config.compute_local_magnitude:
        sspec_output.summary_spectral_parameters.Ml =\
            _param_summary_statistics(
                config, sspec_output,
                param_id='Ml', name='local magnitude', format_spec='{:.2f}',
                logarithmic=False
            )

    # Log info on summary statistics for each parameter
    # Note that numpy.float64 are cast to float to display them in the log
    # in a more readable way
    params_name = ('Mw', 'fc', 't_star')
    means = sspec_output.mean_values()
    sourcepar_mean = {par: float(means[par]) for par in params_name}
    logger.info(f'params_mean: {sourcepar_mean}')
    means_weight = sspec_output.weighted_mean_values()
    sourcepar_mean_weight = {
        par: float(means_weight[par]) for par in params_name}
    logger.info(f'params_mean_weighted: {sourcepar_mean_weight}')
    percentiles = sspec_output.percentiles_values()
    sourcepar_percentiles = {
        par: float(percentiles[par]) for par in params_name}
    logger.info(f'params_percentiles: {sourcepar_percentiles}')
    logger.info('Computing summary statistics: done')
    logger.info('---------------------------------------------------')
