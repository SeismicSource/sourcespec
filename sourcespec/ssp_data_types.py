# -*- coding: utf8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Classes for spectral inversion routines.

:copyright:
    2017-2022 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
import logging
import numpy as np
from collections import OrderedDict
logger = logging.getLogger(__name__.split('.')[-1])


class InitialValues():
    """Initial values for spectral inversion."""

    def __init__(self, Mw_0=None, fc_0=None, t_star_0=None):
        self.Mw_0 = Mw_0
        self.fc_0 = fc_0
        self.t_star_0 = t_star_0

    def __str__(self):
        """String representation."""
        s = 'Mw_0: %s; ' % round(self.Mw_0, 4)
        s += 'fc_0: %s; ' % round(self.fc_0, 4)
        s += 't_star_0: %s' % round(self.t_star_0, 4)
        return s

    def get_params0(self):
        return (self.Mw_0, self.fc_0, self.t_star_0)


class Bounds(object):
    """Bounds for bounded spectral inversion."""

    def __init__(self, config, spec, initial_values):
        self.config = config
        self.spec = spec
        self.hd = spec.stats.hypo_dist
        self.ini_values = initial_values
        self.Mw_min = self.Mw_max = None
        self._set_fc_min_max(config)
        if config.Qo_min_max is None:
            self.t_star_min, self.t_star_max =\
                self._check_minmax(config.t_star_min_max)
        else:
            self.t_star_min, self.t_star_max = self._Qo_to_t_star()
        self._fix_initial_values_t_star()

    def __str__(self):
        """String representation."""
        s = 'Mw: {}, {}; '.format(
            *[round(x, 4) if x is not None else x for x in self.bounds[0]])
        s += 'fc: {}, {}; '.format(
            *[round(x, 4) if x is not None else x for x in self.bounds[1]])
        s += 't_star: {}, {}'.format(
            *[round(x, 4) if x is not None else x for x in self.bounds[2]])
        return s

    def _set_fc_min_max(self, config):
        fc_0 = self.ini_values.fc_0
        if config.fc_min_max is None:
            # If no bound is given, set it to fc_0 +/- a decade
            scale = 10.  # a decade
            self.fc_min = fc_0/scale
            self.fc_max = fc_0*scale
        else:
            self.fc_min, self.fc_max = config.fc_min_max
        if self.fc_min > fc_0:
            logger.warning(
                '{} {}: fc_min ({}) larger than fc_0 ({}). '
                'Using fc_0 instead.'.format(
                    self.spec.id, self.spec.stats.instrtype,
                    self.fc_min, round(fc_0, 4))
            )
            self.fc_min = fc_0
        if self.fc_max < fc_0:
            logger.warning(
                '{} {}: fc_max ({}) smaller than fc_0 ({}). '
                'Using fc_0 instead.'.format(
                    self.spec.id, self.spec.stats.instrtype,
                    self.fc_max, round(fc_0, 4))
            )
            self.fc_max = fc_0

    def _check_minmax(self, minmax):
        if minmax is None:
            return (None, None)
        else:
            return minmax

    def _Qo_to_t_star(self):
        phase = self.config.wave_type[0]
        travel_time = self.spec.stats.travel_times[phase]
        t_star_bounds = travel_time/self.config.Qo_min_max
        return sorted(t_star_bounds)

    def _fix_initial_values_t_star(self):
        if self.ini_values.t_star_0 is not None:
            return
        if None in self.bounds[2]:
            return
        if self.t_star_min < self.ini_values.t_star_0 < self.t_star_max:
            return
        t_star_0 = (self.t_star_max + self.t_star_min) / 2.
        logger.warning(
            '{} {}: initial t_star value ({}) outside '
            'bounds. Using bound average ({})'.format(
                self.spec.id, self.spec.stats.instrtype,
                self.ini_values.t_star_0, round(t_star_0, 4))
        )
        self.ini_values.t_star_0 = t_star_0

    def __call__(self, **kwargs):
        """Interface for basin-hopping."""
        params = kwargs['x_new']
        params_min = np.array(
            (self.Mw_min, self.fc_min, self.t_star_min)).astype(float)
        params_max = np.array(
            (self.Mw_max, self.fc_max, self.t_star_max)).astype(float)
        params_min[np.isnan(params_min)] = 1e-99
        params_max[np.isnan(params_min)] = 1e+99
        tmin = bool(np.all(params >= params_min))
        tmax = bool(np.all(params <= params_max))
        return tmin and tmax

    @property
    def bounds(self):
        """Get bounds for minimize() as sequence of (min, max) pairs."""
        self._bounds = ((self.Mw_min, self.Mw_max),
                        (self.fc_min, self.fc_max),
                        (self.t_star_min, self.t_star_max))
        return self._bounds

    @bounds.setter
    def bounds(self, value):
        """Set bounds from a sequence of three (min, max) pairs."""
        self._bounds = value
        self.Mw_min, self.Mw_max = value[0]
        self.fc_min, self.fc_max = value[1]
        self.t_star_min, self.t_star_max = value[2]

    def get_bounds_curve_fit(self):
        """Get bounds for curve-fit()."""
        bnds = np.array(self.bounds, dtype=float).T
        if np.all(np.isnan(bnds)):
            return None
        bnds[0, np.isnan(bnds[0])] = -1e100
        bnds[1, np.isnan(bnds[1])] = 1e100
        return bnds


# Classes for SourceSpec output

class OrderedAttribDict(OrderedDict):
    """
    An ordered dictionary whose values can be accessed as classattributes.
    """

    def __getattr__(self, attr):
        return self[attr]

    def __setattr__(self, attr, value):
        self[attr] = value


class SpectralParameter(OrderedAttribDict):
    """A spectral parameter measured at one station."""

    def __init__(self, id, name=None, units=None, value=None, uncertainty=None,
                 lower_uncertainty=None, upper_uncertainty=None,
                 confidence_level=None, format=None):
        self._id = id
        self._format = format
        self.name = name
        self.units = units
        self.value = value
        self.uncertainty = uncertainty
        if (lower_uncertainty is not None and
                lower_uncertainty == upper_uncertainty):
            self.uncertainty = lower_uncertainty
            self.lower_uncertainty = self.upper_uncertainty = None
        else:
            self.lower_uncertainty = lower_uncertainty
            self.upper_uncertainty = upper_uncertainty
        self.confidence_level = confidence_level
        self.outlier = False

    def value_uncertainty(self):
        """Return value and uncertainty as 3-element tuple."""
        if self.lower_uncertainty is not None:
            uncertainty = (self.lower_uncertainty, self.upper_uncertainty)
        else:
            uncertainty = (self.uncertainty, self.uncertainty)
        return (self.value, *uncertainty)


class StationParameters(OrderedAttribDict):
    """
    The parameters describing a given station (e.g., its id and location) and
    the spectral parameters measured at that station.

    Spectral parameters are provided as attributes, using SpectralParameter()
    objects.
    """

    def __init__(self, id, instrument_type=None, latitude=None, longitude=None,
                 hypo_dist_in_km=None, epi_dist_in_km=None, azimuth=None):
        self._id = id
        self.instrument_type = instrument_type
        self.latitude = latitude
        self.longitude = longitude
        self.hypo_dist_in_km = hypo_dist_in_km
        self.epi_dist_in_km = epi_dist_in_km
        self.azimuth = azimuth
        self._params = dict()
        self._params_err = dict()
        self._is_outlier = dict()

    def __setattr__(self, attr, value):
        if isinstance(value, SpectralParameter):
            parname = attr
            par = value
            self._params[parname] = par.value
            if par.uncertainty is not None:
                self._params_err[parname] = (par.uncertainty, par.uncertainty)
            else:
                self._params_err[parname] = (
                    par.lower_uncertainty, par.upper_uncertainty
                )
            self._is_outlier[parname] = par.outlier
        self[attr] = value

    def rebuild_dictionaries(self):
        for key, value in self.items():
            if not isinstance(value, SpectralParameter):
                continue
            parname = key
            par = value
            self._params[parname] = par.value
            if par.uncertainty is not None:
                self._params_err[parname] = (par.uncertainty, par.uncertainty)
            elif par.lower_uncertainty is not None:
                self._params_err[parname] = (
                    par.lower_uncertainty, par.upper_uncertainty
                )
            else:
                self._params_err[parname] = (np.nan, np.nan)
            self._is_outlier[parname] = par.outlier


class SummaryStatistics(OrderedAttribDict):
    """
    A summary statistics (e.g., mean, weighted_mean, percentile), along with
    its uncertainty.
    """

    def __init__(self, type, value=None, uncertainty=None,
                 lower_uncertainty=None, upper_uncertainty=None,
                 confidence_level=None, lower_percentage=None,
                 mid_percentage=None, upper_percentage=None,
                 nobs=None, message=None,
                 format=None):
        # type of statistics: e.g., mean, median
        self._type = type
        self.value = value
        self.uncertainty = uncertainty
        if (lower_uncertainty is not None and
                lower_uncertainty == upper_uncertainty):
            self.uncertainty = lower_uncertainty
            self.lower_uncertainty = self.upper_uncertainty = None
        else:
            self.lower_uncertainty = lower_uncertainty
            self.upper_uncertainty = upper_uncertainty
        self.confidence_level = confidence_level
        self.lower_percentage = lower_percentage
        self.mid_percentage = mid_percentage
        self.upper_percentage = upper_percentage
        self.nobs = nobs
        self.message = message
        self._format = format

    def compact_uncertainty(self):
        """Return uncertainty in a compact form."""
        if self.lower_uncertainty is not None:
            return (self.lower_uncertainty, self.upper_uncertainty)
        else:
            return (self.uncertainty, self.uncertainty)


class SummarySpectralParameter(OrderedAttribDict):
    """
    A summary spectral parameter comprising one ore more summary statistics.
    """

    def __init__(self, id, name=None, units=None, format=None):
        self._id = id
        self.name = name
        self.units = units
        # number formatting string
        self._format = format

    def __setattr__(self, attr, value):
        if isinstance(value, SummaryStatistics):
            value._format = self._format
        self[attr] = value


class SourceSpecOutput(OrderedAttribDict):
    """The output of SourceSpec."""

    def __init__(self):
        self.run_info = OrderedAttribDict()
        self.event_info = OrderedAttribDict()
        self.summary_spectral_parameters = OrderedAttribDict()
        self.station_parameters = OrderedAttribDict()

    def value_array(self, key, filter_outliers=False):
        vals = np.array([
            x._params.get(key, np.nan)
            for x in self.station_parameters.values()
        ])
        if filter_outliers:
            outliers = self.outlier_array(key)
            vals = vals[~outliers]
        return vals

    def error_array(self, key, filter_outliers=False):
        errs = np.array([
            x._params_err.get(key, np.nan)
            for x in self.station_parameters.values()
        ])
        if filter_outliers:
            outliers = self.outlier_array(key)
            errs = errs[~outliers]
        return errs

    def outlier_array(self, key):
        outliers = np.array([
            # if we cannot find the given key, we assume outlier=True
            x._is_outlier.get(key, True)
            for x in self.station_parameters.values()
        ])
        return outliers

    def find_outliers(self, key, n):
        """
        Find outliers using the IQR method.

        .. code-block::

                Q1-n*IQR   Q1   median  Q3    Q3+n*IQR
                            |-----:-----|
            o      |--------|     :     |--------|    o  o
                            |-----:-----|
            outlier         <----------->            outliers
                                 IQR

        If ``n`` is ``None``, then the above check is skipped.
        ``Nan`` and ``inf`` values are also marked as outliers.
        """
        values = self.value_array(key)
        station_ids = np.array([id for id in self.station_parameters.keys()])
        naninf = np.logical_or(np.isnan(values), np.isinf(values))
        _values = values[~naninf]
        if n is not None and len(_values) > 0:
            Q1, _, Q3 = np.percentile(_values, [25, 50, 75])
            IQR = Q3-Q1
            outliers = np.logical_or(values < Q1 - n*IQR, values > Q3 + n*IQR)
            outliers = np.logical_or(outliers, naninf)
        else:
            outliers = naninf
        for stat_id, outl in zip(station_ids, outliers):
            stat_par = self.station_parameters[stat_id]
            stat_par[key].outlier = outl
            stat_par.rebuild_dictionaries()

    def mean_values(self):
        """Return a dictionary of mean values."""
        return {
            parname: par.mean.value
            for parname, par in self.summary_spectral_parameters.items()
            if isinstance(par, SummarySpectralParameter)
            and 'mean' in par
        }

    def mean_uncertainties(self):
        """Return a dictionary of mean uncertainties."""
        return {
            parname: par.mean.compact_uncertainty()
            for parname, par in self.summary_spectral_parameters.items()
            if isinstance(par, SummarySpectralParameter)
            and 'mean' in par
        }

    def weighted_mean_values(self):
        """Return a dictionary of weighted mean values."""
        return {
            parname: par.weighted_mean.value
            for parname, par in self.summary_spectral_parameters.items()
            if isinstance(par, SummarySpectralParameter)
            and 'weighted_mean' in par
        }

    def weighted_mean_uncertainties(self):
        """Return a dictionary of weighted mean uncertainties."""
        return {
            parname: par.weighted_mean.compact_uncertainty()
            for parname, par in self.summary_spectral_parameters.items()
            if isinstance(par, SummarySpectralParameter)
            and 'weighted_mean' in par
        }

    def percentiles_values(self):
        """Return a dictionary of percentile values."""
        return {
            parname: par.percentiles.value
            for parname, par in self.summary_spectral_parameters.items()
            if isinstance(par, SummarySpectralParameter)
            and 'percentiles' in par
        }

    def percentiles_uncertainties(self):
        """Return a dictionary of percentile uncertainties."""
        return {
            parname: par.percentiles.compact_uncertainty()
            for parname, par in self.summary_spectral_parameters.items()
            if isinstance(par, SummarySpectralParameter)
            and 'percentiles' in par
        }

    def reference_values(self):
        """Return a dictionary of reference values."""
        try:
            ref_stat = self.summary_spectral_parameters.reference_statistics
        except KeyError:
            raise ValueError('No reference statistics defined')
        if ref_stat == 'mean':
            return self.mean_values()
        elif ref_stat == 'weighted_mean':
            return self.weighted_mean_values()
        elif ref_stat == 'percentiles':
            return self.percentiles_values()
        else:
            msg = 'Invalid reference statistics: {}'.format(ref_stat)
            raise ValueError(msg)

    def reference_uncertainties(self):
        """Return a dictionary of reference uncertainties."""
        try:
            ref_stat = self.summary_spectral_parameters.reference_statistics
        except KeyError:
            raise ValueError('No reference statistics defined')
        if ref_stat == 'mean':
            return self.mean_uncertainties()
        elif ref_stat == 'weighted_mean':
            return self.weighted_mean_uncertainties()
        elif ref_stat == 'percentiles':
            return self.percentiles_uncertainties()
        else:
            msg = 'Invalid reference statistics: {}'.format(ref_stat)
            raise ValueError(msg)

    def reference_summary_parameters(self):
        """
        Return a dictionary of reference summary parameters,
        each being a SummaryStatistics() object.
        """
        try:
            ref_stat = self.summary_spectral_parameters.reference_statistics
        except KeyError:
            raise ValueError('No reference statistics defined')
        return {
            parname: par[ref_stat]
            for parname, par in self.summary_spectral_parameters.items()
            if isinstance(par, SummarySpectralParameter)
            and ref_stat in par
        }


sspec_out_comments = {
    'begin': 'SourceSpec output in YAML format',
    'run_info': 'Information on the SourceSpec run',
    'event_info': 'Information on the event',
    'summary_spectral_parameters':
        'Summary spectral parameters, computed using different statistics',
    'station_parameters':
        'Parameters describing each station and spectral measurements\n'
        'performed at that station'
}
