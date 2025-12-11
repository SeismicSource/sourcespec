# -*- coding: utf8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Classes for spectral inversion routines.

:copyright:
    2017-2025 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
import logging
from collections import OrderedDict
import numpy as np
logger = logging.getLogger(__name__.rsplit('.', maxsplit=1)[-1])


class InitialValues():
    """Initial values for spectral inversion."""

    def __init__(self, Mw_0=None, fc_0=None, t_star_0=None):
        self.Mw_0 = Mw_0
        self.fc_0 = fc_0
        self.t_star_0 = t_star_0

    def __str__(self):
        """String representation."""
        t_star_0 = (
            'None' if self.t_star_0 is None else round(self.t_star_0, 4))
        return (
            f'Mw_0: {round(self.Mw_0, 4)}; '
            f'fc_0: {round(self.fc_0, 4)}; '
            f't_star_0: {t_star_0}'
        )

    def get_params0(self):
        """Get initial values as a tuple."""
        if self.t_star_0 is None:
            return (self.Mw_0, self.fc_0)
        return (self.Mw_0, self.fc_0, self.t_star_0)


class Bounds():
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
        # pylint: disable=consider-using-f-string
        s = 'Mw: {}, {}; '.format(
            *[round(x, 4) if x is not None else x for x in self.bounds[0]])
        s += 'fc: {}, {}; '.format(
            *[round(x, 4) if x is not None else x for x in self.bounds[1]])
        try:
            s += 't_star: {}, {}'.format(
                *[round(x, 4) if x is not None else x for x in self.bounds[2]])
        except IndexError:
            s += 't_star: None, None'
        return s

    def _set_fc_min_max(self, config):
        fc_0 = self.ini_values.fc_0
        if config.fc_min_max is None:
            # If no bound is given, set it to fc_0 +/- a decade
            scale = 10.  # a decade
            self.fc_min = fc_0 / scale
            self.fc_max = fc_0 * scale
        else:
            self.fc_min, self.fc_max = config.fc_min_max
        if self.fc_min > fc_0:
            logger.warning(
                f'{self.spec.id} {self.spec.stats.instrtype}: '
                f'fc_min ({self.fc_min}) larger than '
                f'fc_0 ({round(fc_0, 4)}). '
                'Using fc_0 instead.'
            )
            self.fc_min = fc_0
        if self.fc_max < fc_0:
            logger.warning(
                f'{self.spec.id} {self.spec.stats.instrtype}: '
                f'fc_max ({self.fc_max}) smaller than '
                f'fc_0 ({round(fc_0, 4)}). '
                'Using fc_0 instead.'
            )
            self.fc_max = fc_0

    def _check_minmax(self, minmax):
        return (None, None) if minmax is None else minmax

    def _Qo_to_t_star(self):
        phase = self.config.wave_type[0]
        travel_time = self.spec.stats.travel_times[phase]
        t_star_bounds = travel_time / np.array(self.config.Qo_min_max)
        return sorted(t_star_bounds)

    def _fix_initial_values_t_star(self):
        if self.ini_values.t_star_0 is None:
            return
        if None in self.bounds[2]:
            return
        if self.t_star_min < self.ini_values.t_star_0 < self.t_star_max:
            return
        t_star_0 = (self.t_star_max + self.t_star_min) / 2.
        logger.warning(
            f'{self.spec.id} {self.spec.stats.instrtype}: '
            f'initial t_star value ({self.ini_values.t_star_0}) '
            f'outside bounds. Using bound average ({round(t_star_0, 4)})'
        )
        self.ini_values.t_star_0 = t_star_0

    def __call__(self, **kwargs):
        """Interface for basin-hopping."""
        params = kwargs['x_new']
        params_min = np.array([b[0] for b in self.bounds]).astype(float)
        params_max = np.array([b[1] for b in self.bounds]).astype(float)
        params_min[np.isnan(params_min)] = 1e-99
        params_max[np.isnan(params_min)] = 1e+99
        tmin = bool(np.all(params >= params_min))
        tmax = bool(np.all(params <= params_max))
        return tmin and tmax

    @property
    def bounds(self):
        """Get bounds for minimize() as sequence of (min, max) pairs."""
        if self.t_star_min is None and self.t_star_max is None:
            self._bounds = ((self.Mw_min, self.Mw_max),
                            (self.fc_min, self.fc_max))
        else:
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

    def __init__(self, param_id, name=None, units=None, value=None,
                 uncertainty=None,
                 lower_uncertainty=None, upper_uncertainty=None,
                 confidence_level=None, format_spec=None):
        self.param_id = param_id
        self._format_spec = format_spec
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

    def compact_uncertainty(self):
        """Return uncertainty in a compact form."""
        if self.lower_uncertainty is not None:
            return (self.lower_uncertainty, self.upper_uncertainty)
        if self.uncertainty is not None:
            return (self.uncertainty, self.uncertainty)
        return (np.nan, np.nan)


class StationParameters(OrderedAttribDict):
    """
    The parameters describing a given station (e.g., its id and location) and
    the spectral parameters measured at that station.

    Spectral parameters are provided as attributes, using SpectralParameter()
    objects.
    """

    def __init__(self, station_id, instrument_type=None,
                 latitude=None, longitude=None,
                 hypo_dist_in_km=None, epi_dist_in_km=None, azimuth=None,
                 spectral_snratio_mean=None, spectral_snratio_max=None,
                 rmsn=None, quality_of_fit=None,
                 ignored=False, ignored_reason=None):
        self.station_id = station_id
        self.instrument_type = instrument_type
        self.latitude = latitude
        self.longitude = longitude
        self.hypo_dist_in_km = hypo_dist_in_km
        self.epi_dist_in_km = epi_dist_in_km
        self.azimuth = azimuth
        self.spectral_snratio_mean = spectral_snratio_mean
        self.spectral_snratio_max = spectral_snratio_max
        self.rmsn = rmsn
        self.quality_of_fit = quality_of_fit
        self.ignored = ignored
        self.ignored_reason = ignored_reason
        # The following parameters are expected to be of type SpectralParameter
        self.Mw = None
        self.fc = None
        self.t_star = None
        self.Mo = None
        self.radius = None
        self.ssd = None
        self.Qo = None
        self.Er = None
        self.sigma_a = None
        self.Ml = None

    def get_spectral_parameters(self):
        """Return a dictionary of spectral parameters."""
        return {
            key: value
            for key, value in self.items()
            if isinstance(value, SpectralParameter)
        }


class SummaryStatistics(OrderedAttribDict):
    """
    A summary statistics (e.g., mean, weighted_mean, percentile), along with
    its uncertainty.
    """

    def __init__(self, stat_type, value=None, uncertainty=None,
                 lower_uncertainty=None, upper_uncertainty=None,
                 confidence_level=None, lower_percentage=None,
                 mid_percentage=None, upper_percentage=None,
                 nobs=None, message=None,
                 format_spec=None):
        # type of statistics: e.g., mean, median
        self._stat_type = stat_type
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
        self._format_spec = format_spec

    def compact_uncertainty(self):
        """Return uncertainty in a compact form."""
        if self.lower_uncertainty is not None:
            return (self.lower_uncertainty, self.upper_uncertainty)
        if self.uncertainty is not None:
            return (self.uncertainty, self.uncertainty)
        return (np.nan, np.nan)


class SummarySpectralParameter(OrderedAttribDict):
    """
    A summary spectral parameter comprising one ore more summary statistics.
    """

    def __init__(self, param_id, name=None, units=None, format_spec=None):
        self.param_id = param_id
        self.name = name
        self.units = units
        # number formatting string
        self._format_spec = format_spec

    def __setattr__(self, attr, value):
        if isinstance(value, SummaryStatistics):
            value._format_spec = self._format_spec
        self[attr] = value


class SourceSpecOutput(OrderedAttribDict):
    """The output of SourceSpec."""

    def __init__(self):
        self.run_info = OrderedAttribDict()
        self.event_info = OrderedAttribDict()
        self.inversion_info = OrderedAttribDict()
        self.quality_info = OrderedAttribDict()
        self.summary_spectral_parameters = OrderedAttribDict()
        self.station_parameters = OrderedAttribDict()
        # comments for each section
        # do not remove the underscore from the attribute name!
        self._comments = {
            'begin': 'SourceSpec output in YAML format',
            'run_info': 'Information on the SourceSpec run',
            'event_info': 'Information on the event',
            'inversion_info': 'Information on the inversion procedure',
            'quality_info': 'Quality information derived from the inversion',
            'summary_spectral_parameters':
                'Summary spectral parameters, computed using different '
                'statistics',
            'station_parameters':
                'Parameters describing each station and spectral '
                'measurements\nperformed at that station'
        }

    def value_array(self, key, filter_outliers=False):
        """Return an array of values for the given key."""
        vals = np.array([
            stat_par[key].value
            if key in stat_par and stat_par[key] is not None
            else np.nan
            for stat_par in self.station_parameters.values()
        ])
        if filter_outliers:
            outliers = self.outlier_array(key)
            vals = vals[~outliers]
        return vals

    def error_array(self, key, filter_outliers=False):
        """Return an array of errors (two columns) for the given key."""
        errs = np.array([
            stat_par[key].compact_uncertainty()
            if key in stat_par and stat_par[key] is not None
            else (np.nan, np.nan)
            for stat_par in self.station_parameters.values()
        ])
        if filter_outliers:
            outliers = self.outlier_array(key)
            errs = errs[~outliers]
        return errs

    def outlier_array(self, key):
        """Return an array of outliers for the given key."""
        return np.array(
            [
                stat_par[key].outlier
                if key in stat_par and stat_par[key] is not None
                # if we cannot find the given key, we assume outlier=True
                else True
                for stat_par in self.station_parameters.values()
            ]
        )

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
        naninf = np.logical_or(np.isnan(values), np.isinf(values))
        _values = values[~naninf]
        # cases for which outliers cannot be computed
        if (
            n is None or
            len(_values) == 0 or
            # _values are all the same within 0.01 %
            np.ptp(_values) < 0.0001 * np.mean(_values)
        ):
            outliers = naninf
        else:
            Q1, _, Q3 = np.percentile(_values, [25, 50, 75])
            IQR = Q3 - Q1
            outliers = np.logical_or(
                values < Q1 - n * IQR, values > Q3 + n * IQR)
            outliers = np.logical_or(outliers, naninf)
        station_ids = np.array(list(self.station_parameters.keys()))
        for stat_id, outl in zip(station_ids, outliers):
            stat_par = self.station_parameters[stat_id]
            if key in stat_par and stat_par[key] is not None:
                stat_par[key].outlier = outl

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

    def mean_nobs(self):
        """
        Return a dictionary of number of observations used for computing mean.
        """
        return {
            parname: par.mean.nobs
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

    def weighted_mean_nobs(self):
        """
        Return a dictionary of number of observations used for computing
        weighted mean.
        """
        return {
            parname: par.weighted_mean.nobs
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

    def percentiles_nobs(self):
        """
        Return a dictionary of number of observations used for computing
        percentiles.
        """
        return {
            parname: par.percentiles.nobs
            for parname, par in self.summary_spectral_parameters.items()
            if isinstance(par, SummarySpectralParameter)
            and 'percentiles' in par
        }

    def reference_values(self):
        """Return a dictionary of reference values."""
        try:
            ref_stat = self.summary_spectral_parameters.reference_statistics
        except KeyError as e:
            raise ValueError('No reference statistics defined') from e
        if ref_stat == 'mean':
            return self.mean_values()
        if ref_stat == 'percentiles':
            return self.percentiles_values()
        if ref_stat == 'weighted_mean':
            return self.weighted_mean_values()
        msg = f'Invalid reference statistics: {ref_stat}'
        raise ValueError(msg)

    def reference_uncertainties(self):
        """Return a dictionary of reference uncertainties."""
        try:
            ref_stat = self.summary_spectral_parameters.reference_statistics
        except KeyError as e:
            raise ValueError('No reference statistics defined') from e
        if ref_stat == 'mean':
            return self.mean_uncertainties()
        if ref_stat == 'percentiles':
            return self.percentiles_uncertainties()
        if ref_stat == 'weighted_mean':
            return self.weighted_mean_uncertainties()
        msg = f'Invalid reference statistics: {ref_stat}'
        raise ValueError(msg)

    def reference_summary_parameters(self):
        """
        Return a dictionary of reference summary parameters,
        each being a SummaryStatistics() object.
        """
        try:
            ref_stat = self.summary_spectral_parameters.reference_statistics
        except KeyError as e:
            raise ValueError('No reference statistics defined') from e
        return {
            parname: par[ref_stat]
            for parname, par in self.summary_spectral_parameters.items()
            if isinstance(par, SummarySpectralParameter)
            and ref_stat in par
        }
