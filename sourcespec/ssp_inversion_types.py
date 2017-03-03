# -*- coding: utf8 -*-
"""
Classes for spectral inversion routines.

:copyright:
    2017 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement, Version 2.1
    (http://www.cecill.info/index.en.html)
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import logging
import numpy as np


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
        self.fc_min, self.fc_max = self._check_minmax(config.fc_min_max)
        if config.Qo_min_max is None:
            self.t_star_min, self.t_star_max =\
                self._check_minmax(config.t_star_min_max)
        else:
            self.t_star_min, self.t_star_max = self._Qo_to_t_star()
        self._fix_initial_values()

    def __str__(self):
        """String representation."""
        s = 'Mw: {}, {}; '.format(
            *[round(x, 4) if x is not None else x for x in self.bounds[0]])
        s += 'fc: {}, {}; '.format(
            *[round(x, 4) if x is not None else x for x in self.bounds[1]])
        s += 't_star: {}, {}'.format(
            *[round(x, 4) if x is not None else x for x in self.bounds[2]])
        return s

    def _check_minmax(self, minmax):
        if minmax is None:
            return (None, None)
        else:
            return minmax

    def _Qo_to_t_star(self):
        t_star_max, t_star_min =\
            self.hd/(self.config.hypo.vs*np.array(self.config.Qo_min_max))
        return self._nan_to_none((t_star_min, t_star_max))

    def _fix_initial_values_Mw(self):
        if self.ini_values.Mw_0 is not None:
            return
        if None in self.bounds[0]:
            return
        if self.Mw_min < self.ini_values.Mw_0 < self.Mw_max:
            return
        Mw_0 = (self.Mw_max + self.Mw_min) / 2.
        logging.warning('%s %s: initial Mw value: %s outside '
                        'bounds. Using bound average: %s' %
                        (self.spec.id, self.spec.stats.instrtype,
                         self.ini_values.Mw_0, round(Mw_0, 4)))
        self.ini_values.Mw_0 = Mw_0

    def _fix_initial_values_fc(self):
        if self.ini_values.fc_0 is not None:
            return
        if None in self.bounds[1]:
            return
        if self.fc_min < self.ini_values.fc_0 < self.fc_max:
            return
        fc_0 = (self.fc_max + self.fc_min) / 2.
        logging.warning('%s %s: initial fc value: %s outside '
                        'bounds. Using bound average: %s' %
                        (self.spec.id, self.spec.stats.instrtype,
                         self.ini_values.fc_0, round(fc_0, 4)))
        self.ini_values.fc_0 = fc_0

    def _fix_initial_values_t_star(self):
        if self.ini_values.t_star_0 is not None:
            return
        if None in self.bounds[2]:
            return
        if self.t_star_min < self.ini_values.t_star_0 < self.t_star_max:
            return
        t_star_0 = (self.t_star_max + self.t_star_min) / 2.
        logging.warning('%s %s: initial t_star value: %s outside '
                        'bounds. Using bound average: %s' %
                        (self.spec.id, self.spec.stats.instrtype,
                         self.ini_values.t_star_0, round(t_star_0, 4)))
        self.ini_values.t_star_0 = t_star_0

    def _fix_initial_values(self):
        self._fix_initial_values_Mw()
        self._fix_initial_values_fc()
        self._fix_initial_values_t_star()

    def __call__(self, **kwargs):
        """Interface for basin-hopping."""
        params = kwargs['x_new']
        tmin = bool(np.all(params >= self.params_min))
        tmax = bool(np.all(params <= self.params_max))
        return tmin and tmax

    def _get_bounds(self):
        """Get bounds for minimize()."""
        self._bounds = ((self.Mw_min, self.Mw_max),
                        (self.fc_min, self.fc_max),
                        (self.t_star_min, self.t_star_max))
        return self._bounds

    def _set_bounds(self, value):
        self._bounds = value
        self.Mw_min, self.Mw_max = value[0]
        self.fc_min, self.fc_max = value[1]
        self.t_star_min, self.t_star_max = value[2]

    bounds = property(_get_bounds, _set_bounds)

    def get_bounds_curve_fit(self):
        """Get bounds for curve-fit()."""
        bnds = np.array(self.bounds, dtype=float).T
        if np.all(np.isnan(bnds)):
            return None
        bnds[0, np.isnan(bnds[0])] = -1e100
        bnds[1, np.isnan(bnds[1])] = 1e100
        return bnds
