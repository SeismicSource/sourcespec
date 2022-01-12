# -*- coding: utf8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
A simple grid search algorithm for sourcespec.

:copyright:
    2022 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
import numpy as np


def grid_search(misfit_func, bounds, nsteps):
    """
    Search the minimum of ``misfit_func`` over a grid.

    bounds : sequence of (min, max) pairs for each dimension.
    nsteps : number of grid steps for each dimension
    """
    values = []
    for bds, ns in zip(bounds, nsteps):
        if None in bds:
            msg = 'All parameters must be bounded for grid search'
            raise RuntimeError(msg)
        values.append(np.linspace(*bds, ns))
    values = np.meshgrid(*values, indexing='ij')

    # small helper function to transform args into a tuple
    def mf(*args):
        return misfit_func(args)
    mf = np.vectorize(mf)
    misfit = mf(*values)
    min_idx = np.unravel_index(np.argmin(misfit), misfit.shape)
    params_opt = np.array([v[min_idx] for v in values])
    return params_opt
