
# -*- coding: utf8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
A class for sampling a parameter space over a grid.

Sampling can be perfomred by several approaches.
The class provides optimal solutions, uncertainties and plotting methods.

:copyright:
    2022 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
import os
import numpy as np
import logging
logger = logging.getLogger(__name__.split('.')[-1])


class GridSampling():
    """
    A class for sampling a parameter space over a grid.

    Sampling can be perfomred by several approaches.
    The class provides optimal solutions, uncertainties and plotting methods.
    """

    def __init__(self, misfit_func, bounds, nsteps):
        """
        Init grid sampling.

        bounds : sequence of (min, max) pairs for each dimension.
        nsteps : number of grid steps for each dimension
        """
        self.misfit_func = misfit_func
        self.bounds = bounds
        self.nsteps = nsteps
        self.misfit = None
        self._values = None
        self._min_idx = None

    @property
    def values(self):
        if self._values is not None:
            return self._values
        values = []
        for bds, ns in zip(self.bounds, self.nsteps):
            if None in bds:
                msg = 'All parameters must be bounded for grid sampling'
                raise RuntimeError(msg)
            values.append(np.linspace(*bds, ns))
        self._values = np.meshgrid(*values, indexing='ij')
        return self._values

    @property
    def min_idx(self):
        if self.misfit is None:
            return None
        if self._min_idx is None:
            return np.unravel_index(np.argmin(self.misfit), self.misfit.shape)
        else:
            return self._min_idx

    @property
    def params_opt(self):
        if self.misfit is None:
            return None
        return np.array([v[self.min_idx] for v in self.values])

    def grid_search(self):
        """Sample the misfit function by simple grid search."""
        # small helper function to transform args into a tuple
        def mf(*args):
            return self.misfit_func(args)
        mf = np.vectorize(mf)
        self.misfit = mf(*self.values)

    def plot_conditional(self, config, params_name, params_unit, label):
        """Plot conditional misfit for each parameter."""
        import matplotlib.pyplot as plt
        ndim = self.misfit.ndim
        fig, ax = plt.subplots(ndim, 1, figsize=(9, 9), dpi=300)
        for dim in range(ndim):
            # `dim` is the dimension to keep
            # we move `dim` to the last axis
            mm = np.moveaxis(self.misfit, dim, -1)
            # we fix ndim-1 coordinates of the minimum
            idx = tuple(v for n, v in enumerate(self.min_idx) if n != dim)
            # we extract from mm a 1-d profile along dim,
            # by fixing all the other dimensions (conditional misfit)
            mm = mm[idx]
            # same thing for values: we extract a 1d array of values along dim
            v = np.moveaxis(self.values[dim], dim, -1)
            v = v[0, 0]
            ax[dim].plot(v, mm)
            ax[dim].axvline(self.params_opt[dim], color='red')
            xlabel = params_name[dim]
            unit = params_unit[dim]
            if unit:
                xlabel += ' ({})'.format(unit)
            ax[dim].set_xlabel(xlabel)
            ax[dim].set_ylabel('misfit')
        ax[0].set_title(label)
        plt.tight_layout()
        outdir = os.path.join(config.options.outdir, 'misfit')
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        evid = config.hypo.evid
        figfile_base = os.path.join(outdir, evid)
        figfile_base += '.misfit_{}.'.format(label.replace(' ', '_'))
        fmt = config.PLOT_SAVE_FORMAT
        if fmt == 'pdf_multipage':
            fmt = 'pdf'
        figfile = figfile_base + fmt
        if config.PLOT_SHOW:
            plt.show()
        if config.PLOT_SAVE:
            fig.savefig(figfile, bbox_inches='tight')
            logger.info(
                '{}: conditional misfit plot saved to: {}'.format(
                    label, figfile))
