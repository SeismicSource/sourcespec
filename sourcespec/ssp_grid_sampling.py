
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

    def __init__(self, misfit_func, bounds, nsteps, sampling_mode):
        """
        Init grid sampling.

        bounds : sequence of (min, max) pairs for each dimension.
        nsteps : number of grid steps for each dimension.
        sampling_mode : sequence of 'lin' or 'log' for each diemesion.
        """
        self.misfit_func = misfit_func
        self.bounds = bounds
        self.nsteps = nsteps
        self.sampling_mode = sampling_mode
        self.misfit = None
        self._values = None
        self._min_idx = None

    @property
    def values(self):
        if self._values is not None:
            return self._values
        values = []
        for bds, ns, mode in zip(self.bounds, self.nsteps, self.sampling_mode):
            if None in bds:
                msg = 'All parameters must be bounded for grid sampling'
                raise RuntimeError(msg)
            if mode == 'log':
                if bds[0] == 0:
                    bds = (bds[1]/ns, bds[1])
                bds = np.log10(bds)
                values.append(np.logspace(*bds, ns))
            else:
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

    def plot_conditional_misfit(self, config, params_name, params_unit, label):
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
            popt = self.params_opt[dim]
            ax[dim].axvline(popt, color='red')
            text = '{:.4f}  '.format(popt)
            ax[dim].text(
                popt, 0.9, text, color='red', ha='right',
                transform=ax[dim].get_xaxis_transform())
            if self.sampling_mode[dim] == 'log':
                ax[dim].set_xscale('log')
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
        figfile_base += '.cond_misfit_{}.'.format(label.replace(' ', '_'))
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
