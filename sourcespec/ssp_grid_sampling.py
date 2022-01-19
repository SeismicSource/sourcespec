
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
import matplotlib.pyplot as plt
import logging
logger = logging.getLogger(__name__.split('.')[-1])


class GridSampling():
    """
    A class for sampling a parameter space over a grid.

    Sampling can be perfomred by several approaches.
    The class provides optimal solutions, uncertainties and plotting methods.
    """

    def __init__(self, misfit_func, bounds, nsteps, sampling_mode,
                 params_name, params_unit):
        """
        Init grid sampling.

        bounds : sequence of (min, max) pairs for each dimension.
        nsteps : number of grid steps for each dimension.
        sampling_mode : sequence of 'lin' or 'log' for each diemesion.
        params_name : sequence of parameter names (str).
        params_name : sequence of parameter units (str).
        """
        self.misfit_func = misfit_func
        self.bounds = bounds
        self.nsteps = nsteps
        self.sampling_mode = sampling_mode
        self.params_name = params_name
        self.params_unit = params_unit
        self.misfit = None
        self._cut_misfit = None
        self._values = None
        self._min_idx = None
        self.truebounds = []
        for bds, ns, mode in zip(self.bounds, self.nsteps, self.sampling_mode):
            if None in bds:
                msg = 'All parameters must be bounded for grid sampling'
                raise RuntimeError(msg)
            if mode == 'log':
                if bds[0] == 0:
                    bds = (bds[1]/ns, bds[1])
                bds = tuple(np.log10(bds))
            self.truebounds.append(bds)

    @property
    def values(self):
        if self._values is not None:
            return self._values
        values = []
        for bds, ns, mode in zip(
                self.truebounds, self.nsteps, self.sampling_mode):
            if mode == 'log':
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

    def values_1d(self, dim):
        """Extract a 1D array of parameter values along one dimension."""
        # same thing for values: we extract a 1d array of values along dim
        v = np.moveaxis(self.values[dim], dim, -1)
        return v[0, 0]

    @property
    def conditional_misfit(self):
        """
        Compute conditional misfit along each dimension.

        Conditional misfit is computed by fixing the other parameters to
        their optimal value.
        """
        if self.misfit is None:
            return None
        ndim = self.misfit.ndim
        cond_misfit = []
        for dim in range(ndim):
            # `dim` is the dimension to keep
            # we move `dim` to the last axis
            mm = np.moveaxis(self.misfit, dim, -1)
            # we fix ndim-1 coordinates of the minimum
            idx = tuple(v for n, v in enumerate(self.min_idx) if n != dim)
            # we extract from mm a 1-d profile along dim,
            # by fixing all the other dimensions (conditional misfit)
            mm = mm[idx]
            cond_misfit.append(mm)
        return tuple(cond_misfit)

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

    def plot_conditional_misfit(self, config, label):
        """Plot conditional misfit for each parameter."""
        ndim = self.misfit.ndim
        fig, ax = plt.subplots(ndim, 1, figsize=(9, 9), dpi=300)
        for dim, mm in enumerate(self.conditional_misfit):
            v = self.values_1d(dim)
            ax[dim].plot(v, mm)
            popt = self.params_opt[dim]
            ax[dim].axvline(popt, color='red')
            text = '{:.4f}  '.format(popt)
            ax[dim].text(
                popt, 0.9, text, color='red', ha='right',
                transform=ax[dim].get_xaxis_transform())
            if self.sampling_mode[dim] == 'log':
                ax[dim].set_xscale('log')
            xlabel = self.params_name[dim]
            unit = self.params_unit[dim]
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
            if not config.PLOT_SHOW:
                plt.close(fig)
            logger.info(
                '{}: conditional misfit plot saved to: {}'.format(
                    label, figfile))

    def plot_misfit_2d(self, config, plot_par_idx, label):
        """Plot a 2D misfit map."""
        # Find the index to extract
        idx = tuple(
            v for n, v in enumerate(self.min_idx) if n not in plot_par_idx)
        if plot_par_idx[0] < plot_par_idx[1]:
            # Move axes to keep at the end
            end_idx = (-2, -1)
            mm = np.moveaxis(self.misfit, plot_par_idx, end_idx)
            # extract a 2D misfit map
            mm = mm[idx].T
        else:
            # Move axes to keep at the end
            end_idx = (-1, -2)
            mm = np.moveaxis(self.misfit, plot_par_idx, end_idx)
            # extract a 2D misfit map
            mm = mm[idx]
        # set extent for imshow
        bds = np.take(np.array(self.truebounds), plot_par_idx, axis=0)
        extent = bds.flatten()
        # take log of parameters, if necessary
        cond = np.array(self.sampling_mode) == 'log'
        params_opt = np.array(self.params_opt)
        params_opt[cond] = np.log10(params_opt[cond])
        # extract parameter info only for the two selected parameters
        params_opt = np.take(params_opt, plot_par_idx)
        params_name = np.take(self.params_name, plot_par_idx)
        params_unit = np.take(self.params_unit, plot_par_idx)
        sampling_mode = np.take(self.sampling_mode, plot_par_idx)
        fig, ax = plt.subplots(1, 1, figsize=(9, 8), dpi=300)
        # imshow: saturate scale at 2 times the minimum
        mmap = ax.imshow(
            mm, vmax=2*mm.min(), origin='lower', cmap='viridis',
            extent=extent, aspect='auto'
        )
        ax.set_xlim(extent[0], extent[1])
        ax.set_ylim(extent[2], extent[3])
        ax.scatter(*params_opt, facecolors='none', edgecolors='w')
        ax.set_title(label)
        xlabel = params_name[0]
        unit = params_unit[0]
        if unit:
            xlabel += ' ({})'.format(unit)
        if sampling_mode[0] == 'log':
            # the grid is plotted by imshow() in linear scale,
            # so we create a fake xaxis with logscale
            ax.xaxis.set_visible(False)
            ax2 = ax.twiny()
            ax2.xaxis.set_label_position('bottom')
            ax2.xaxis.set_ticks_position('bottom')
            ax2.set_xlim(10**extent[0], 10**extent[1])
            ax2.set_xscale('log')
            ax2.set_xlabel(xlabel)
        else:
            ax.set_xlabel(xlabel)
        ylabel = params_name[1]
        unit = params_unit[1]
        if unit:
            ylabel += ' ({})'.format(unit)
        if sampling_mode[1] == 'log':
            # the grid is plotted by imshow() in linear scale,
            # so we create a fake yaxis with logscale
            ax.yaxis.set_visible(False)
            ax2 = ax.twinx()
            ax2.yaxis.set_label_position('left')
            ax2.yaxis.set_ticks_position('left')
            ax2.set_ylim(10**extent[2], 10**extent[3])
            ax2.set_yscale('log')
            ax2.set_ylabel(ylabel)
        else:
            ax.set_ylabel(ylabel)
        cbar = plt.colorbar(mmap, ax=ax, extend='max')
        cbar.set_label('misfit')
        outdir = os.path.join(config.options.outdir, 'misfit')
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        evid = config.hypo.evid
        figfile_base = os.path.join(outdir, evid)
        suffix = '{}-{}'.format(*params_name)
        figfile_base += '.misfit_{}_{}.'.format(
            suffix, label.replace(' ', '_'))
        fmt = config.PLOT_SAVE_FORMAT
        if fmt == 'pdf_multipage':
            fmt = 'pdf'
        figfile = figfile_base + fmt
        if config.PLOT_SHOW:
            plt.show()
        if config.PLOT_SAVE:
            fig.savefig(figfile, bbox_inches='tight')
            if not config.PLOT_SHOW:
                plt.close(fig)
            logger.info(
                '{}: misfit map saved to: {}'.format(
                    label, figfile))
