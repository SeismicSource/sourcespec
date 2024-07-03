
# -*- coding: utf8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
A class for sampling a parameter space over a grid.

Sampling can be performed by several approaches.
The class provides optimal solutions, uncertainties and plotting methods.

:copyright:
    2022-2024 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
import os
import warnings
import logging
import numpy as np
from scipy.signal import peak_widths as find_peak_widths
# pylint: disable=no-name-in-module
from scipy.signal._peak_finding_utils import PeakPropertyWarning
import matplotlib.pyplot as plt
from .kdtree import KDTree
from .savefig import savefig
logger = logging.getLogger(__name__.rsplit('.', maxsplit=1)[-1])


def find_peak_width(x, peak_idx, rel_height, negative=False):
    """
    Find width of a single peak at a given relative height.

    rel_height: float parameter between 0 and 1
                0 means the base of the curve and 1 the peak value
                (Note: this is the opposite of scipy.peak_widths)
    """
    if rel_height < 0 or rel_height > 1:
        msg = 'rel_height must be between 0 and 1'
        raise ValueError(msg)
    sign = -1 if negative else 1
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', category=PeakPropertyWarning)
        _, width_height, idx_left, idx_right =\
            find_peak_widths(sign * x, [peak_idx, ], 1 - rel_height)
    idx_left = int(idx_left)
    idx_right = int(idx_right)
    width_height = sign * width_height[0]
    # fall back approach if the previous one fails
    if idx_left == idx_right:
        height = x.max() - x.min()
        if not negative:
            rel_height = 1 - rel_height
        width_height = x.max() - rel_height * height
        # Search for the indexes of the misfit curve points which are
        # closest to width_height, on the left and on the right
        #   Note: This assumes that the misfit function is monotonic.
        #         A safer but less precise approach is the commented one,
        #         based on "iii"
        # iii = np.where(np.isclose(x, width_height, rtol=0.1))
        try:
            # idx_left = np.min(iii[iii < peak_idx])
            x2 = x.copy()
            x2[peak_idx:] = np.inf
            idx_left = np.argmin(np.abs(x2 - width_height))
        except ValueError:
            idx_left = 0
        try:
            # idx_right = np.max(iii[iii > peak_idx])
            x2 = x.copy()
            x2[:peak_idx] = np.inf
            idx_right = np.argmin(np.abs(x2 - width_height))
        except ValueError:
            idx_right = -1
    return width_height, idx_left, idx_right


class GridSampling():
    """
    A class for sampling a parameter space over a grid.

    Sampling can be performed by several approaches.
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
        self._conditional_misfit = None
        self._conditional_peak_widths = None
        self._values = None
        self._min_idx = None
        self.truebounds = []
        for bds, ns, mode in zip(self.bounds, self.nsteps, self.sampling_mode):
            if None in bds:
                msg = 'All parameters must be bounded for grid sampling'
                raise RuntimeError(msg)
            if mode == 'log':
                if bds[0] == 0:
                    bds = (bds[1] / ns, bds[1])
                bds = tuple(np.log10(bds))
            self.truebounds.append(bds)
        self.kdt = None
        self.extent = None

    @property
    def values(self):
        """Return a meshgrid of parameter values."""
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
        """Find the index of the minimum of the misfit function."""
        if self.misfit is None:
            return None
        if self._min_idx is None:
            return np.unravel_index(
                np.nanargmin(self.misfit), self.misfit.shape)
        return self._min_idx

    @property
    def values_1d(self):
        """Extract a 1D array of parameter values along one dimension."""
        # same thing for values: we extract a 1d array of values along dim
        ndim = len(self.values)
        values_1d = []
        for dim in range(ndim):
            v = np.moveaxis(self.values[dim], dim, -1)
            values_1d.append(v[0, 0])
        return tuple(values_1d)

    @property
    def conditional_misfit(self):
        """
        Compute conditional misfit along each dimension.

        Conditional misfit is computed by fixing the other parameters to
        their optimal value.
        """
        if self.misfit is None:
            return None
        if self._conditional_misfit is not None:
            return self._conditional_misfit
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
        self._conditional_misfit = tuple(cond_misfit)
        return self._conditional_misfit

    @property
    def params_opt(self):
        """Return optimal parameters."""
        if self.misfit is None:
            return None
        return np.array([v[self.min_idx] for v in self.values])

    @property
    def params_err(self):
        """Return optimal parameters uncertainties."""
        if self.misfit is None:
            return None
        error = []
        for p, w in zip(self.params_opt, self.conditional_peak_widths):
            err_left = p - w[1]
            err_right = w[2] - p
            error.append((err_left, err_right))
        return tuple(error)

    @property
    def conditional_peak_widths(self):
        """Find width of conditional misfit around its minimum."""
        if self.misfit is None:
            return None
        if self._conditional_peak_widths is not None:
            return self._conditional_peak_widths
        peak_widths = []
        rel_height = np.exp(-0.5)  # height of a gaussian for x=sigma
        for mm, idx, values in zip(
                self.conditional_misfit, self.min_idx, self.values_1d):
            width_height, idx_left, idx_right = find_peak_width(
                mm, idx, rel_height, negative=True)
            peak_widths.append(
                (width_height, values[idx_left], values[idx_right]))
        self._conditional_peak_widths = tuple(peak_widths)
        return self._conditional_peak_widths

    def grid_search(self):
        """Sample the misfit function by simple grid search."""
        # small helper function to transform args into a tuple
        def mf(*args):
            return self.misfit_func(args)
        mf = np.vectorize(mf)
        self.misfit = mf(*self.values)

    def kdtree_search(self):
        """Sample the misfit function using kdtree search."""
        # small helper function to transform misfit to pdf and manage logscale
        def mf(args):
            newargs = []
            for a, mode in zip(args, self.sampling_mode):
                if mode == 'log':
                    a = 10**a
                newargs.append(a)
            return np.exp(-self.misfit_func(newargs))
        extent = sum(self.truebounds, ())
        maxdiv = (20, 2000, 200)
        kdt = KDTree(extent, 2, mf, maxdiv=maxdiv)
        while kdt.ncells <= np.prod(maxdiv):
            oldn = kdt.ncells
            kdt.divide()
            if kdt.ncells == oldn:
                break
        deltas = []
        for bds, ns in zip(self.truebounds, self.nsteps):
            deltas.append((bds[1] - bds[0]) / ns)
        pdf, extent = kdt.get_pdf(deltas)
        self.kdt = kdt
        self.misfit = -np.log(pdf)
        self.nsteps = self.misfit.shape
        self.extent = extent

    def plot_conditional_misfit(self, config, label):
        """Plot conditional misfit for each parameter."""
        # Check config, if we need to plot at all
        if not config.plot_show and not config.plot_save:
            return
        ndim = self.misfit.ndim
        fig, ax = plt.subplots(ndim, 1, figsize=(5, 5), dpi=300)
        for dim, mm in enumerate(self.conditional_misfit):
            v = self.values_1d[dim]
            ax[dim].plot(v, mm)
            popt = self.params_opt[dim]
            w = self.conditional_peak_widths[dim]
            ax[dim].plot((w[1], w[2]), (w[0], w[0]), color='red', ls='dashed')
            ax[dim].axvline(popt, color='red')
            text = f'  {popt:.3f}  '
            # set horizontal alignment based on whether we are
            # on the left or right part of the plot
            ha = 'left' if popt < np.nanmean(v) else 'right'
            ax[dim].text(
                popt, 0.9, text, color='red', ha=ha, va='top',
                transform=ax[dim].get_xaxis_transform())
            if self.sampling_mode[dim] == 'log':
                ax[dim].set_xscale('log')
            xlabel = self.params_name[dim]
            if self.params_unit[dim]:
                xlabel += f' ({self.params_unit[dim]})'
            # add some margin if values are too close
            tolerance = 1e-3
            if np.isclose(v.min(), v.max(), rtol=tolerance):
                vmean = np.mean(v)
                vmin = vmean - tolerance * np.abs(vmean)
                vmax = vmean + tolerance * np.abs(vmean)
                ax[dim].set_xlim(vmin, vmax)
                ax[dim].set_ylim(0, 1)
            ax[dim].set_xlabel(xlabel)
            ax[dim].set_ylabel('misfit')
        ax[0].set_title(label)
        plt.tight_layout()
        figfile_base = self._get_figfile_base(config)
        figfile_base += f".cond_misfit_{label.replace(' ', '_')}."
        fmt = config.plot_save_format
        if fmt == 'pdf_multipage':
            fmt = 'pdf'
        figfile = figfile_base + fmt
        if config.plot_show:
            plt.show()
        if config.plot_save:
            savefig(fig, figfile, fmt, bbox_inches='tight')
            if not config.plot_show:
                plt.close(fig)
            logger.info(
                f'{label}: conditional misfit plot saved to: {figfile}')
            config.figures['misfit_1d'].append(figfile)

    def plot_misfit_2d(self, config, plot_par_idx, label):
        """Plot a 2D conditional misfit map."""
        # Check config, if we need to plot at all
        if not config.plot_show and not config.plot_save:
            return
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
        params_opt_all = np.array(self.params_opt)
        params_opt_all[cond] = np.log10(params_opt_all[cond])
        # extract parameter info only for the two selected parameters
        params_opt = np.take(params_opt_all, plot_par_idx)
        params_name = np.take(self.params_name, plot_par_idx)
        params_unit = np.take(self.params_unit, plot_par_idx)
        sampling_mode = np.take(self.sampling_mode, plot_par_idx)
        fig, ax = plt.subplots(1, 1, figsize=(5, 4), dpi=300)
        # imshow: saturate scale at 2 times the minimum
        mmap = ax.imshow(
            mm, vmax=2 * mm.min(), origin='lower', cmap='viridis',
            extent=extent, aspect='auto'
        )
        xmin, xmax, ymin, ymax = extent
        # add some margin if values are too close
        tolerance = 1e-3
        if np.isclose(xmin, xmax, rtol=tolerance):
            xmean = np.mean((xmin, xmax))
            xmin = xmean - tolerance * np.abs(xmean)
            xmax = xmean + tolerance * np.abs(xmean)
        if np.isclose(ymin, ymax, rtol=tolerance):
            ymean = np.mean((ymin, ymax))
            ymin = ymean - tolerance * np.abs(ymean)
            ymax = ymean + tolerance * np.abs(ymean)
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
        if self.kdt is not None:
            self._plot_kdtree_structure(ax, plot_par_idx, params_opt_all)
        ax.scatter(*params_opt, facecolors='none', edgecolors='w')
        ax.set_title(label)
        xlabel = params_name[0]
        if params_unit[0]:
            xlabel += f' ({params_unit[0]})'
        if sampling_mode[0] == 'log':
            self._set_xlogscale(ax, extent, xlabel)
        else:
            ax.set_xlabel(xlabel)
        ylabel = params_name[1]
        if params_unit[1]:
            ylabel += f' ({params_unit[1]})'
        if sampling_mode[1] == 'log':
            self._set_ylogscale(ax, extent, ylabel)
        else:
            ax.set_ylabel(ylabel)
        cbar = plt.colorbar(mmap, ax=ax, extend='max')
        cbar.set_label('misfit')
        figfile_base = self._get_figfile_base(config)
        params_string = f'{params_name[0]}-{params_name[1]}'
        figfile_base += f".misfit_{params_string}_{label.replace(' ', '_')}."
        fmt = config.plot_save_format
        if fmt == 'pdf_multipage':
            fmt = 'pdf'
        figfile = figfile_base + fmt
        if config.plot_show:
            plt.show()
        if config.plot_save:
            savefig(fig, figfile, fmt, bbox_inches='tight')
            if not config.plot_show:
                plt.close(fig)
            logger.info(f'{label}: conditional misfit map saved to: {figfile}')
            config.figures[f'misfit_{params_string}'].append(figfile)

    def _get_figfile_base(self, config):
        outdir = os.path.join(config.options.outdir, 'misfit')
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        evid = config.event.event_id
        return os.path.join(outdir, evid)

    def _set_ylogscale(self, ax, extent, ylabel):
        # the grid is plotted by imshow() in linear scale,
        # so we create a fake yaxis with logscale
        ax.yaxis.set_visible(False)
        ax2 = ax.twinx()
        ax2.yaxis.set_label_position('left')
        ax2.yaxis.set_ticks_position('left')
        ax2.set_ylim(10**extent[2], 10**extent[3])
        ax2.set_yscale('log')
        ax2.set_ylabel(ylabel)

    def _set_xlogscale(self, ax, extent, xlabel):
        # the grid is plotted by imshow() in linear scale,
        # so we create a fake xaxis with logscale
        ax.xaxis.set_visible(False)
        ax2 = ax.twiny()
        ax2.xaxis.set_label_position('bottom')
        ax2.xaxis.set_ticks_position('bottom')
        ax2.set_xlim(10**extent[0], 10**extent[1])
        ax2.set_xscale('log')
        ax2.set_xlabel(xlabel)

    def _plot_kdtree_structure(self, ax, plot_par_idx, params_opt_all):
        coords = np.array([cell.coords for cell in self.kdt.cells])
        # find the complement to plot_par_idx
        allidx = np.arange(len(self.nsteps))
        ii = allidx[~np.isin(allidx, plot_par_idx)]
        # find coordinates that will be discarded
        cc = np.take(coords, ii, axis=1)
        # find optimal parameters that will be discarded
        pp = np.take(params_opt_all, ii)
        # find discarded coordinates that are closer to
        # discarded optimal parameters, i.e. the coordinates that
        # lay on the 2D plot plane
        dist = np.abs(cc - pp)
        # init condition to False
        cond = np.zeros_like(dist[:, 0]).astype(bool)
        for d in dist.T:
            _cond = np.isclose(d, d.min())
            cond = np.logical_or(cond, _cond)
        coords_2d = coords[cond]
        # -- End of code to find coordinates that lay on the 2D plot plane
        # now take coordinates that will be plotted
        # we plot all coords in gray
        coords = np.take(coords, plot_par_idx, axis=1)
        ax.scatter(
            coords[:, 0], coords[:, 1], s=2,
            facecolor='gray', edgecolors='none'
        )
        # coords which lay on the plane are plot in black and larger symbol
        coords_2d = np.take(coords_2d, plot_par_idx, axis=1)
        ax.scatter(
            coords_2d[:, 0], coords_2d[:, 1], s=8,
            facecolor='k', edgecolors='none'
        )
