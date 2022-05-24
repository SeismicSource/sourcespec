# -*- coding: utf-8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Plot parameter statistics.

:copyright:
    2022 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patheffects as mpe
import logging
from sourcespec.savepng import savepng
matplotlib.rcParams['pdf.fonttype'] = 42  # to edit text in Illustrator
logger = logging.getLogger(__name__.split('.')[-1])
# Reduce logging level for Matplotlib to avoid DEBUG messages
mpl_logger = logging.getLogger('matplotlib')
mpl_logger.setLevel(logging.WARNING)


def box_plots(config, sourcepar):
    """Show parameter statistics through box plots."""
    # Check config, if we need to plot at all
    if not config.PLOT_SHOW and not config.PLOT_SAVE:
        return

    param_names_units_colors = {
        'Mw': ('Mw', None, '#EE5835'),
        'fc': ('Corner Frequency', 'Hz', '#6FBA6C'),
        't_star': ('t*', 's', '#9EBAE2'),
        'ra': ('Source Radius', 'm', '#FAAC64'),
        'bsd': ('Brune Stress Drop', 'MPa', '#D4ADD2'),
        'Qo': ('Quality Factor', None, '#C07131'),
        'Er': ('Radiated Energy', 'N.m', '#00E3E9'),
        'Ml': ('Ml', None, '#FC8384')
    }
    npars = len(param_names_units_colors)

    path_effect = [
        mpe.Stroke(linewidth=5, foreground='black'),
        mpe.Stroke(foreground='white', alpha=1),
        mpe.Normal()
    ]

    fig, axes = plt.subplots(npars, 1, figsize=(8, 8))
    fig.set_tight_layout(True)
    for ax, (param, (name, unit, color)) in zip(
            axes, param_names_units_colors.items()):
        values = sourcepar.value_array(param)
        # remove nan and inf values
        values = values[~np.isnan(values)]
        values = values[~np.isinf(values)]
        if len(values) == 0:
            ax.set_visible(False)
            continue
        vmean = np.mean(values)
        # remove values that are very far from the mean
        # (otherwhise plot is too compressed)
        values = values[np.abs(values-vmean)/vmean < 10]
        whis = config.nIQR
        if whis is None:
            whis = (0, 100)
        bplot = ax.boxplot(
            values, vert=False, whis=whis,
            widths=0.7, patch_artist=True)
        ax.scatter(
            values, np.ones_like(values), color='dimgray', edgecolor='white',
            zorder=10)
        for box in bplot['boxes']:
            box.set_facecolor(color)
        for line in bplot['medians']:
            line.set_color('white')
            line.set_linewidth(2)
            line.set_path_effects(path_effect)
        if unit is None:
            ax.set_xlabel('{}'.format(name))
        else:
            ax.set_xlabel('{} ({})'.format(name, unit))
        ax.tick_params(left=False, labelleft=False)

    evid = config.hypo.evid
    figfile_base = os.path.join(config.options.outdir, evid)
    figfile_base += '.boxplot.'
    fmt = config.PLOT_SAVE_FORMAT
    if fmt == 'pdf_multipage':
        fmt = 'pdf'
    figfile = figfile_base + fmt
    if config.PLOT_SHOW:
        plt.show()
    if config.PLOT_SAVE:
        if fmt == 'png':
            savepng(fig, figfile, bbox_inches='tight')
        else:
            fig.savefig(figfile, bbox_inches='tight')
        logger.info('Parameters box plot saved to: ' + figfile)
