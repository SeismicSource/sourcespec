# -*- coding: utf-8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Plot parameter statistics.

:copyright:
    2022-2023 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
import os
import contextlib
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patheffects as mpe
import logging

from sourcespec._version import get_versions
from sourcespec.savefig import savefig
matplotlib.rcParams['pdf.fonttype'] = 42  # to edit text in Illustrator
logger = logging.getLogger(__name__.split('.')[-1])
# Reduce logging level for Matplotlib to avoid DEBUG messages
mpl_logger = logging.getLogger('matplotlib')
mpl_logger.setLevel(logging.WARNING)


class PlotParam():
    def __init__(self, name, unit, color):
        self.name = name
        self.unit = unit
        self.color = color


def box_plots(config, sspec_output):
    """Show parameter statistics through box plots."""
    # Check config, if we need to plot at all
    if not config.plot_show and not config.plot_save:
        return

    plot_params = {
        'Mw': PlotParam('Mw', None, '#EE5835'),
        'fc': PlotParam('Corner Frequency', 'Hz', '#6FBA6C'),
        't_star': PlotParam('t-star', 's', '#9EBAE2'),
        'radius': PlotParam('Source Radius', 'm', '#FAAC64'),
        'bsd': PlotParam('Brune Stress Drop', 'MPa', '#D4ADD2'),
        'Qo': PlotParam('Quality Factor', None, '#C07131'),
        'Er': PlotParam('Radiated Energy', 'N·m', '#00E3E9'),
        'Ml': PlotParam('Ml', None, '#FC8384')
    }
    npars = len(plot_params)

    path_effect = [
        mpe.Stroke(linewidth=5, foreground='black'),
        mpe.Stroke(foreground='white', alpha=1),
        mpe.Normal()
    ]

    fig, axes = plt.subplots(npars, 1, figsize=(8, 9), dpi=300)
    fig.set_tight_layout(True)
    # Create an invisible axis and use it for title and footer
    ax0 = fig.add_subplot(111, label='ax0')
    ax0.set_axis_off()
    for ax, param in zip(axes, plot_params.keys()):
        plot_param = plot_params[param]
        values = sspec_output.value_array(param)
        # remove nan and inf values
        values = values[~np.isnan(values)]
        values = values[~np.isinf(values)]
        if len(values) == 0:
            ax.set_visible(False)
            continue
        vmean = np.mean(values)
        # remove values that are very far from the mean
        # (otherwise plot is too compressed)
        values = values[np.abs(values-vmean)/vmean < 10]
        # use logscale if values span a large interval
        if max(values)/min(values) > 50:
            ax.set_xscale('log')
        whiskers = config.nIQR
        if whiskers is None:
            whiskers = (0, 100)
        bplot = ax.boxplot(
            values, vert=False, whis=whiskers,
            widths=0.7, patch_artist=True)
        ax.scatter(
            values, np.ones_like(values), color='dimgray', edgecolor='white',
            zorder=10)
        for box in bplot['boxes']:
            box.set_facecolor(plot_param.color)
        for line in bplot['medians']:
            line.set_color('white')
            line.set_linewidth(2)
            line.set_path_effects(path_effect)
        if plot_param.unit is None:
            ax.set_xlabel(f'{plot_param.name}')
        else:
            ax.set_xlabel(f'{plot_param.name} ({plot_param.unit})')
        ax.tick_params(left=False, labelleft=False)
        ax.minorticks_on()

    # Add event information as a title
    hypo = config.hypo
    textstr = (
        f'evid: {hypo.evid}\nlon: {hypo.longitude:.3f} '
        f'lat: {hypo.latitude:.3f} depth: {hypo.depth:.1f} km'
    )
    with contextlib.suppress(AttributeError):
        textstr += f' time: {hypo.origin_time.format_iris_web_service()}'
    ax0.text(
        0., 1.08, textstr, fontsize=10, linespacing=1.5,
        ha='left', va='top', transform=ax0.transAxes)
    textstr = (
        f'SourceSpec v{get_versions()["version"]} '
        f'– {config.end_of_run.strftime("%Y-%m-%d %H:%M:%S")} '
        f'{config.end_of_run_tz} '
    )
    textstr2 = ''
    if config.author_name is not None:
        textstr2 += config.author_name
    elif config.author_email is not None:
        textstr2 += config.author_email
    if config.agency_short_name is not None:
        if textstr2 != '':
            textstr2 += ' - '
        textstr2 += config.agency_short_name
    elif config.agency_full_name is not None:
        if textstr2 != '':
            textstr2 += ' - '
        textstr2 += config.agency_full_name
    if textstr2 != '':
        textstr = f'{textstr}\n{textstr2} '
    ypos = -0.08 if axes[-1].get_visible() else 0.04
    ax0.text(
        1., ypos, textstr, fontsize=8, linespacing=1.5,
        ha='right', va='top', transform=ax0.transAxes)

    if config.plot_show:
        plt.show()
    if config.plot_save:
        _savefig(config, fig)


def _savefig(config, fig):
    evid = config.hypo.evid
    figfile_base = os.path.join(config.options.outdir, f'{evid}.boxplot.')
    fmt = config.plot_save_format
    if fmt == 'pdf_multipage':
        fmt = 'pdf'
    figfile = f'{figfile_base}{fmt}'
    savefig(fig, figfile, fmt, bbox_inches='tight')
    config.figures['boxplots'].append(figfile)
    logger.info(f'Parameters box plot saved to: {figfile}')
