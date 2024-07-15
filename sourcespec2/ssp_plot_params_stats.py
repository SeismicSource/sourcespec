# -*- coding: utf-8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Plot parameter statistics.

:copyright:
    2022-2026 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
import os
import contextlib
import logging
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patheffects as mpe
from .setup import config
from ._version import get_versions
from .savefig import savefig
logger = logging.getLogger(__name__.rsplit('.', maxsplit=1)[-1])
# Reduce logging level for Matplotlib to avoid DEBUG messages
logging.getLogger('matplotlib').setLevel(logging.WARNING)


class PlotParam():
    """A plot parameter."""
    def __init__(self, name, unit, color):
        self.name = name
        self.unit = unit
        self.color = color


def box_plots(sspec_output):
    """
    Show parameter statistics through box plots.

    :param sspec_output: SourceSpec output.
    :type sspec_output: :class:`sourcespec.ssp_data_types.SourceSpecOutput`
    """
    # Check config, if we need to plot at all
    if not config.plot_show and not config.plot_save:
        return
    matplotlib.rcParams['pdf.fonttype'] = 42  # to edit text in Illustrator

    plot_params = {
        'Mw': PlotParam('Mw', None, '#EE5835'),
        'fc': PlotParam('Corner Frequency', 'Hz', '#6FBA6C'),
        't_star': PlotParam('t-star', 's', '#9EBAE2'),
        'radius': PlotParam('Source Radius', 'm', '#FAAC64'),
        'ssd': PlotParam('Static Stress Drop', 'MPa', '#D4ADD2'),
        'Qo': PlotParam('Quality Factor', None, '#C07131'),
        'Er': PlotParam('Radiated Energy', 'N·m', '#00E3E9'),
        'sigma_a': PlotParam('Apparent Stress', 'MPa', '#943B99'),
        'Ml': PlotParam('Ml', None, '#FC8384')
    }
    npars = len(plot_params)

    path_effect = [
        mpe.Stroke(linewidth=5, foreground='black'),
        mpe.Stroke(foreground='white', alpha=1),
        mpe.Normal()
    ]

    fig, axes = plt.subplots(npars, 1, figsize=(8, 10), dpi=300)
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
        _newvalues = values[np.abs(values - vmean) / vmean < 10]
        if len(_newvalues) > 0:
            values = _newvalues
        # recompute mean after removing outliers
        vmean = np.mean(values)
        min_values = np.min(values)
        max_values = np.max(values)
        # use logscale if values span a large interval, avoid zero values
        logscale = False
        if max_values / min_values > 50:
            if min_values == 0:
                min_values = 1e-10
            values = values[values != 0]
            ax.set_xscale('log')
            logscale = True
        whiskers = config.nIQR
        if whiskers is None:
            whiskers = (0, 100)
        bplot = ax.boxplot(
            values, vert=False, whis=whiskers,
            widths=0.7, patch_artist=True)
        ax.scatter(
            values, np.ones_like(values), color='dimgray', edgecolor='white',
            zorder=10)
        # if values are very close to each other, add some margin
        margin = 0.001 * vmean
        if max_values - min_values < margin:
            max_values += margin
            # avoid zero or negative min_values in logscale
            if logscale:
                _min_values = min_values - margin
                while _min_values <= 0:
                    margin /= 2
                    _min_values = min_values - margin
            min_values -= margin
            ax.set_xlim(min_values, max_values)
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
    evid = config.event.event_id
    hypo = config.event.hypocenter
    ev_lon = hypo.longitude.value_in_deg
    ev_lat = hypo.latitude.value_in_deg
    ev_depth = hypo.depth.value_in_km
    textstr = f'{config.options.evname} — ' if config.options.evname else ''
    textstr += (
        f'evid: {evid}\n'
        f'lon: {ev_lon:.3f} lat: {ev_lat:.3f} depth: {ev_depth:.1f} km'
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
        _savefig(fig)


def _savefig(fig):
    """
    Save the figure to a file.

    :param fig: Figure to save.
    :type fig: :class:`matplotlib.figure.Figure`
    """
    evid = config.event.event_id
    figfile_base = os.path.join(config.options.outdir, f'{evid}.boxplot.')
    fmt = config.plot_save_format
    if fmt == 'pdf_multipage':
        fmt = 'pdf'
    figfile = f'{figfile_base}{fmt}'
    savefig(fig, figfile, fmt, bbox_inches='tight')
    config.figures['boxplots'].append(figfile)
    logger.info(f'Parameters box plot saved to: {figfile}')
