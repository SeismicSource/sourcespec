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

from sourcespec._version import get_versions
from sourcespec.savefig import savefig
matplotlib.rcParams['pdf.fonttype'] = 42  # to edit text in Illustrator
logger = logging.getLogger(__name__.split('.')[-1])
# Reduce logging level for Matplotlib to avoid DEBUG messages
mpl_logger = logging.getLogger('matplotlib')
mpl_logger.setLevel(logging.WARNING)


def box_plots(config, sourcepar):
    """Show parameter statistics through box plots."""
    # Check config, if we need to plot at all
    if not config.plot_show and not config.plot_save:
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

    fig, axes = plt.subplots(npars, 1, figsize=(8, 9), dpi=300)
    fig.set_tight_layout(True)
    # Create an invisible axis and use it for title and footer
    ax0 = fig.add_subplot(111, label='ax0')
    ax0.set_axis_off()
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

    # Add event information as a title
    hypo = config.hypo
    textstr = 'evid: {}\nlon: {:.3f} lat: {:.3f} depth: {:.1f} km'
    textstr = textstr.format(
        hypo.evid, hypo.longitude, hypo.latitude, hypo.depth)
    try:
        textstr += ' time: {}'.format(
            hypo.origin_time.format_iris_web_service())
    except AttributeError:
        pass
    ax0.text(0., 1.08, textstr, fontsize=10, linespacing=1.5,
             ha='left', va='top', transform=ax0.transAxes)
    # Add code and author information at the figure bottom
    textstr = 'SourceSpec v{} '.format(get_versions()['version'])
    textstr += '– {} {} '.format(
        config.end_of_run.strftime('%Y-%m-%d %H:%M:%S'),
        config.end_of_run_tz)
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
        textstr = '{}\n{} '.format(textstr, textstr2)
    if not axes[-1].get_visible():
        ypos = 0.04
    else:
        ypos = -0.08
    ax0.text(1., ypos, textstr, fontsize=8, linespacing=1.5,
             ha='right', va='top', transform=ax0.transAxes)

    evid = config.hypo.evid
    figfile_base = os.path.join(config.options.outdir, evid)
    figfile_base += '.boxplot.'
    fmt = config.plot_save_format
    if fmt == 'pdf_multipage':
        fmt = 'pdf'
    figfile = figfile_base + fmt
    if config.plot_show:
        plt.show()
    if config.plot_save:
        savefig(fig, figfile, fmt, bbox_inches='tight')
        config.figures['boxplots'].append(figfile)
        logger.info('Parameters box plot saved to: ' + figfile)
