# -*- coding: utf-8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Plot stacked spectra, along with summary inverted spectrum.

:copyright:
    2023 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
import os
import logging
import numpy as np
from sourcespec.ssp_util import moment_to_mag, mag_to_moment
from sourcespec.ssp_spectral_model import spectral_model
from sourcespec.savefig import savefig
from sourcespec._version import get_versions
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patheffects as PathEffects
logger = logging.getLogger(__name__.split('.')[-1])
# Reduce logging level for Matplotlib to avoid DEBUG messages
mpl_logger = logging.getLogger('matplotlib')
mpl_logger.setLevel(logging.WARNING)


def _summary_synth_spec(sspec_output, fmin, fmax):
    npts = 100
    freq = np.logspace(np.log10(fmin), np.log10(fmax), npts)
    summary_values = sspec_output.reference_values()
    params_name = ('Mw', 'fc', 't_star')
    summary_params = {p: summary_values[p] for p in params_name}
    synth_model_mag = spectral_model(freq, **summary_params)
    synth_model = mag_to_moment(synth_model_mag)
    return freq, synth_model


def _make_fig(config):
    matplotlib.rcParams['pdf.fonttype'] = 42  # to edit text in Illustrator
    figsize = (7, 7)
    if config.plot_save_format in ['pdf', 'pdf_multipage', 'svg']:
        dpi = 72
    else:
        dpi = 200
    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    ax.grid(True, which='both', linestyle='solid', color='#DDDDDD', zorder=0)
    ax.set_xlabel('Frequency (Hz)')
    ax.set_ylabel('Seismic moment (Nm)')
    return fig, ax


def _make_ax2(ax):
    ax2 = ax.twinx()
    moment_minmax = ax.get_ylim()
    mag_minmax = moment_to_mag(moment_minmax)
    ax2.set_ylim(mag_minmax)
    ax2.yaxis.set_tick_params(which='both', labelright=True, pad=0, width=2)
    ax2.set_ylabel('Magnitude')


def _nspectra_text(spec_st, ax):
    nspectra = len(spec_st)
    if nspectra == 1:
        text = '{} inverted spectrum'.format(nspectra)
    else:
        text = '{} inverted spectra'.format(nspectra)
    # Path effect to contour text in white
    path_effects = [PathEffects.withStroke(linewidth=3, foreground='white')]
    color = 'red'
    ax.text(
        0.95, 0.95, text,
        horizontalalignment='right',
        verticalalignment='top',
        color=color,
        fontsize=9,
        transform=ax.transAxes,
        zorder=50,
        path_effects=path_effects)


def _summary_params_text(sspec_output, ax):
    summary_values = sspec_output.reference_values()
    params_name = ('Mo', 'Mw', 'fc', 't_star')
    summary_params = {p: summary_values[p] for p in params_name}
    params_text = 'Summary parameters:\n'
    params_text += 'Mo: {:.1e} Mw: {:.2f}\n'.format(
        summary_params['Mo'], summary_params['Mw'])
    params_text += 'fc: {:.2f} t*: {:.2f}\n'.format(
        summary_params['fc'], summary_params['t_star'])
    color = 'black'
    # Path effect to contour text in white
    path_effects = [PathEffects.withStroke(linewidth=3, foreground='white')]
    ax.text(
        0.05, 0.01, params_text,
        horizontalalignment='left',
        verticalalignment='bottom',
        color=color,
        fontsize=9,
        transform=ax.transAxes,
        zorder=50,
        path_effects=path_effects)


def _add_title(config, ax):
    # Add event information as a title
    hypo = config.hypo
    textstr = 'evid: {}\nlon: {:.3f} lat: {:.3f} depth: {:.1f} km'
    textstr = textstr.format(
        hypo.evid, hypo.longitude, hypo.latitude, hypo.depth)
    try:
        textstr += '\ntime: {}'.format(
            hypo.origin_time.format_iris_web_service())
    except AttributeError:
        pass
    ax.text(0., 1.15, textstr, fontsize=10, linespacing=1.5,
            ha='left', va='top', transform=ax.transAxes)


def _add_code_author(config, ax):
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
    ypos = -0.15
    ax.text(1., ypos, textstr, fontsize=8, linespacing=1.5,
            ha='right', va='top', transform=ax.transAxes)


def _savefig(config, fig):
    evid = config.hypo.evid
    figfile_base = os.path.join(config.options.outdir, evid)
    figfile_base += '.stacked_spectra.'
    fmt = config.plot_save_format
    if fmt == 'pdf_multipage':
        fmt = 'pdf'
    figfile = figfile_base + fmt
    if config.plot_show:
        plt.show()
    if config.plot_save:
        savefig(fig, figfile, fmt, bbox_inches='tight')
        config.figures['stacked_spectra'].append(figfile)
        logger.info('Stacked spectra plot saved to: ' + figfile)


def plot_stacked_spectra(config, spec_st, sspec_output):
    """
    Plot stacked spectra, along with summary inverted spectrum.
    """
    # Check config, if we need to plot at all
    if not config.plot_show and not config.plot_save:
        return
    # Select "H" spectra
    spec_st = spec_st.select(channel='??H')
    # plotting
    fig, ax = _make_fig(config)
    color = 'red'
    alpha = 0.9
    linewidth = 2
    fmins = list()
    fmaxs = list()
    for spec in spec_st:
        freqs = spec.get_freq()
        ax.loglog(
            freqs, spec.data, color=color, lw=linewidth,
            alpha=alpha, zorder=20)
        fmins.append(freqs.min())
        fmaxs.append(freqs.max())
    fmin = min(fmins)
    fmax = max(fmaxs)
    freq, synth_model = _summary_synth_spec(sspec_output, fmin, fmax)
    color = 'black'
    ax.loglog(
        freq, synth_model, color=color, lw=linewidth,
        alpha=alpha, zorder=40)
    _make_ax2(ax)
    _nspectra_text(spec_st, ax)
    _summary_params_text(sspec_output, ax)
    _add_title(config, ax)
    _add_code_author(config, ax)
    _savefig(config, fig)
