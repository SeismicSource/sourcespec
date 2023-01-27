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
    summary_params_no_att = summary_params.copy()
    summary_params_no_att['t_star'] = 0
    synth_model_mag_no_att = spectral_model(freq, **summary_params_no_att)
    synth_model_no_att = mag_to_moment(synth_model_mag_no_att)
    summary_params_no_fc = summary_params.copy()
    summary_params_no_fc['fc'] = 1e999
    synth_model_mag_no_fc = spectral_model(freq, **summary_params_no_fc)
    synth_model_no_fc = mag_to_moment(synth_model_mag_no_fc)
    return freq, synth_model, synth_model_no_att, synth_model_no_fc


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


def _plot_fc_and_mw(sspec_output, ax, ax2):
    fc = sspec_output.reference_values()['fc']
    fc_err_left, fc_err_right = sspec_output.reference_uncertainties()['fc']
    fc_min = fc - fc_err_left
    if fc_min < 0:
        fc_min = 0.01
    ax.axvspan(
        fc_min, fc+fc_err_right, color='#bbbbbb', alpha=0.3, zorder=9)
    ax.axvline(fc, color='#999999', linewidth=2., zorder=10)
    Mw = sspec_output.reference_values()['Mw']
    Mw_err_left, Mw_err_right = sspec_output.reference_uncertainties()['Mw']
    ax2.axhspan(
        Mw-Mw_err_left, Mw+Mw_err_right, color='#bbbbbb',
        alpha=0.3, zorder=8)


def _make_ax2(ax):
    ax2 = ax.twinx()
    moment_minmax = ax.get_ylim()
    mag_minmax = moment_to_mag(moment_minmax)
    ax2.set_ylim(mag_minmax)
    ax2.yaxis.set_tick_params(which='both', labelright=True, pad=0, width=2)
    ax2.set_ylabel('Magnitude')
    return ax2


def _nspectra_text(spec_st):
    nspectra = len(spec_st)
    if nspectra == 1:
        text = 'Inverted spectrum ({})'.format(nspectra)
    else:
        text = 'Inverted spectra ({})'.format(nspectra)
    return text


def _summary_params_text(sspec_output, ax):
    summary_values = sspec_output.reference_values()
    summary_uncertainties = sspec_output.reference_uncertainties()
    params_name = ('Mo', 'Mw', 'fc', 't_star')
    summary_params = {p: summary_values[p] for p in params_name}
    summary_errors = {p: summary_uncertainties[p] for p in params_name}
    Mo_text = 'Mo: {:.1e} [- {:.1e}, + {:.1e}] Nm'.format(
        summary_params['Mo'], *summary_errors['Mo'])
    Mw_text = 'Mw: {:.2f} [- {:.2f}, + {:.2f}]'.format(
        summary_params['Mw'], *summary_errors['Mw'])
    fc_text = 'fc: {:.2f} [- {:.2f}, + {:.2f}] Hz'.format(
        summary_params['fc'], *summary_errors['fc'])
    t_star_text = 't*: {:.2f} [- {:.2f}, + {:.2f}] s'.format(
        summary_params['t_star'], *summary_errors['t_star'])
    params_text = 'Summary parameters:\n'
    params_text += '{}\n{}\n{}\n{}'.format(
        Mo_text, Mw_text, fc_text, t_star_text)
    color = 'black'
    # Path effect to contour text in white
    path_effects = [PathEffects.withStroke(linewidth=3, foreground='white')]
    ax.text(
        0.05, 0.02, params_text,
        horizontalalignment='left',
        verticalalignment='bottom',
        color=color,
        fontsize=9,
        linespacing=1.5,
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
    alpha = 0.5
    linewidth = 2
    fmins = list()
    fmaxs = list()
    for spec in spec_st:
        freqs = spec.get_freq()
        spec_handle, = ax.loglog(
            freqs, spec.data, color=color, lw=linewidth,
            alpha=alpha, zorder=20)
        fmins.append(freqs.min())
        fmaxs.append(freqs.max())
    fmin = min(fmins)
    fmax = max(fmaxs)
    freq, synth_model, synth_model_no_att, synth_model_no_fc =\
        _summary_synth_spec(sspec_output, fmin, fmax)
    color = 'black'
    alpha = 0.9
    linewidth = 2
    synth_handle, = ax.loglog(
        freq, synth_model, color=color, lw=linewidth,
        alpha=alpha, zorder=40)
    legend_handles = [spec_handle, synth_handle]
    legend_labels = [
        _nspectra_text(spec_st),
        'Summary parameters'
    ]
    if config.plot_spectra_no_attenuation:
        color = 'gray'
        linewidth = 2
        synth_no_att_handle, = ax.loglog(
            freq, synth_model_no_att, color=color, lw=linewidth,
            alpha=alpha, zorder=39)
        legend_handles.append(synth_no_att_handle)
        legend_labels.append('Summary parameters,\nno attenuation')
    if config.plot_spectra_no_fc:
        color = 'gray'
        linewidth = 2
        linestyle = 'dashed'
        synth_no_fc_handle, = ax.loglog(
            freq, synth_model_no_fc, color=color, lw=linewidth, ls=linestyle,
            alpha=alpha, zorder=39)
        legend_handles.append(synth_no_fc_handle)
        legend_labels.append('Summary parameters,\nno fc')
    ax.legend(
        legend_handles, legend_labels,
        loc='upper right', framealpha=1
    ).set_zorder(100)
    ax2 = _make_ax2(ax)
    _plot_fc_and_mw(sspec_output, ax, ax2)
    _summary_params_text(sspec_output, ax)
    _add_title(config, ax)
    _add_code_author(config, ax)
    _savefig(config, fig)
