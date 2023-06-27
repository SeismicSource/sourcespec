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
import contextlib
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


def _spectral_model_moment(freq, Mw, fc, t_star):
    return mag_to_moment(spectral_model(freq, Mw=Mw, fc=fc, t_star=t_star))


def _summary_synth_spec(sspec_output, fmin, fmax):
    npts = 100
    freq = np.logspace(np.log10(fmin), np.log10(fmax), npts)
    summary_values = sspec_output.reference_values()
    Mw = summary_values['Mw']
    fc = summary_values['fc']
    t_star = summary_values['t_star']
    synth_model = _spectral_model_moment(freq, Mw, fc, t_star)
    synth_model_no_att = _spectral_model_moment(freq, Mw, fc, 0)
    synth_model_no_fc = _spectral_model_moment(freq, Mw, 1e999, t_star)
    return freq, synth_model, synth_model_no_att, synth_model_no_fc


def _make_fig(config):
    matplotlib.rcParams['pdf.fonttype'] = 42  # to edit text in Illustrator
    figsize = (7, 7)
    dpi = (
        72 if config.plot_save_format in ['pdf', 'pdf_multipage', 'svg']
        else 200
    )
    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    ax.grid(True, which='both', linestyle='solid', color='#DDDDDD', zorder=0)
    ax.set_xlabel('Frequency (Hz)')
    ax.set_ylabel('Seismic moment (N·m)')
    return fig, ax


def _plot_fc_and_mw(sspec_output, ax, ax2):
    fc = sspec_output.reference_values()['fc']
    fc_err_left, fc_err_right = sspec_output.reference_uncertainties()['fc']
    fc_min = fc - fc_err_left
    if fc_min < 0:
        fc_min = 0.01
    ax.axvspan(
        fc_min, fc + fc_err_right, color='#bbbbbb', alpha=0.3, zorder=9)
    ax.axvline(fc, color='#999999', linewidth=2., zorder=10)
    Mw = sspec_output.reference_values()['Mw']
    Mw_err_left, Mw_err_right = sspec_output.reference_uncertainties()['Mw']
    ax2.axhspan(
        Mw - Mw_err_left, Mw + Mw_err_right, color='#bbbbbb',
        alpha=0.3, zorder=8)


def _make_ax2(ax):
    ax2 = ax.twinx()
    # Move ax2 below ax
    ax.zorder = 2
    ax2.zorder = 1
    ax.patch.set_visible(False)
    ax2.patch.set_visible(True)
    moment_minmax = ax.get_ylim()
    mag_minmax = moment_to_mag(moment_minmax)
    ax2.set_ylim(mag_minmax)
    ax2.yaxis.set_tick_params(which='both', labelright=True, pad=0, width=2)
    ax2.set_ylabel('Magnitude')
    return ax2


def _nspectra_text(spec_list):
    nspectra = len(spec_list)
    return (
        'Inverted spectrum'
        if nspectra == 1
        else f'Inverted spectra ({nspectra})'
    )


def _summary_params_text(sspec_output, ax):
    summary_values = sspec_output.reference_values()
    summary_uncertainties = sspec_output.reference_uncertainties()
    Mo_value = summary_values['Mo']
    Mo_err_left, Mo_err_right = summary_uncertainties['Mo']
    Mo_text = (
        f'Mo: {Mo_value:.1e} [- {Mo_err_left:.1e}, + {Mo_err_right:.1e}] N·m')
    Mw_value = summary_values['Mw']
    Mw_err_left, Mw_err_right = summary_uncertainties['Mw']
    Mw_text = (
        f'Mw: {Mw_value:.2f} [- {Mw_err_left:.2f}, + {Mw_err_right:.2f}]')
    fc_value = summary_values['fc']
    fc_err_left, fc_err_right = summary_uncertainties['fc']
    fc_text = (
        f'fc: {fc_value:.2f} [- {fc_err_left:.2f}, + {fc_err_right:.2f}] Hz')
    t_star_value = summary_values['t_star']
    t_star_err_left, t_star_err_right = summary_uncertainties['t_star']
    t_star_text = (
        f't*: {t_star_value:.2f} '
        f'[- {t_star_err_left:.2f}, + {t_star_err_right:.2f}] s')
    params_text = (
        f'Summary parameters:\n{Mo_text}\n{Mw_text}\n{fc_text}\n{t_star_text}'
    )
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
    evid = config.event.event_id
    hypo = config.event.hypocenter
    ev_lon = hypo.longitude.value_in_deg
    ev_lat = hypo.latitude.value_in_deg
    ev_depth = hypo.depth.value_in_km
    textstr = (
        f'evid: {evid}\nlon: {ev_lon:.3f} lat: {ev_lat:.3f} '
        f'depth: {ev_depth:.1f} km'
    )
    with contextlib.suppress(AttributeError):
        textstr += f'\ntime: {hypo.origin_time.format_iris_web_service()}'
    ax.text(
        0., 1.15, textstr, fontsize=10, linespacing=1.5,
        ha='left', va='top', transform=ax.transAxes)


def _add_code_author(config, ax):
    # Add code and author information to the figure bottom
    textstr = (
        f'SourceSpec v{get_versions()["version"]} '
        f'- {config.end_of_run.strftime("%Y-%m-%d %H:%M:%S")} '
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
    ypos = -0.15
    ax.text(
        1., ypos, textstr, fontsize=8, linespacing=1.5,
        ha='right', va='top', transform=ax.transAxes)


def _savefig(config, fig):
    evid = config.event.event_id
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
        logger.info(f'Stacked spectra plot saved to: {figfile}')


def plot_stacked_spectra(config, spec_st, sspec_output):
    """
    Plot stacked spectra, along with summary inverted spectrum.
    """
    # Check config, if we need to plot at all
    if not config.plot_show and not config.plot_save:
        return
    # Select "H" spectra (note: spec_st.select cannot be used because it is
    # case unsensitive)
    selected_specs = [
        spec for spec in spec_st if spec.stats.channel[-1] == 'H']
    # plotting
    fig, ax = _make_fig(config)
    color = 'red'
    alpha = 0.5
    linewidth = 2
    fmins = []
    fmaxs = []
    for spec in selected_specs:
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
        _nspectra_text(selected_specs),
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
