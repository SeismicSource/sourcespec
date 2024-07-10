# -*- coding: utf-8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Plot stacked spectra, along with summary inverted spectrum.

:copyright:
    2023-2025 Claudio Satriano <satriano@ipgp.fr>
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
import matplotlib.patheffects as PathEffects
from matplotlib.collections import LineCollection
from .config import config
from .ssp_util import moment_to_mag, mag_to_moment
from .ssp_spectral_model import spectral_model
from .savefig import savefig
from ._version import get_versions
logger = logging.getLogger(__name__.rsplit('.', maxsplit=1)[-1])
# Reduce logging level for Matplotlib to avoid DEBUG messages
logging.getLogger('matplotlib').setLevel(logging.WARNING)


def _spectral_model_moment(freq, Mw, fc, t_star):
    """
    Compute spectral model in moment units.

    :param freq: Frequency array
    :type freq: :class:`numpy.ndarray`
    :param Mw: Moment magnitude
    :type Mw: float
    :param fc: Corner frequency
    :type fc: float
    :param t_star: Attenuation parameter
    :type t_star: float

    :return: Spectral model in moment units
    :rtype: :class:`numpy.ndarray`
    """
    return mag_to_moment(spectral_model(freq, Mw=Mw, fc=fc, t_star=t_star))


def _summary_synth_spec(sspec_output, fmin, fmax):
    """
    Compute synthetic spectrum for summary parameters.

    :param sspec_output: Output of spectral inversion
    :type sspec_output: :class:`~sourcespec.ssp_data_types.SourceSpecOutput`
    :param fmin: Minimum frequency
    :type fmin: float
    :param fmax: Maximum frequency
    :type fmax: float

    :return: Frequency array, synthetic spectrum, synthetic spectrum without
        attenuation, synthetic spectrum without corner frequency
    :rtype: tuple of :class:`numpy.ndarray`
    """
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


def _make_fig():
    """
    Create figure and axes for stacked spectra plot.

    :return: Figure and axes
    :rtype: tuple of :class:`matplotlib.figure.Figure`,
        :class:`matplotlib.axes.Axes`
    """
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
    """
    Plot corner frequency and moment magnitude on the axes.

    :param sspec_output: Output of spectral inversion
    :type sspec_output: :class:`~sourcespec.ssp_data_types.SourceSpecOutput`
    :param ax: Axes for stacked spectra plot
    :type ax: :class:`matplotlib.axes.Axes`
    :param ax2: Axes for magnitude plot
    :type ax2: :class:`matplotlib.axes.Axes`
    """
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
    """
    Create a second y-axis for magnitude units.

    :param ax: Axes for stacked spectra plot
    :type ax: :class:`matplotlib.axes.Axes`

    :return: Axes for magnitude units
    :rtype: :class:`matplotlib.axes.Axes`
    """
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
    """
    Return text for the number of spectra.

    :param spec_list: List of spectra
    :type spec_list: list

    :return: Text for the number of spectra
    :rtype: str
    """
    nspectra = len(spec_list)
    return (
        'Inverted spectrum'
        if nspectra == 1
        else f'Inverted spectra ({nspectra})'
    )


def _summary_params_text(sspec_output, ax):
    """
    Add summary parameters to the plot.

    :param sspec_output: Output of spectral inversion
    :type sspec_output: :class:`~sourcespec.ssp_data_types.SourceSpecOutput`
    :param ax: Axes for stacked spectra plot
    :type ax: :class:`matplotlib.axes.Axes`
    """
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
    Er_value = summary_values['Er']
    Er_err_left, Er_err_right = summary_uncertainties['Er']
    Er_text = (
        f'Er: {Er_value:.1e} [- {Er_err_left:.1e}, + {Er_err_right:.1e}] N·m')
    params_text = (
        'Summary parameters:\n'
        f'{Mo_text}\n{Mw_text}\n{fc_text}\n{t_star_text}\n{Er_text}'
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


def _add_title(ax):
    """
    Add event information as a title to the plot.

    :param ax: Axes for stacked spectra plot
    :type ax: :class:`matplotlib.axes.Axes`
    """
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
    ax.text(
        0., 1.10, textstr, fontsize=10, linespacing=1.5,
        ha='left', va='top', transform=ax.transAxes)


def _add_code_author(ax):
    """
    Add code and author information to the figure bottom.

    :param ax: Axes for stacked spectra plot
    :type ax: :class:`matplotlib.axes.Axes`
    """
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


def _savefig(fig):
    """
    Save figure to file.

    :param fig: Figure to save
    :type fig: :class:`matplotlib.figure.Figure`
    """
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
        savefig(fig, figfile, fmt, bbox_inches='tight', dpi=300)
        config.figures['stacked_spectra'].append(figfile)
        logger.info(f'Stacked spectra plot saved to: {figfile}')


def _truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    """
    Truncate a colormap to a specific range.
    """
    return matplotlib.colors.LinearSegmentedColormap.from_list(
        f'trunc({cmap.name},{minval:.2f},{maxval:.2f})',
        cmap(np.linspace(minval, maxval, n))
    )


def plot_stacked_spectra(spec_st, weight_st, sspec_output):
    """
    Plot stacked spectra, along with summary inverted spectrum.

    :param spec_st: Stream of spectra
    :type spec_st: :class:`~sourcespec.spectrum.SpectrumStream`
    :param weight_st: Stream of spectral weights
    :type weight_st: :class:`~sourcespec.spectrum.SpectrumStream`
    :param sspec_output: Output of spectral inversion
    :type sspec_output: :class:`~sourcespec.ssp_data_types.SourceSpecOutput`
    """
    # Check config, if we need to plot at all
    if not config.plot_show and not config.plot_save:
        return
    # Select "H" spectra (note: spec_st.select cannot be used because it is
    # case unsensitive). Also, ignore spectra with ignore=True
    selected_specs = [
        spec for spec in spec_st
        if spec.stats.channel[-1] == 'H'
        and not spec.stats.ignore
    ]
    # plotting
    fig, ax = _make_fig()
    ax.set_xscale('log')
    ax.set_yscale('log')
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
        'w_red', ['white', 'red'], N=100)
    cmap = _truncate_colormap(cmap, 0.2, 1)
    norm = matplotlib.colors.Normalize(vmin=0, vmax=1)
    alpha = 0.5
    linewidth = 2
    fmins = []
    fmaxs = []
    specmins = []
    specmaxs = []
    for spec in selected_specs:
        # find the same spec in weight_st
        try:
            weight_spec = weight_st.select(id=spec.id)[0]
            weight = weight_spec.data
        except IndexError:
            # this should not happen, but if it does, use a weight of 1
            weight = np.ones_like(spec.data)
        color = cmap(norm(weight))
        freqs = spec.freq
        # create a multi-segment line to use a colormap
        points = np.array([freqs, spec.data]).T.reshape((-1, 1, 2))
        segments = np.concatenate([points[:-1], points[1:]], axis=1)
        # generate a rasterized LineCollection to produce smaller files
        lc = LineCollection(
            segments, cmap=cmap, norm=norm, alpha=alpha, rasterized=True,
            zorder=20)
        lc.set_array(weight)
        lc.set_linewidth(linewidth)
        ax.add_collection(lc)
        # store min/max values for axes limits
        fmins.append(freqs.min())
        fmaxs.append(freqs.max())
        specmins.append(spec.data.min())
        specmaxs.append(spec.data.max())
    fmin = min(fmins)
    fmax = max(fmaxs)
    ax.set_xlim(fmin, fmax)
    specmin = min(specmins)
    specmax = max(specmaxs)
    padding = 0.05*(np.log10(specmax) - np.log10(specmin))
    ax.set_ylim(
        10**(np.log10(specmin)-padding),
        10**(np.log10(specmax)+padding)
    )
    freq, synth_model, synth_model_no_att, synth_model_no_fc =\
        _summary_synth_spec(sspec_output, fmin, fmax)
    color = 'black'
    alpha = 0.9
    linewidth = 2
    synth_handle, = ax.plot(
        freq, synth_model, color=color, lw=linewidth,
        alpha=alpha, zorder=40)
    # draw an invisible line to add to the legend
    spec_handle, = ax.plot(
        [1], [1], color=cmap(norm(0.5)), lw=linewidth,
        alpha=alpha, zorder=20)
    legend_handles = [spec_handle, synth_handle]
    legend_labels = [
        _nspectra_text(selected_specs),
        'Summary parameters'
    ]
    if config.plot_spectra_no_attenuation:
        color = 'gray'
        linewidth = 2
        synth_no_att_handle, = ax.plot(
            freq, synth_model_no_att, color=color, lw=linewidth,
            alpha=alpha, zorder=39)
        legend_handles.append(synth_no_att_handle)
        legend_labels.append('Summary parameters,\nno attenuation')
    if config.plot_spectra_no_fc:
        color = 'gray'
        linewidth = 2
        linestyle = 'dashed'
        synth_no_fc_handle, = ax.plot(
            freq, synth_model_no_fc, color=color, lw=linewidth, ls=linestyle,
            alpha=alpha, zorder=39)
        legend_handles.append(synth_no_fc_handle)
        legend_labels.append('Summary parameters,\nno fc')
    ax.legend(
        legend_handles, legend_labels,
        loc='upper right', framealpha=1
    ).set_zorder(100)
    axins = ax.inset_axes([0.05, 0.27, 0.3, 0.04])
    cbar = matplotlib.colorbar.ColorbarBase(
        axins, cmap=cmap, norm=norm, orientation='horizontal')
    cbar.ax.text(
        0.5, 0.5, 'Weight', fontsize=10, color='w', ha='center', va='center')
    cbar.ax.set_xticks([0, 0.5, 1], ['0', '0.5', '1'], size=8)
    ax2 = _make_ax2(ax)
    _plot_fc_and_mw(sspec_output, ax, ax2)
    _summary_params_text(sspec_output, ax)
    _add_title(ax)
    _add_code_author(ax)
    _savefig(fig)
