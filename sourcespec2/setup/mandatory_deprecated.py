# -*- coding: utf8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Mandatory and deprecated config parameters for sourcespec.

:copyright:
    2013-2025 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""

# Mandatory parameters, which cannot be None
mandatory_config_params = [
    'p_arrival_tolerance',
    's_arrival_tolerance',
    'signal_pre_time',
    'win_length',
    'taper_halfwidth',
    'spectral_smooth_width_decades',
    'bp_freqmin_acc',
    'bp_freqmax_acc',
    'bp_freqmin_shortp',
    'bp_freqmax_shortp',
    'bp_freqmin_broadb',
    'bp_freqmax_broadb',
    'freq1_acc',
    'freq2_acc',
    'freq1_shortp',
    'freq2_shortp',
    'freq1_broadb',
    'freq2_broadb',
    'rmsmin',
    'sn_min',
    'spectral_sn_min',
    'rpp',
    'rps',
    'geom_spread_n_exponent',
    'geom_spread_cutoff_distance',
    'f_weight',
    'weight',
    't_star_0',
    't_star_0_variability',
    'a', 'b', 'c',
    'ml_bp_freqmin',
    'ml_bp_freqmax',
    'n_sigma',
    'lower_percentage',
    'mid_percentage',
    'upper_percentage',
    'nIQR',
    'plot_spectra_maxrows',
    'plot_traces_maxrows',
    'plot_station_text_size'
]


def check_deprecated_config_params(config_obj):
    """
    Check if deprecated config parameters are used in the given config object.

    Return a list of deprecation messages for each deprecated parameter found.
    """
    deprecation_msgs = []
    if 's_win_length' in config_obj or 'noise_win_length' in config_obj:
        deprecation_msgs.append(
            '> "s_win_length" and "noise_win_length" config parameters '
            'are no more\n'
            '   supported. Both are replaced by "win_length".\n'
        )
    if 'traceids' in config_obj:
        deprecation_msgs.append(
            '> "traceids" config parameter has been renamed to '
            '"traceid_mapping_file".\n'
        )
    if 'ignore_stations' in config_obj or 'use_stations' in config_obj:
        deprecation_msgs.append(
            '> "ignore_stations" and "use_stations" config parameters '
            'have been renamed to\n'
            '  "ignore_traceids" and "use_traceids", respectively.\n'
        )
    if 'dataless' in config_obj:
        deprecation_msgs.append(
            '> "dataless" config parameter has been renamed to '
            '"station_metadata".\n'
        )
    if 'clip_nmax' in config_obj:
        deprecation_msgs.append(
            '> "clip_nmax" config parameter has been renamed to '
            '"clip_max_percent".\n'
            '   Note that the new default is 5% (current value in your config '
            f'file: {config_obj["clip_nmax"]}%)\n'
        )
    if 'trace_format' in config_obj:
        deprecation_msgs.append(
            '> "trace_format" config parameter is no more supported.\n'
            '   Use "sensitivity" to manually specify how sensor sensitivity '
            'should be computed.\n'
        )
    if 'PLOT_SHOW' in config_obj:
        deprecation_msgs.append(
            '> "PLOT_SHOW" config parameter has been renamed to "plot_show".\n'
        )
    if 'PLOT_SAVE' in config_obj:
        deprecation_msgs.append(
            '> "PLOT_SAVE" config parameter has been renamed to "plot_save".\n'
        )
    if 'PLOT_SAVE_FORMAT' in config_obj:
        deprecation_msgs.append(
            '> "PLOT_SAVE_FORMAT" config parameter has been renamed to '
            '"plot_save_format".\n'
        )
    if 'vp' in config_obj and 'vp_source' not in config_obj:
        deprecation_msgs.append(
            '> "vp" config parameter has been renamed to "vp_source".\n'
        )
    if 'vs' in config_obj and 'vs_source' not in config_obj:
        deprecation_msgs.append(
            '> "vs" config parameter has been renamed to "vs_source".\n'
        )
    if 'rho' in config_obj and 'rho_source' not in config_obj:
        deprecation_msgs.append(
            '> "rho" config parameter has been renamed to "rho_source".\n'
        )
    if 'pre_p_time' in config_obj:
        deprecation_msgs.append(
            '> "pre_p_time" config parameter has been renamed to '
            '"noise_pre_time".\n'
        )
    if 'pre_s_time' in config_obj:
        deprecation_msgs.append(
            '> "pre_s_time" config parameter has been renamed to '
            '"signal_pre_time".\n'
        )
    if 'rps_from_focal_mechanism' in config_obj:
        deprecation_msgs.append(
            '> "rps_from_focal_mechanism" config parameter has been renamed '
            'to "rp_from_focal_mechanism".\n'
        )
    if 'paz' in config_obj:
        deprecation_msgs.append(
            '> "paz" config parameter has been removed and merged with '
            '"station_metadata".\n'
        )
    if 'max_epi_dist' in config_obj:
        deprecation_msgs.append(
            '> "max_epi_dist" config parameter has been removed and replaced '
            'by "epi_dist_ranges".\n'
        )
    if 'max_freq_Er' in config_obj:
        deprecation_msgs.append(
            '> "max_freq_Er" config parameter has been removed and replaced '
            'by "Er_freq_min_max".\n'
        )
    if 'pi_misfit_max' in config_obj:
        deprecation_msgs.append(
            '> "pi_misfit_max" config parameter has been renamed to '
            '"pi_quality_of_fit_min" and its meaning has been reversed.\n'
            '   pi_quality_of_fit_min is the minimum acceptable quality of '
            'fit in percent (0-100).\n'
        )
    return deprecation_msgs
