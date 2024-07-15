# -*- coding: utf8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Mandatory and deprecated config parameters for sourcespec.

:copyright:
    2013-2024 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""

# Mandatory parameters, which cannot be None
mandatory_config_params = [
    'p_arrival_tolerance',
    's_arrival_tolerance',
    'noise_pre_time',
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


# Deprecated parameters (keys) and messages to be shown to the user (values)
# The messages are lists of strings, each string is a line of the message.
deprecated_config_params = {
    's_win_length':
        ['"s_win_length" has been removed and replaced by "win_length".'],
    'noise_win_length':
        ['"noise_win_length" has been removed and replaced by "win_length".'],
    'traceids':
        ['"traceids" has been renamed to "traceid_mapping_file".'],
    'ignore_stations':
        ['"ignore_stations" has been replaced by "ignore_traceids".'],
    'use_stations':
        ['"use_stations" has been replaced by "use_traceids".'],
    'dataless':
        ['"dataless" has been renamed to "station_metadata".'],
    'clip_nmax':
        [
            '"clip_nmax" has been renamed to "clip_max_percent".',
            'Note that the new default is 5%.'
        ],
    'trace_format':
        [
            '"trace_format" is no more supported.',
            'Use "sensitivity" to manually specify how sensor sensitivity '
            'should be computed.'
        ],
    'PLOT_SHOW':
        ['"PLOT_SHOW" has been renamed to "plot_show".'],
    'PLOT_SAVE':
        ['"PLOT_SAVE" has been renamed to "plot_save".'],
    'PLOT_SAVE_FORMAT':
        ['"PLOT_SAVE_FORMAT" has been renamed to "plot_save_format".'],
    'vp':
        ['"vp" has been renamed to "vp_source".'],
    'vs':
        ['"vs" has been renamed to "vs_source".'],
    'rho':
        ['"rho" has been renamed to "rho_source".'],
    'pre_p_time':
        ['"pre_p_time" has been renamed to "noise_pre_time".'],
    'pre_s_time':
        ['"pre_s_time" has been renamed to "signal_pre_time".'],
    'rps_from_focal_mechanism':
        [
            '"rps_from_focal_mechanism" has been renamed to '
            '"rp_from_focal_mechanism".'
        ],
    'paz':
        ['"paz" has been removed and merged with "station_metadata".'],
    'max_epi_dist':
        ['"max_epi_dist" has been removed and replaced by "epi_dist_ranges".'],
    'max_freq_Er':
        ['"max_freq_Er" has been removed and replaced by "Er_freq_min_max".'],
}
