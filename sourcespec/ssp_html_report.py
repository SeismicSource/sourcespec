
# -*- coding: utf-8 -*-
"""
Generate an HTML report for source_spec.

:copyright:
    2021 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement, Version 2.1
    (http://www.cecill.info/index.en.html)
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import os
import logging
import shutil
import re
import numpy as np
from glob import glob
from sourcespec._version import get_versions
logger = logging.getLogger(__name__.split('.')[-1])


def _multireplace(string, replacements, ignore_case=False):
    """
    Given a string and a replacement map, it returns the replaced string.

    :param str string: string to execute replacements on
    :param dict replacements: replacement dictionary
                              {value to find: value to replace}
    :param bool ignore_case: whether the match should be case insensitive
    :rtype: str

    Source: https://gist.github.com/bgusach/a967e0587d6e01e889fd1d776c5f3729
    """
    if not replacements:
        # Edge case that'd produce a funny regex and cause a KeyError
        return string

    # If case insensitive, we need to normalize the old string so that later a
    # replacement can be found. For instance with {"HEY": "lol"} we should
    # match and find a replacement for "hey",
    # "HEY", "hEy", etc.
    if ignore_case:
        def normalize_old(s):
            return s.lower()

        re_mode = re.IGNORECASE

    else:
        def normalize_old(s):
            return s

        re_mode = 0

    replacements = {
        normalize_old(key): val for key, val in replacements.items()
    }

    # Place longer ones first to keep shorter substrings from matching where
    # the longer ones should take place For instance given the replacements
    # {'ab': 'AB', 'abc': 'ABC'} against the string 'hey abc', it should
    # produce 'hey ABC' and not 'hey ABc'
    rep_sorted = sorted(replacements, key=len, reverse=True)
    rep_escaped = map(re.escape, rep_sorted)

    # Create a big OR regex that matches any of the substrings to replace
    pattern = re.compile("|".join(rep_escaped), re_mode)

    # For each match, look up the new string in the replacements, being the key
    # the normalized old string
    return pattern.sub(
        lambda match: replacements[normalize_old(match.group(0))], string)


def _format_exponent(value, reference):
    """Format `value` to a string having the same exponent than `reference`."""
    # get the exponent of reference value
    xp = np.int(np.floor(np.log10(np.abs(reference))))
    # format value to print it with the same exponent of reference value
    return '{:5.3f}e{:+03d}'.format(value/10**xp, xp)


def html_report(config, sourcepar, sourcepar_err):
    """Generate an HTML report."""
    # Read template files
    template_dir = os.path.join(
        os.path.dirname(os.path.realpath(__file__)),
        'html_report_template'
    )
    style_css = os.path.join(template_dir, 'style.css')
    index_html = os.path.join(template_dir, 'index.html')
    index_html_out = os.path.join(config.options.outdir, 'index.html')
    traces_plot_html = os.path.join(template_dir, 'traces_plot.html')
    spectra_plot_html = os.path.join(template_dir, 'spectra_plot.html')
    station_table_row_html = os.path.join(
        template_dir, 'station_table_row.html')

    # Copy CSS to output dir
    shutil.copy(style_css, config.options.outdir)

    # Version and run completed
    ssp_version = get_versions()['version']
    run_completed = '{} {}'.format(
        config.end_of_run.strftime('%Y-%m-%d %H:%M:%S'),
        config.end_of_run_tz
    )

    # Output files and maps
    hypo = config.hypo
    evid = hypo.evid
    config_file = '{}.ssp.conf'.format(evid)
    out_file = '{}.ssp.out'.format(evid)
    log_file = '{}.ssp.log'.format(evid)
    map_mag = '{}.map_mag.png'.format(evid)
    map_fc = '{}.map_fc.png'.format(evid)

    # Trace plot files
    traces_plot = open(traces_plot_html).read()
    traces_plot_files = glob(
        os.path.join(config.options.outdir, '*.traces.png'))
    traces_plot_files += glob(
        os.path.join(config.options.outdir, '*.traces.??.png'))
    traces_plots = ''
    traces_plot_id = ''
    for traces_plot_file in sorted(traces_plot_files):
        traces_plot_file = os.path.basename(traces_plot_file)
        traces_plots += traces_plot.replace(
            '{TRACES_PLOT_ID}', traces_plot_id).replace(
            '{TRACES_PLOT_FILE}', traces_plot_file)
        traces_plot_id = ' id="print"'

    # Spectral plot files
    spectra_plot = open(spectra_plot_html).read()
    spectra_plot_files = glob(
        os.path.join(config.options.outdir, '*.ssp.png'))
    spectra_plot_files += glob(
        os.path.join(config.options.outdir, '*.ssp.??.png'))
    spectra_plots = ''
    spectra_plot_id = ''
    for spectra_plot_file in sorted(spectra_plot_files):
        spectra_plot_file = os.path.basename(spectra_plot_file)
        spectra_plots += spectra_plot.replace(
            '{SPECTRA_PLOT_ID}', spectra_plot_id).replace(
            '{SPECTRA_PLOT_FILE}', spectra_plot_file)
        spectra_plot_id = ' id="print"'

    # Station table
    station_table_row = open(station_table_row_html).read()
    station_table_rows = ''
    for statId in sorted(sourcepar.keys()):
        if statId in ['means', 'errors', 'means_weight', 'errors_weight']:
            continue
        par = sourcepar[statId]
        err = sourcepar_err[statId]
        id, type = statId.split()
        Mw = par['Mw']
        Mw_err = err['Mw']
        fc = par['fc']
        fc_err = err['fc']
        t_star = par['t_star']
        t_star_err = err['t_star']
        Mo = par['Mo']
        hyp_dist = par['hyp_dist']
        az = par['az']
        Er = par['Er']
        replacements = {
            '{STATION_ID}': id,
            '{STATION_TYPE}': type,
            '{STATION_MW}': '{:.3f}'.format(Mw),
            '{STATION_MW_ERR}': '{:.3f}'.format(Mw_err),
            '{STATION_FC}': '{:.3f}'.format(fc),
            '{STATION_FC_ERR}': '{:.3f}'.format(fc_err),
            '{STATION_TSTAR}': '{:.3f}'.format(t_star),
            '{STATION_TSTAR_ERR}': '{:.3f}'.format(t_star_err),
            '{STATION_M0}': '{:.3e}'.format(Mo),
            '{STATION_DIST}': '{:.3f}'.format(hyp_dist),
            '{STATION_AZ}': '{:.3f}'.format(az),
            '{STATION_ER}': '{:.3e}'.format(Er),
        }
        station_table_rows += _multireplace(station_table_row, replacements)

    # Main HTML page
    means = sourcepar['means']
    errors = sourcepar['errors']
    means_weight = sourcepar['means_weight']
    errors_weight = sourcepar['errors_weight']
    Mw_mean = means['Mw']
    Mw_error = errors['Mw']
    Mw_mean_weight = means_weight['Mw']
    Mw_error_weight = errors_weight['Mw']
    Mo_mean = means['Mo']
    Mo_minus, Mo_plus = errors['Mo']
    Mo_mean_weight = means_weight['Mo']
    Mo_minus_weight, Mo_plus_weight = errors_weight['Mo']
    fc_mean = means['fc']
    fc_minus, fc_plus = errors['fc']
    fc_mean_weight = means_weight['fc']
    fc_minus_weight, fc_plus_weight = errors_weight['fc']
    t_star_mean = means['t_star']
    t_star_error = errors['t_star']
    t_star_mean_weight = means_weight['t_star']
    t_star_error_weight = errors_weight['t_star']
    ra_mean = means['ra']
    ra_minus, ra_plus = errors['ra']
    bsd_mean = means['bsd']
    bsd_minus, bsd_plus = errors['bsd']
    Er_mean = means['Er']
    Er_minus, Er_plus = errors['Er']
    replacements = {
        '{VERSION}': ssp_version,
        '{RUN_COMPLETED}': run_completed,
        '{EVENTID}': evid,
        '{EVENT_LONGITUDE}': '{:8.3f}'.format(hypo.longitude),
        '{EVENT_LATITUDE}': '{:7.3f}'.format(hypo.latitude),
        '{EVENT_DEPTH}': '{:5.1f}'.format(hypo.depth),
        '{ORIGIN_TIME}': '{}'.format(hypo.origin_time),
        '{MW}': '{:.2f}'.format(Mw_mean),
        '{MW_ERR}': '{:.2f}'.format(Mw_error),
        '{MW_WEIGHT}': '{:.2f}'.format(Mw_mean_weight),
        '{MW_WEIGHT_ERR}': '{:.2f}'.format(Mw_error_weight),
        '{M0}': '{:.3e}'.format(Mo_mean),
        '{M0_ERR_MINUS}': '{}'.format(_format_exponent(Mo_minus, Mo_mean)),
        '{M0_ERR_PLUS}': '{}'.format(_format_exponent(Mo_plus, Mo_mean)),
        '{M0_WEIGHT}': '{:.3e}'.format(Mo_mean_weight),
        '{M0_WEIGHT_ERR_MINUS}': '{}'.format(
            _format_exponent(Mo_minus_weight, Mo_mean_weight)),
        '{M0_WEIGHT_ERR_PLUS}': '{}'.format(
            _format_exponent(Mo_plus_weight, Mo_mean_weight)),
        '{FC}': '{:.3f}'.format(fc_mean),
        '{FC_ERR_MINUS}': '{:.3f}'.format(fc_minus),
        '{FC_ERR_PLUS}': '{:.3f}'.format(fc_plus),
        '{FC_WEIGHT}': '{:.3f}'.format(fc_mean_weight),
        '{FC_WEIGHT_ERR_MINUS}': '{:.3f}'.format(fc_minus_weight),
        '{FC_WEIGHT_ERR_PLUS}': '{:.3f}'.format(fc_plus_weight),
        '{TSTAR}': '{:.3f}'.format(t_star_mean),
        '{TSTAR_ERR}': '{:.3f}'.format(t_star_error),
        '{TSTAR_WEIGHT}': '{:.3f}'.format(t_star_mean_weight),
        '{TSTAR_WEIGHT_ERR}': '{:.3f}'.format(t_star_error_weight),
        '{RADIUS}': '{:.3f}'.format(ra_mean),
        '{RADIUS_ERR_MINUS}': '{:.3f}'.format(ra_minus),
        '{RADIUS_ERR_PLUS}': '{:.3f}'.format(ra_plus),
        '{BSD}': '{:.3e}'.format(bsd_mean),
        '{BSD_ERR_MINUS}': '{}'.format(_format_exponent(bsd_minus, bsd_mean)),
        '{BSD_ERR_PLUS}': '{}'.format(_format_exponent(bsd_plus, bsd_mean)),
        '{ER}': '{:.3e}'.format(Er_mean),
        '{ER_ERR_MINUS}': '{}'.format(_format_exponent(Er_minus, Er_mean)),
        '{ER_ERR_PLUS}': '{}'.format(_format_exponent(Er_plus, Er_mean)),
        '{CONF_FILE_BNAME}': config_file,
        '{CONF_FILE}': config_file,
        '{OUT_FILE_BNAME}': out_file,
        '{OUT_FILE}': out_file,
        '{LOG_FILE_BNAME}': log_file,
        '{LOG_FILE}': log_file,
        '{MAP_MAG}': map_mag,
        '{MAP_FC}': map_fc,
        '{TRACES_PLOTS}': traces_plots,
        '{SPECTRA_PLOTS}': spectra_plots,
        '{STATION_TABLE_ROWS}': station_table_rows,
    }
    index = open(index_html).read()
    index = _multireplace(index, replacements)
    with open(index_html_out, 'w') as fp:
        fp.write(index)
    logger.info('HTML report written to file: ' + index_html_out)
