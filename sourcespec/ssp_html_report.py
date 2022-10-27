# -*- coding: utf-8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Generate an HTML report for source_spec.

:copyright:
    2021-2022 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
import os
import logging
import shutil
import re
import numpy as np
from urllib.parse import urlparse
from sourcespec._version import get_versions
from sourcespec.ssp_data_types import SpectralParameter
from sourcespec.ssp_setup import ssp_exit
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


def _agency_logo(config):
    agency_logo = config.agency_logo
    if config.agency_logo is None:
        return ""
    # check if agency_logo is a URL
    parsed_url = urlparse(agency_logo)
    if not parsed_url.scheme:
        # check if it is a file
        if not os.path.exists(agency_logo):
            logger.error('Cannot find the agency logo file: {}'.format(
                agency_logo))
            ssp_exit(1)
        shutil.copy(agency_logo, config.options.outdir)
        agency_logo = os.path.basename(agency_logo)
    agency_logo_img = '<img class="logo" src="{}"/>'.format(agency_logo)
    indent5 = 5*'  '
    indent6 = 6*'  '
    if config.agency_url is not None:
        agency_logo_html =\
            '{}<a href="{}" target="_blank">\n{}{}\n{}</a>'.format(
                indent5, config.agency_url,
                indent6, agency_logo_img,
                indent5
            )
    else:
        agency_logo_html = indent5 + agency_logo_img
    agency_logo_html = agency_logo_html +\
        '\n{}<hr class="solid">'.format(indent5)
    return agency_logo_html


def _logo_file_url():
    cdn_baseurl = 'https://cdn.jsdelivr.net/gh/SeismicSource/sourcespec@1.6'
    logo_file = cdn_baseurl + '/imgs/SourceSpec_logo.svg'
    return logo_file


def _version_and_run_completed(config):
    ssp_version = get_versions()['version']
    run_completed = '{} {}'.format(
        config.end_of_run.strftime('%Y-%m-%d %H:%M:%S'),
        config.end_of_run_tz
    )
    return ssp_version, run_completed


def _author_and_agency(config):
    author = ''
    if config.author_name is not None:
        author = config.author_name
    elif config.author_email is not None:
        author = config.author_email
    if config.author_email is not None:
        author = '<a href="mailto:{}">{}</a>'.format(
            config.author_email, author)
    agency = ''
    if config.agency_full_name is not None:
        agency = config.agency_full_name
        if config.agency_short_name is not None:
            agency += ' ({})'.format(config.agency_short_name)
    elif config.agency_short_name is not None:
        agency = config.agency_short_name
    elif config.agency_url is not None:
        agency = config.agency_url
    if config.agency_url is not None:
        agency = '<a href="{}" target="_blank">{}</a>'.format(
            config.agency_url, agency)
    if author != '':
        author = '<br/><br/>' + author
    if author == '' and agency != '':
        agency = '<br/><br/>' + agency
    if author != '' and agency != '':
        agency = '<br/>' + agency
    return author, agency


def _format_exponent(value, reference):
    """Format `value` to a string having the same exponent than `reference`."""
    # get the exponent of reference value
    xp = np.int(np.floor(np.log10(np.abs(reference))))
    # format value to print it with the same exponent of reference value
    return '{:5.3f}e{:+03d}'.format(value/10**xp, xp)


def _summary_value_and_err_text(value, error, fmt):
    """Format summary value and error text."""
    if error[0] == error[1]:
        text = fmt + '<br>&#177;' + fmt
        text = text.format(value, error[0])
    else:
        text = fmt + '<br>-' + fmt + '<br>+' + fmt
        text = text.format(value, error[0], error[1])
    return text


def _station_value_and_err_text(par, key, fmt):
    """Format station value and error text."""
    try:
        _par = par[key]
    except KeyError:
        return '', ''
    if isinstance(_par, SpectralParameter):
        value = _par.value
        outlier = _par.outlier
    else:
        value = _par
        outlier = False
    value_text = '<nobr>{}</nobr>'.format(fmt.format(value))
    err_text = None
    if isinstance(_par, SpectralParameter):
        if _par.uncertainty is not None:
            # use HTML code for ±, for compatibility with Edge
            err_text = '<nobr>&#177;{}</nobr>'.format(
                fmt.format(_par.uncertainty))
        elif _par.lower_uncertainty is not None:
            err_text = '<nobr>-{}</nobr><br/><nobr>+{}</nobr>'.format(
                fmt.format(_par.lower_uncertainty),
                fmt.format(_par.upper_uncertainty)
            )
    if outlier:
        value_text = '<span style="color:#979A9A">' + value_text + '</span>'
        if err_text is not None:
            err_text = '<span style="color:#979A9A">' + err_text + '</span>'
    return value_text, err_text


def _misfit_table_rows(misfit_plot_files):
    template_dir = os.path.join(
        os.path.dirname(os.path.realpath(__file__)),
        'html_report_template'
    )
    misfit_table_column_html = os.path.join(
        template_dir, 'misfit_table_column.html')
    misfit_table_column = open(misfit_table_column_html).read()
    misfit_table_rows = ''
    for n, misfit_plot_file in enumerate(sorted(misfit_plot_files)):
        if n % 3 == 0:
            misfit_table_rows += '<tr>\n'
        misfit_plot_file = os.path.join(
            'misfit', os.path.basename(misfit_plot_file))
        misfit_table_rows += misfit_table_column.replace(
            '{MISFIT_PLOT}', misfit_plot_file)
        misfit_table_rows += '\n'
        if n % 3 == 2:
            misfit_table_rows += 10*' '
            misfit_table_rows += '</tr>\n'
            misfit_table_rows += 10*' '
    misfit_table_rows += 10*' '
    misfit_table_rows += '</tr>'
    return misfit_table_rows


def _misfit_page(config):
    """Generate an HTML page with misfit plots."""
    # Read template files
    template_dir = os.path.join(
        os.path.dirname(os.path.realpath(__file__)),
        'html_report_template'
    )
    misfit_html = os.path.join(template_dir, 'misfit.html')
    misfit_html_out = os.path.join(config.options.outdir, 'misfit.html')

    # Logo file
    logo_file = _logo_file_url()

    # Version and run completed
    ssp_version, run_completed = _version_and_run_completed(config)

    # Author and agency
    author, agency = _author_and_agency(config)

    # event info
    hypo = config.hypo
    evid = hypo.evid

    # 1d conditional misfit plots
    misfit_plot_files = config.figures['misfit_1d']
    one_d_misfit_table_rows = _misfit_table_rows(misfit_plot_files)

    # 2d conditional misfit plots: fc-Mw
    misfit_plot_files = config.figures['misfit_fc-Mw']
    two_d_misfit_table_rows_fc_mw = _misfit_table_rows(misfit_plot_files)

    # 2d conditional misfit plots: fc-tstar
    misfit_plot_files = config.figures['misfit_fc-t_star']
    two_d_misfit_table_rows_fc_tstar = _misfit_table_rows(misfit_plot_files)

    # 2d conditional misfit plots: tstar-Mw
    misfit_plot_files = config.figures['misfit_t_star-Mw']
    two_d_misfit_table_rows_tstar_mw = _misfit_table_rows(misfit_plot_files)

    # Main HTML page
    replacements = {
        '{LOGO_FILE}': logo_file,
        '{VERSION}': ssp_version,
        '{RUN_COMPLETED}': run_completed,
        '{AUTHOR}': author,
        '{AGENCY}': agency,
        '{EVENTID}': evid,
        '{1D_MISFIT_TABLE_ROWS}': one_d_misfit_table_rows,
        '{2D_MISFIT_TABLE_ROWS_FC_MW}': two_d_misfit_table_rows_fc_mw,
        '{2D_MISFIT_TABLE_ROWS_FC_TSTAR}': two_d_misfit_table_rows_fc_tstar,
        '{2D_MISFIT_TABLE_ROWS_TSTAR_MW}': two_d_misfit_table_rows_tstar_mw
    }
    misfit = open(misfit_html).read()
    misfit = _multireplace(misfit, replacements)
    with open(misfit_html_out, 'w') as fp:
        fp.write(misfit)


def html_report(config, sspec_output):
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
    quakeml_file_link_html = os.path.join(
        template_dir, 'quakeml_file_link.html')

    # Copy CSS to output dir
    shutil.copy(style_css, config.options.outdir)

    # Logo file
    logo_file = _logo_file_url()

    # HTML for agency logo
    agency_logo = _agency_logo(config)

    # Version and run completed
    ssp_version, run_completed = _version_and_run_completed(config)

    # Author and agency
    author, agency = _author_and_agency(config)

    # Output files and maps
    hypo = config.hypo
    evid = hypo.evid
    run_id = config.options.run_id
    config_file = '{}.ssp.conf'.format(evid)
    yaml_file = '{}.ssp.yaml'.format(evid)
    log_file = '{}.ssp.log'.format(evid)
    station_maps = config.figures['station_maps']
    map_mag = [mapfile for mapfile in station_maps if 'map_mag' in mapfile][0]
    map_mag = os.path.basename(map_mag)
    map_fc = [mapfile for mapfile in station_maps if 'map_fc' in mapfile][0]
    map_fc = os.path.basename(map_fc)
    box_plots = config.figures['boxplots'][0]
    box_plots = os.path.basename(box_plots)

    # Trace plot files
    traces_plot = open(traces_plot_html).read()
    traces_plot_files = config.figures['traces']
    traces_plots = ''
    traces_plot_class = ''
    n_traces_plot_files = len(traces_plot_files)
    for n, traces_plot_file in enumerate(sorted(traces_plot_files)):
        if n_traces_plot_files > 1:
            traces_plot_counter =\
                '<span class="print_inline">&nbsp;({} of {})</span>'.format(
                    n+1, n_traces_plot_files)
        else:
            traces_plot_counter = ''
        traces_plot_file = os.path.basename(traces_plot_file)
        traces_plots += traces_plot.\
            replace('{TRACES_PLOT_CLASS}', traces_plot_class).\
            replace('{TRACES_PLOT_COUNTER}', traces_plot_counter).\
            replace('{TRACES_PLOT_FILE}', traces_plot_file)
        traces_plot_class = ' class="print"'

    # Spectral plot files
    spectra_plot = open(spectra_plot_html).read()
    spectra_plot_files = config.figures['spectra_regular']
    spectra_plots = ''
    spectra_plot_class = ''
    n_spectra_plot_files = len(spectra_plot_files)
    for n, spectra_plot_file in enumerate(sorted(spectra_plot_files)):
        if n_spectra_plot_files > 1:
            spectra_plot_counter =\
                '<span class="print_inline">&nbsp;({} of {})</span>'.format(
                    n+1, n_spectra_plot_files)
        else:
            spectra_plot_counter = ''
        spectra_plot_file = os.path.basename(spectra_plot_file)
        spectra_plots += spectra_plot.\
            replace('{SPECTRA_PLOT_CLASS}', spectra_plot_class).\
            replace('{SPECTRA_PLOT_COUNTER}', spectra_plot_counter).\
            replace('{SPECTRA_PLOT_FILE}', spectra_plot_file)
        spectra_plot_class = ' class="print"'

    # Station table
    station_table_row = open(station_table_row_html).read()
    station_table_rows = ''
    stationpar = sspec_output.station_parameters
    for statId in sorted(stationpar.keys()):
        par = stationpar[statId]
        instrument_type = par.instrument_type
        Mw_text, Mw_err_text = _station_value_and_err_text(par, 'Mw', '{:.3f}')
        fc_text, fc_err_text = _station_value_and_err_text(par, 'fc', '{:.3f}')
        t_star_text, t_star_err_text =\
            _station_value_and_err_text(par, 't_star', '{:.3f}')
        Qo_text, Qo_err_text = _station_value_and_err_text(par, 'Qo', '{:.1f}')
        Mo_text, Mo_err_text = _station_value_and_err_text(par, 'Mo', '{:.3e}')
        bsd_text, bsd_err_text =\
            _station_value_and_err_text(par, 'bsd', '{:.3e}')
        ra_text, ra_err_text =\
            _station_value_and_err_text(par, 'radius', '{:.3f}')
        Er_text, _ = _station_value_and_err_text(par, 'Er', '{:.3e}')
        Ml_text, _ = _station_value_and_err_text(par, 'Ml', '{:.3f}')
        hyp_dist_text, _ =\
            _station_value_and_err_text(par, 'hypo_dist_in_km', '{:.3f}')
        az_text, _ = _station_value_and_err_text(par, 'azimuth', '{:.3f}')
        replacements = {
            '{STATION_ID}': statId,
            '{STATION_TYPE}': instrument_type,
            '{STATION_MW}': Mw_text,
            '{STATION_MW_ERR}': Mw_err_text,
            '{STATION_FC}': fc_text,
            '{STATION_FC_ERR}': fc_err_text,
            '{STATION_TSTAR}': t_star_text,
            '{STATION_TSTAR_ERR}': t_star_err_text,
            '{STATION_Q0}': Qo_text,
            '{STATION_Q0_ERR}': Qo_err_text,
            '{STATION_M0}': Mo_text,
            '{STATION_M0_ERR}': Mo_err_text,
            '{STATION_BSD}': bsd_text,
            '{STATION_BSD_ERR}': bsd_err_text,
            '{STATION_RA}': ra_text,
            '{STATION_RA_ERR}': ra_err_text,
            '{STATION_ER}': Er_text,
            '{STATION_ML}': Ml_text,
            '{STATION_DIST}': hyp_dist_text,
            '{STATION_AZ}': az_text,
        }
        # Local magnitude, if computed
        if config.compute_local_magnitude:
            Ml_comment_begin = ''
            Ml_comment_end = ''
        else:
            Ml_comment_begin = '<!--'
            Ml_comment_end = '-->'
        replacements.update({
            '{ML_COMMENT_BEGIN}': Ml_comment_begin,
            '{ML_COMMENT_END}': Ml_comment_end
        })
        station_table_rows += _multireplace(station_table_row, replacements)

    # Event and run info
    replacements = {
        '{AGENCY_LOGO}': agency_logo,
        '{LOGO_FILE}': logo_file,
        '{VERSION}': ssp_version,
        '{RUN_COMPLETED}': run_completed,
        '{AUTHOR}': author,
        '{AGENCY}': agency,
        '{EVENTID}': evid,
        '{RUNID}': run_id,
        '{EVENT_LONGITUDE}': '{:8.3f}'.format(hypo.longitude),
        '{EVENT_LATITUDE}': '{:7.3f}'.format(hypo.latitude),
        '{EVENT_DEPTH}': '{:5.1f}'.format(hypo.depth),
        '{ORIGIN_TIME}': '{}'.format(hypo.origin_time),
    }
    # Link to event page, if defined
    event_url = config.event_url
    if event_url is not None:
        event_url = event_url.replace('$EVENTID', evid)
        parsed_url = urlparse(event_url)
        if not parsed_url.scheme:
            logger.warning(
                '{} is not a valid URL and will not be used'.format(event_url))
            event_url = None
    if event_url is not None:
        event_url_comment_begin = ''
        event_url_comment_end = ''
    else:
        event_url = ''
        event_url_comment_begin = '<!--'
        event_url_comment_end = '-->'
    replacements.update({
        '{EVENT_URL}': event_url,
        '{EVENT_URL_COMMENT_BEGIN}': event_url_comment_begin,
        '{EVENT_URL_COMMENT_END}': event_url_comment_end
    })
    # Only show Run ID if it is not empty
    if run_id:
        run_id_comment_begin = ''
        run_id_comment_end = ''
    else:
        run_id_comment_begin = '<!--'
        run_id_comment_end = '-->'
    replacements.update({
        '{RUNID_COMMENT_BEGIN}': run_id_comment_begin,
        '{RUNID_COMMENT_END}': run_id_comment_end
    })

    # Summary spectral parameters
    ref_stat = sspec_output.summary_spectral_parameters.reference_statistics
    col_mean_highlighted = col_wmean_highlighted = col_perc_highlighted = ''
    if ref_stat == 'mean':
        col_mean_highlighted = 'class="highlighted_column"'
    elif ref_stat == 'weighted_mean':
        col_wmean_highlighted = 'class="highlighted_column"'
    elif ref_stat == 'percentiles':
        col_perc_highlighted = 'class="highlighted_column"'
    replacements.update({
        '{COL_MEAN_HIGHLIGHTED}': col_mean_highlighted,
        '{COL_WMEAN_HIGHLIGHTED}': col_wmean_highlighted,
        '{COL_PERC_HIGHLIGHTED}': col_perc_highlighted,
    })

    summary_values = sspec_output.reference_values()
    summary_uncertainties = sspec_output.reference_uncertainties()
    Mw_summary = summary_values['Mw']
    try:
        Mw_summary_error_minus, Mw_summary_error_plus =\
            summary_uncertainties['Mw']
        Mw_summary_str = '{:.2f} [- {:.2f}, + {:.2f}]'.format(
            Mw_summary, Mw_summary_error_minus, Mw_summary_error_plus)
    except TypeError:
        Mw_summary_error = summary_uncertainties['Mw']
        Mw_summary_str = '{:.2f} ± {:.2f}'.format(Mw_summary, Mw_summary_error)
    replacements.update({'{MW_SUMMARY}': Mw_summary_str})
    fc_summary = summary_values['fc']
    try:
        fc_summary_error_minus, fc_summary_error_plus =\
            summary_uncertainties['fc']
        fc_summary_str = '{:.3f} [- {:.3f}, + {:.3f}]'.format(
            fc_summary, fc_summary_error_minus, fc_summary_error_plus)
    except TypeError:
        fc_summary_error = summary_uncertainties['fc']
        fc_summary_str = '{:.3f} ± {:.3f}'.format(fc_summary, fc_summary_error)
    replacements.update({'{FC_SUMMARY}': fc_summary_str})

    means = sspec_output.mean_values()
    mean_errors = sspec_output.mean_uncertainties()
    wmeans = sspec_output.weighted_mean_values()
    wmean_errors = sspec_output.weighted_mean_uncertainties()
    percentiles = sspec_output.percentiles_values()
    percentile_errors = sspec_output.percentiles_uncertainties()

    n_sigma = config.n_sigma
    n_sigma = int(n_sigma) if n_sigma.is_integer() else n_sigma
    n_sigma = '{} sigma'.format(n_sigma)
    mid_pct, lower_pct, upper_pct =\
        config.mid_percentage, config.lower_percentage, config.upper_percentage
    mid_pct = int(mid_pct) if mid_pct.is_integer() else mid_pct
    lower_pct = int(lower_pct) if lower_pct.is_integer() else lower_pct
    upper_pct = int(upper_pct) if upper_pct.is_integer() else upper_pct
    percentages = '{}%, [{}%, {}%]'.format(mid_pct, lower_pct, upper_pct)
    replacements.update({
        '{N_SIGMA}': n_sigma,
        '{PERCENTAGES}': percentages
    })

    Mw_mean = means['Mw']
    Mw_mean_error = mean_errors['Mw']
    Mw_wmean = wmeans['Mw']
    Mw_wmean_error = wmean_errors['Mw']
    Mw_perc = percentiles['Mw']
    Mw_perc_error = percentile_errors['Mw']
    replacements.update({
        '{MW_MEAN_AND_ERR}': _summary_value_and_err_text(
            Mw_mean, Mw_mean_error, '{:.2f}'),
        '{MW_WMEAN_AND_ERR}': _summary_value_and_err_text(
            Mw_wmean, Mw_wmean_error, '{:.2f}'),
        '{MW_PERC_AND_ERR}': _summary_value_and_err_text(
            Mw_perc, Mw_perc_error, '{:.2f}'),
    })
    Mo_mean = means['Mo']
    Mo_mean_error = mean_errors['Mo']
    Mo_wmean = wmeans['Mo']
    Mo_wmean_error = wmean_errors['Mo']
    Mo_perc = percentiles['Mo']
    Mo_perc_error = percentile_errors['Mo']
    replacements.update({
        '{M0_MEAN_AND_ERR}': _summary_value_and_err_text(
            Mo_mean, Mo_mean_error, '{:.3e}'),
        '{M0_WMEAN_AND_ERR}': _summary_value_and_err_text(
            Mo_wmean, Mo_wmean_error, '{:.3e}'),
        '{M0_PERC_AND_ERR}': _summary_value_and_err_text(
            Mo_perc, Mo_perc_error, '{:.3e}'),
    })
    fc_mean = means['fc']
    fc_mean_error = mean_errors['fc']
    fc_wmean = wmeans['fc']
    fc_wmean_error = wmean_errors['fc']
    fc_perc = percentiles['fc']
    fc_perc_error = percentile_errors['fc']
    replacements.update({
        '{FC_MEAN_AND_ERR}': _summary_value_and_err_text(
            fc_mean, fc_mean_error, '{:.3f}'),
        '{FC_WMEAN_AND_ERR}': _summary_value_and_err_text(
            fc_wmean, fc_wmean_error, '{:.3f}'),
        '{FC_PERC_AND_ERR}': _summary_value_and_err_text(
            fc_perc, fc_perc_error, '{:.3f}'),
    })
    t_star_mean = means['t_star']
    t_star_mean_error = mean_errors['t_star']
    t_star_wmean = wmeans['t_star']
    t_star_wmean_error = wmean_errors['t_star']
    t_star_perc = percentiles['t_star']
    t_star_perc_error = percentile_errors['t_star']
    replacements.update({
        '{TSTAR_MEAN_AND_ERR}': _summary_value_and_err_text(
            t_star_mean, t_star_mean_error, '{:.3f}'),
        '{TSTAR_WMEAN_AND_ERR}': _summary_value_and_err_text(
            t_star_wmean, t_star_wmean_error, '{:.3f}'),
        '{TSTAR_PERC_AND_ERR}': _summary_value_and_err_text(
            t_star_perc, t_star_perc_error, '{:.3f}'),
    })
    Qo_mean = means['Qo']
    Qo_mean_error = mean_errors['Qo']
    Qo_wmean = wmeans['Qo']
    Qo_wmean_error = wmean_errors['Qo']
    Qo_perc = percentiles['Qo']
    Qo_perc_error = percentile_errors['Qo']
    replacements.update({
        '{Q0_MEAN_AND_ERR}': _summary_value_and_err_text(
            Qo_mean, Qo_mean_error, '{:.1f}'),
        '{Q0_WMEAN_AND_ERR}': _summary_value_and_err_text(
            Qo_wmean, Qo_wmean_error, '{:.1f}'),
        '{Q0_PERC_AND_ERR}': _summary_value_and_err_text(
            Qo_perc, Qo_perc_error, '{:.1f}'),
    })
    ra_mean = means['radius']
    ra_mean_error = mean_errors['radius']
    ra_wmean = wmeans['radius']
    ra_wmean_error = wmean_errors['radius']
    ra_perc = percentiles['radius']
    ra_perc_error = percentile_errors['radius']
    replacements.update({
        '{RADIUS_MEAN_AND_ERR}': _summary_value_and_err_text(
            ra_mean, ra_mean_error, '{:.3f}'),
        '{RADIUS_WMEAN_AND_ERR}': _summary_value_and_err_text(
            ra_wmean, ra_wmean_error, '{:.3f}'),
        '{RADIUS_PERC_AND_ERR}': _summary_value_and_err_text(
            ra_perc, ra_perc_error, '{:.3f}'),
    })
    bsd_mean = means['bsd']
    bsd_mean_error = mean_errors['bsd']
    bsd_wmean = wmeans['bsd']
    bsd_wmean_error = wmean_errors['bsd']
    bsd_perc = percentiles['bsd']
    bsd_perc_error = percentile_errors['bsd']
    replacements.update({
        '{BSD_MEAN_AND_ERR}': _summary_value_and_err_text(
            bsd_mean, bsd_mean_error, '{:.3e}'),
        '{BSD_WMEAN_AND_ERR}': _summary_value_and_err_text(
            bsd_wmean, bsd_wmean_error, '{:.3e}'),
        '{BSD_PERC_AND_ERR}': _summary_value_and_err_text(
            bsd_perc, bsd_perc_error, '{:.3e}'),
    })
    Er_mean = means['Er']
    Er_mean_error = mean_errors['Er']
    Er_perc = percentiles['Er']
    Er_perc_error = percentile_errors['Er']
    replacements.update({
        '{ER_MEAN_AND_ERR}': _summary_value_and_err_text(
            Er_mean, Er_mean_error, '{:.3e}'),
        '{ER_PERC_AND_ERR}': _summary_value_and_err_text(
            Er_perc, Er_perc_error, '{:.3e}'),
    })
    # Local magnitude, if computed
    if config.compute_local_magnitude:
        Ml_mean = means['Ml']
        Ml_mean_error = mean_errors['Ml']
        Ml_perc = percentiles['Ml']
        Ml_perc_error = percentile_errors['Ml']
        Ml_comment_begin = ''
        Ml_comment_end = ''
    else:
        Ml_mean = Ml_perc = np.nan
        Ml_mean_error = Ml_perc_error = (np.nan, np.nan)
        Ml_comment_begin = '<!--'
        Ml_comment_end = '-->'
    replacements.update({
        '{ML_MEAN_AND_ERR}': _summary_value_and_err_text(
            Ml_mean, Ml_mean_error, '{:.2f}'),
        '{ML_PERC_AND_ERR}': _summary_value_and_err_text(
            Ml_perc, Ml_perc_error, '{:.2f}'),
        '{ML_COMMENT_BEGIN}': Ml_comment_begin,
        '{ML_COMMENT_END}': Ml_comment_end,
    })

    # Output files and plots
    replacements.update({
        '{CONF_FILE_BNAME}': config_file,
        '{CONF_FILE}': config_file,
        '{YAML_FILE_BNAME}': yaml_file,
        '{YAML_FILE}': yaml_file,
        '{LOG_FILE_BNAME}': log_file,
        '{LOG_FILE}': log_file,
        '{MAP_MAG}': map_mag,
        '{MAP_FC}': map_fc,
        '{TRACES_PLOTS}': traces_plots,
        '{SPECTRA_PLOTS}': spectra_plots,
        '{BOX_PLOTS}': box_plots,
        '{STATION_TABLE_ROWS}': station_table_rows,
    })
    # Misfit plots (when using grid search)
    if 'misfit_1d' in config.figures:
        misfit_plot_comment_begin = ''
        misfit_plot_comment_end = ''
        _misfit_page(config)
    else:
        misfit_plot_comment_begin = '<!--'
        misfit_plot_comment_end = '-->'
    replacements.update({
        '{MISFIT_PLOT_COMMENT_BEGIN}': misfit_plot_comment_begin,
        '{MISFIT_PLOT_COMMENT_END}': misfit_plot_comment_end
    })
    # QuakeML file (if produced)
    if config.qml_file_out is not None:
        quakeml_file = os.path.basename(config.qml_file_out)
        quakeml_file_link = open(quakeml_file_link_html).read()
        quakeml_file_link = quakeml_file_link\
            .replace('{QUAKEML_FILE}', quakeml_file)\
            .replace('{QUAKEML_FILE_BNAME}', quakeml_file)
    else:
        quakeml_file_link = ''
    replacements.update({
        '{QUAKEML_FILE_LINK}': quakeml_file_link
    })

    index = open(index_html).read()
    index = _multireplace(index, replacements)
    with open(index_html_out, 'w') as fp:
        fp.write(index)
    logger.info('HTML report written to file: ' + index_html_out)
