# -*- coding: utf-8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Generate an HTML report for source_spec.

:copyright:
    2021-2025 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
import os
import logging
import shutil
import re
import contextlib
from urllib.parse import urlparse
import numpy as np
from .setup import config
from ._version import get_versions
from .ssp_data_types import SpectralParameter
logger = logging.getLogger(__name__.rsplit('.', maxsplit=1)[-1])
VALID_FIGURE_FORMATS = ('.png', '.svg')


def _multireplace(string, replacements, ignore_case=False):
    """
    Given a string and a replacement map, it returns the replaced string.

    :param tring: string to execute replacements on
    :type string: str
    :param replacements: replacement dictionary
        {value to find: value to replace}
    :type replacements: dict
    :param ignore_case: whether the match should be case insensitive
    :type ignore_case: bool

    :return: replaced string
    :rtype: str

    .. note::

        Source:
        https://gist.github.com/bgusach/a967e0587d6e01e889fd1d776c5f3729
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


def _agency_logo_path():
    """
    Return the path to the agency logo file.

    :return: path to the agency logo file
    :rtype: str
    """
    agency_logo_path = config.agency_logo
    if agency_logo_path is None:
        return None
    # check if agency_logo is a URL
    parsed_url = urlparse(agency_logo_path)
    if not parsed_url.scheme:
        # check if it is a file
        if not os.path.exists(agency_logo_path):
            logger.warning(
                f'Cannot find the agency logo file: {agency_logo_path}')
            return None
        bname = os.path.basename(agency_logo_path)
        dest = os.path.join(config.options.outdir, bname)
        if not os.path.exists(dest):
            shutil.copy(agency_logo_path, config.options.outdir)
        agency_logo_path = bname
    return agency_logo_path


def _agency_logo():
    """
    Return the HTML code for the agency logo.

    :return: HTML code for the agency logo
    :rtype: str
    """
    agency_logo_path = _agency_logo_path()
    if agency_logo_path is None:
        return ''
    agency_logo_img = f'<img class="logo" src="{agency_logo_path}"/>'
    indent5 = 5 * '  '
    indent6 = 6 * '  '
    if config.agency_url is not None:
        agency_logo_html = (
            f'{indent5}<a href="{config.agency_url}" target="_blank">\n'
            f'{indent6}{agency_logo_img}\n'
            f'{indent5}</a>'
        )
    else:
        agency_logo_html = indent5 + agency_logo_img
    agency_logo_html = f'{agency_logo_html}\n{indent5}<hr class="solid">'
    return agency_logo_html


def _logo_file_url():
    """
    Return the URL of the SourceSpec logo file.

    :return: URL of the SourceSpec logo file
    :rtype: str
    """
    cdn_baseurl = 'https://cdn.jsdelivr.net/gh/SeismicSource/sourcespec@1.6'
    return f'{cdn_baseurl}/imgs/SourceSpec_logo.svg'


def _version_and_run_completed():
    """
    Return the SourceSpec version and the run completion date.

    :return: SourceSpec version and run completion date
    :rtype: tuple of str
    """
    ssp_version = get_versions()['version']
    run_completed = (
        f'{config.end_of_run.strftime("%Y-%m-%d %H:%M:%S")} '
        f'{config.end_of_run_tz}'
    )
    return ssp_version, run_completed


def _author_html():
    """
    Return the HTML code for the author.

    :return: HTML code for the author
    :rtype: str
    """
    author = ''
    if config.author_name is not None:
        author = config.author_name
    elif config.author_email is not None:
        author = config.author_email
    if config.author_email is not None:
        author = f'<a href="mailto:{config.author_email}">{author}</a>'
    return author


def _agency_html():
    """
    Return the HTML code for the agency.

    :return: HTML code for the agency
    :rtype: str
    """
    agency = ''
    if config.agency_full_name is not None:
        agency = config.agency_full_name
        if config.agency_short_name is not None:
            agency += f' ({config.agency_short_name})'
    elif config.agency_short_name is not None:
        agency = config.agency_short_name
    elif config.agency_url is not None:
        agency = config.agency_url
    if config.agency_url is not None:
        agency = f'<a href="{config.agency_url}" target="_blank">{agency}</a>'
    return agency


def _author_and_agency_html(author, agency):
    """
    Return the HTML code for the author and the agency.

    :param author: HTML code for the author
    :type author: str
    :param agency: HTML code for the agency
    :type agency: str

    :return: HTML code for the author and the agency
    :rtype: str
    """
    if author != '':
        author = f'<br/><br/>{author}'
    if author == '' and agency != '':
        agency = f'<br/><br/>{agency}'
    if author != '' and agency != '':
        agency = f'<br/>{agency}'
    return author + agency


def _page_footer():
    """
    Return the HTML code for the page footer.

    :return: HTML code for the page footer
    :rtype: str
    """
    footer_html = ''
    indent3 = 3 * '  '
    indent4 = 4 * '  '
    footer_html += f'{indent3}<div class="text_footer">\n'
    author = _author_html()
    agency = _agency_html()
    auth_agen_text = ''
    if author != '':
        if agency != '':
            auth_agen_text +=\
                f'{indent4}{author}\n{indent4}-\n{indent4}{agency}\n'
        if agency == '':
            auth_agen_text += f'{indent4}{author}\n'
    if author == '' and agency != '':
        auth_agen_text += f'{indent4}{agency}\n'
    run_completed = config.end_of_run.strftime('%Y-%m-%d')
    footer_html += (
        f'{auth_agen_text}\n{indent4}-\n{indent4}{run_completed}\n'
        if auth_agen_text
        else f'{indent4}{run_completed}\n'
    )
    footer_html += f'{indent3}</div>\n'
    agency_logo_path = _agency_logo_path()
    if agency_logo_path is not None:
        if config.agency_url is not None:
            a_agency = f'<a href="{config.agency_url}" target="_blank">'
            a_agency_close = '</a>'
        else:
            a_agency = a_agency_close = ''
        footer_html += (
            f'{indent3}<div class="logo_footer">\n{indent4}{a_agency}'
            f'<img src="{agency_logo_path}"/>{a_agency_close}\n'
            f'{indent3}</div>\n'
        )
    return footer_html


def _format_exponent(value, reference):
    """
    Format `value` to a string having the same exponent than `reference`.

    :param value: value to format
    :type value: float
    :param reference: reference value
    :type reference: float

    :return: formatted value
    :rtype: str
    """
    # get the exponent of reference value
    xp = np.int(np.floor(np.log10(np.abs(reference))))
    # format value to print it with the same exponent of reference value
    base = 10**xp
    return f'{value/base:5.3f}e{xp:+03d}'


def _summary_value_and_err_text(value, error, fmt):
    """
    Format summary value and error text.

    :param value: value
    :type value: float
    :param error: error
    :type error: float or tuple of float
    :param fmt: format string
    :type fmt: str

    :return: formatted text
    :rtype: str
    """
    if error[0] == error[1]:
        text = f'{fmt}<br>&#177;{fmt}'
        text = text.format(value, error[0])
    else:
        text = f'{fmt}<br>-{fmt}<br>+{fmt}'
        text = text.format(value, error[0], error[1])
    return text


def _station_value_and_err_text(par, key, fmt):
    """
    Format station value and error text.

    :param par: dictionary of station parameters
    :type par: dict
    :param key: parameter key
    :type key: str
    :param fmt: format string
    :type fmt: str

    :return: formatted value and error text
    :rtype: tuple of str
    """
    # par[key] can be None even if key is in par (e.g., if key is 'Ml' and
    # local magnitude is not computed)
    _par = par[key] if key in par else None
    if _par is None:
        return '', ''
    if isinstance(_par, SpectralParameter):
        value = _par.value
        outlier = _par.outlier
    else:
        value = _par
        outlier = False
    value_text = f'<nobr>{fmt.format(value)}</nobr>'
    err_text = None
    if isinstance(_par, SpectralParameter):
        if _par.uncertainty is not None:
            # use HTML code for ±, for compatibility with Edge
            err_text = f'<nobr>&#177;{fmt.format(_par.uncertainty)}</nobr>'
        elif _par.lower_uncertainty is not None:
            err_text = (
                f'<nobr>-{fmt.format(_par.lower_uncertainty)}</nobr><br/>'
                f'<nobr>+{fmt.format(_par.upper_uncertainty)}</nobr>'
            )
    if outlier:
        value_text = f'<span style="color:#979A9A">{value_text}</span>'
        if err_text is not None:
            err_text = f'<span style="color:#979A9A">{err_text}</span>'
    return value_text, err_text


def _misfit_table_rows(misfit_plot_files):
    """
    Return the HTML code for the misfit table rows.

    :param misfit_plot_files: list of misfit plot files
    :type misfit_plot_files: list of str

    :return: HTML code for the misfit table rows
    :rtype: str
    """
    template_dir = os.path.join(
        os.path.dirname(os.path.realpath(__file__)),
        'html_report_template'
    )
    misfit_table_column_html = os.path.join(
        template_dir, 'misfit_table_column.html')
    with open(
        misfit_table_column_html, encoding='utf-8'
    ) as fp:
        misfit_table_column = fp.read()
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
            misfit_table_rows += 10 * ' '
            misfit_table_rows += '</tr>\n'
            misfit_table_rows += 10 * ' '
    misfit_table_rows += 10 * ' '
    misfit_table_rows += '</tr>'
    return misfit_table_rows


def _misfit_page():
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
    ssp_version, run_completed = _version_and_run_completed()

    # Author and agency
    author = _author_html()
    agency = _agency_html()
    author_and_agency = _author_and_agency_html(author, agency)

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
        '{AUTHOR_AND_AGENCY}': author_and_agency,
        '{EVENTID}': config.event.event_id,
        '{1D_MISFIT_TABLE_ROWS}': one_d_misfit_table_rows,
        '{2D_MISFIT_TABLE_ROWS_FC_MW}': two_d_misfit_table_rows_fc_mw,
        '{2D_MISFIT_TABLE_ROWS_FC_TSTAR}': two_d_misfit_table_rows_fc_tstar,
        '{2D_MISFIT_TABLE_ROWS_TSTAR_MW}': two_d_misfit_table_rows_tstar_mw
    }
    with open(misfit_html, encoding='utf-8') as fp:
        misfit = fp.read()
    misfit = _multireplace(misfit, replacements)
    with open(misfit_html_out, 'w', encoding='utf-8') as fp:
        fp.write(misfit)


def _add_run_info_to_html(replacements):
    """
    Add run info to HTML report.

    :param replacements: replacement dictionary
    :type replacements: dict
    """
    ssp_url = 'https://sourcespec.seismicsource.org'
    logo_file = _logo_file_url()
    agency_logo = _agency_logo()
    ssp_version, run_completed = _version_and_run_completed()
    author = _author_html()
    if not author:
        author_comment_begin = '<!--'
        author_comment_end = '-->'
    else:
        author_comment_begin = ''
        author_comment_end = ''
    agency = _agency_html()
    if not agency:
        agency_comment_begin = '<!--'
        agency_comment_end = '-->'
    else:
        agency_comment_begin = ''
        agency_comment_end = ''
    author_and_agency = _author_and_agency_html(author, agency)
    page_footer = _page_footer()
    replacements.update({
        '{AGENCY_LOGO}': agency_logo,
        '{LOGO_FILE}': logo_file,
        '{VERSION}': ssp_version,
        '{SSP_URL}': ssp_url,
        '{RUN_COMPLETED}': run_completed,
        '{AUTHOR}': author,
        '{AUTHOR_COMMENT_BEGIN}': author_comment_begin,
        '{AUTHOR_COMMENT_END}': author_comment_end,
        '{AGENCY}': agency,
        '{AGENCY_COMMENT_BEGIN}': agency_comment_begin,
        '{AGENCY_COMMENT_END}': agency_comment_end,
        '{AUTHOR_AND_AGENCY}': author_and_agency,
        '{PAGE_FOOTER}': page_footer,
    })


def _add_event_info_to_html(replacements):
    """
    Add event info to HTML report.

    :param replacements: replacement dictionary
    :type replacements: dict
    """
    evid = config.event.event_id
    evname = config.event.name
    hypo = config.event.hypocenter
    run_id = config.options.run_id
    replacements.update({
        '{EVENTID}': evid,
        '{EVENT_NAME}': evname,
        '{RUNID}': run_id,
        '{EVENT_LONGITUDE}': f'{hypo.longitude.value_in_deg:8.3f}',
        '{EVENT_LATITUDE}': f'{hypo.latitude.value_in_deg:7.3f}',
        '{EVENT_DEPTH}': f'{hypo.depth.value_in_km:5.1f}',
        '{ORIGIN_TIME}': f'{hypo.origin_time}',
    })
    # Link to event page, if defined
    event_url = config.event_url
    if event_url is not None:
        event_url = event_url.replace('$EVENTID', evid)
        parsed_url = urlparse(event_url)
        if not parsed_url.scheme:
            logger.warning(
                f'{event_url} is not a valid URL and will not be used')
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
    # Only show Event Name if it is not empty
    if evname:
        evname_comment_begin = evname_comment_end = ''
    else:
        evname_comment_begin = '<!--'
        evname_comment_end = '-->'
    replacements.update({
        '{EVENT_NAME_COMMENT_BEGIN}': evname_comment_begin,
        '{EVENT_NAME_COMMENT_END}': evname_comment_end
    })
    # Only show Run ID if it is not empty
    if run_id:
        run_id_comment_begin = run_id_comment_end = ''
    else:
        run_id_comment_begin = '<!--'
        run_id_comment_end = '-->'
    replacements.update({
        '{RUNID_COMMENT_BEGIN}': run_id_comment_begin,
        '{RUNID_COMMENT_END}': run_id_comment_end
    })


def _add_maps_to_html(replacements):
    """
    Add maps to HTML report.

    :param replacements: replacement dictionary
    :type replacements: dict
    """
    try:
        station_maps = [
            m for m in config.figures['station_maps']
            if m.endswith(VALID_FIGURE_FORMATS)]
    except KeyError:
        station_maps = []
    map_mag = ''
    map_fc = ''
    with contextlib.suppress(IndexError):
        map_mag = [
            mapfile for mapfile in station_maps if 'map_mag' in mapfile][0]
        map_mag = os.path.basename(map_mag)
    with contextlib.suppress(IndexError):
        map_fc = [
            mapfile for mapfile in station_maps if 'map_fc' in mapfile][0]
        map_fc = os.path.basename(map_fc)
    if map_mag:
        map_mag_comment_begin = map_mag_comment_end = ''
    else:
        map_mag_comment_begin = '<!--'
        map_mag_comment_end = '-->'
    if map_fc:
        map_fc_comment_begin = map_fc_comment_end = ''
    else:
        map_fc_comment_begin = '<!--'
        map_fc_comment_end = '-->'
    replacements.update({
        '{MAP_MAG}': map_mag,
        '{MAP_MAG_COMMENT_BEGIN}': map_mag_comment_begin,
        '{MAP_MAG_COMMENT_END}': map_mag_comment_end,
        '{MAP_FC}': map_fc,
        '{MAP_FC_COMMENT_BEGIN}': map_fc_comment_begin,
        '{MAP_FC_COMMENT_END}': map_fc_comment_end
    })


def _add_traces_plots_to_html(templates, replacements):
    """
    Add trace plots to HTML report.

    :param templates: template files
    :type templates: :class:`HTMLTemplates`
    :param replacements: replacement dictionary
    :type replacements: dict
    """
    with open(templates.traces_plot_html, encoding='utf-8') as fp:
        traces_plot = fp.read()
    traces_plot_files = [
        t for t in config.figures['traces']
        if t.endswith(VALID_FIGURE_FORMATS)]
    traces_plots = ''
    traces_plot_class = ''
    n_traces_plot_files = len(traces_plot_files)
    for n, traces_plot_file in enumerate(sorted(traces_plot_files)):
        if n_traces_plot_files > 1:
            traces_plot_counter = (
                f'<span class="print_inline">&nbsp;({n+1} of '
                f'{n_traces_plot_files})</span>'
            )
        else:
            traces_plot_counter = ''
        traces_plot_file = os.path.basename(traces_plot_file)
        traces_plots += traces_plot.\
            replace('{TRACES_PLOT_CLASS}', traces_plot_class).\
            replace('{TRACES_PLOT_COUNTER}', traces_plot_counter).\
            replace('{TRACES_PLOT_FILE}', traces_plot_file)
        traces_plot_class = ' class="print"'
    if traces_plots:
        traces_plots_comment_begin = traces_plots_comment_end = ''
    else:
        traces_plots_comment_begin = '<!--'
        traces_plots_comment_end = '-->'
    replacements.update({
        '{TRACES_PLOTS}': traces_plots,
        '{TRACES_PLOTS_COMMENT_BEGIN}': traces_plots_comment_begin,
        '{TRACES_PLOTS_COMMENT_END}': traces_plots_comment_end
    })


def _add_spectra_plots_to_html(templates, replacements):
    """
    Add spectra plots to HTML report.

    :param templates: template files
    :type templates: :class:`HTMLTemplates`
    :param replacements: replacement dictionary
    :type replacements: dict
    """
    with open(templates.spectra_plot_html, encoding='utf-8') as fp:
        spectra_plot = fp.read()
    spectra_plot_files = [
        s for s in config.figures['spectra_regular']
        if s.endswith(VALID_FIGURE_FORMATS)]
    spectra_plots = ''
    spectra_plot_class = ''
    n_spectra_plot_files = len(spectra_plot_files)
    for n, spectra_plot_file in enumerate(sorted(spectra_plot_files)):
        if n_spectra_plot_files > 1:
            spectra_plot_counter = (
                f'<span class="print_inline">&nbsp;({n+1} of '
                f'{n_spectra_plot_files})</span>'
            )
        else:
            spectra_plot_counter = ''
        spectra_plot_file = os.path.basename(spectra_plot_file)
        spectra_plots += spectra_plot.\
            replace('{SPECTRA_PLOT_CLASS}', spectra_plot_class).\
            replace('{SPECTRA_PLOT_COUNTER}', spectra_plot_counter).\
            replace('{SPECTRA_PLOT_FILE}', spectra_plot_file)
        spectra_plot_class = ' class="print"'
    if spectra_plots:
        spectra_plots_comment_begin = spectra_plots_comment_end = ''
    else:
        spectra_plots_comment_begin = '<!--'
        spectra_plots_comment_end = '-->'
    replacements.update({
        '{SPECTRA_PLOTS}': spectra_plots,
        '{SPECTRA_PLOTS_COMMENT_BEGIN}': spectra_plots_comment_begin,
        '{SPECTRA_PLOTS_COMMENT_END}': spectra_plots_comment_end
    })


def _add_inversion_info_to_html(sspec_output, replacements):
    """
    Add inversion info to HTML report.

    :param sspec_output: SourceSpec output
    :type sspec_output: :class:`~sourcespec.ssp_data_types.SourceSpecOutput`
    :param replacements: replacement dictionary
    :type replacements: dict
    """
    # Inversion information
    inversion_algorithms = {
        'TNC': 'Truncated Newton',
        'LM': 'Levenberg-Marquardt',
        'BH': 'Basin-hopping',
        'GS': 'Grid search',
        'IS': 'K-d tree importance sampling',
    }
    weightings = {
        'noise': 'Noise weighting',
        'frequency': 'Frequency weighting',
        'inv_frequency': 'Inverse frequency weighting',
        'no_weight': 'No weighting',
    }
    inversion_algorithm = inversion_algorithms[
        sspec_output.inversion_info.algorithm]
    inversion_wave_type = sspec_output.inversion_info.wave_type
    inversion_weighting = weightings[sspec_output.inversion_info.weighting]
    inversion_t_star_0 = f'{sspec_output.inversion_info.t_star_0} s'
    inversion_invert_t_star_0 =\
        str(sspec_output.inversion_info.invert_t_star_0)
    inversion_t_star_0_variability =\
        f'{sspec_output.inversion_info.t_star_0_variability * 100:.1f} %'
    if sspec_output.inversion_info.t_star_min_max == 'null':
        inversion_t_star_min_max = '-'
    else:
        inversion_t_star_min_max =\
            f'{sspec_output.inversion_info.t_star_min_max} s'
    if sspec_output.inversion_info.fc_min_max == 'null':
        inversion_fc_min_max = '-'
    else:
        inversion_fc_min_max =\
            f'{sspec_output.inversion_info.fc_min_max} Hz'
    if sspec_output.inversion_info.Qo_min_max == 'null':
        inversion_Qo_min_max = '-'
    else:
        inversion_Qo_min_max =\
            str(sspec_output.inversion_info.Qo_min_max)
    replacements.update({
        '{INVERSION_ALGORITHM}': inversion_algorithm,
        '{INVERSION_WAVE_TYPE}': inversion_wave_type,
        '{INVERSION_WEIGHTING}': inversion_weighting,
        '{INVERSION_T_STAR_0}': inversion_t_star_0,
        '{INVERSION_INVERT_T_STAR_0}': inversion_invert_t_star_0,
        '{INVERSION_T_STAR_0_VARIABILITY}': inversion_t_star_0_variability,
        '{INVERSION_T_STAR_MIN_MAX}': inversion_t_star_min_max,
        '{INVERSION_FC_MIN_MAX}': inversion_fc_min_max,
        '{INVERSION_Q0_MIN_MAX}': inversion_Qo_min_max,
    })


def _add_summary_spectral_params_to_html(sspec_output, replacements):
    """
    Add summary spectral parameters to HTML report.

    :param sspec_output: SourceSpec output
    :type sspec_output: :class:`~sourcespec.ssp_data_types.SourceSpecOutput`
    :param replacements: replacement dictionary
    :type replacements: dict
    """
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
        Mw_summary_str = (
            f'{Mw_summary:.2f} '
            f'[- {Mw_summary_error_minus:.2f}, + {Mw_summary_error_plus:.2f}]'
        )
    except TypeError:
        Mw_summary_error = summary_uncertainties['Mw']
        Mw_summary_str = f'{Mw_summary:.2f} ± {Mw_summary_error:.2f}'
    replacements.update({'{MW_SUMMARY}': Mw_summary_str})
    fc_summary = summary_values['fc']
    try:
        fc_summary_error_minus, fc_summary_error_plus =\
            summary_uncertainties['fc']
        fc_summary_str = (
            f'{fc_summary:.3f} '
            f'[- {fc_summary_error_minus:.3f}, + {fc_summary_error_plus:.3f}]'
        )
    except TypeError:
        fc_summary_error = summary_uncertainties['fc']
        fc_summary_str = f'{fc_summary:.3f} ± {fc_summary_error:.3f}'
    replacements.update({'{FC_SUMMARY}': fc_summary_str})

    means = sspec_output.mean_values()
    mean_errors = sspec_output.mean_uncertainties()
    wmeans = sspec_output.weighted_mean_values()
    wmean_errors = sspec_output.weighted_mean_uncertainties()
    percentiles = sspec_output.percentiles_values()
    percentile_errors = sspec_output.percentiles_uncertainties()

    n_sigma = config.n_sigma
    n_sigma = int(n_sigma) if float(n_sigma).is_integer() else n_sigma
    n_sigma = f'{n_sigma} sigma'
    mid_pct, lower_pct, upper_pct =\
        config.mid_percentage, config.lower_percentage, config.upper_percentage
    mid_pct = int(mid_pct) if float(mid_pct).is_integer() else mid_pct
    lower_pct = int(lower_pct) if float(lower_pct).is_integer() else lower_pct
    upper_pct = int(upper_pct) if float(upper_pct).is_integer() else upper_pct
    percentages = f'{mid_pct}%, [{lower_pct}%, {upper_pct}%]'
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
    ssd_mean = means['ssd']
    ssd_mean_error = mean_errors['ssd']
    ssd_wmean = wmeans['ssd']
    ssd_wmean_error = wmean_errors['ssd']
    ssd_perc = percentiles['ssd']
    ssd_perc_error = percentile_errors['ssd']
    replacements.update({
        '{SSD_MEAN_AND_ERR}': _summary_value_and_err_text(
            ssd_mean, ssd_mean_error, '{:.3e}'),
        '{SSD_WMEAN_AND_ERR}': _summary_value_and_err_text(
            ssd_wmean, ssd_wmean_error, '{:.3e}'),
        '{SSD_PERC_AND_ERR}': _summary_value_and_err_text(
            ssd_perc, ssd_perc_error, '{:.3e}'),
    })
    Er_mean = means['Er']
    Er_mean_error = mean_errors['Er']
    Er_wmean = wmeans['Er']
    Er_wmean_error = wmean_errors['Er']
    Er_perc = percentiles['Er']
    Er_perc_error = percentile_errors['Er']
    replacements.update({
        '{ER_MEAN_AND_ERR}': _summary_value_and_err_text(
            Er_mean, Er_mean_error, '{:.3e}'),
        '{ER_WMEAN_AND_ERR}': _summary_value_and_err_text(
            Er_wmean, Er_wmean_error, '{:.3e}'),
        '{ER_PERC_AND_ERR}': _summary_value_and_err_text(
            Er_perc, Er_perc_error, '{:.3e}'),
    })
    sigma_a_mean = means['sigma_a']
    sigma_a_mean_error = mean_errors['sigma_a']
    sigma_a_wmean = wmeans['sigma_a']
    sigma_a_wmean_error = wmean_errors['sigma_a']
    sigma_a_perc = percentiles['sigma_a']
    sigma_a_perc_error = percentile_errors['sigma_a']
    replacements.update({
        '{SIGMA_A_MEAN_AND_ERR}': _summary_value_and_err_text(
            sigma_a_mean, sigma_a_mean_error, '{:.3e}'),
        '{SIGMA_A_WMEAN_AND_ERR}': _summary_value_and_err_text(
            sigma_a_wmean, sigma_a_wmean_error, '{:.3e}'),
        '{SIGMA_A_PERC_AND_ERR}': _summary_value_and_err_text(
            sigma_a_perc, sigma_a_perc_error, '{:.3e}'),
    })
    # Local magnitude, if computed
    if config.compute_local_magnitude:
        Ml_mean = means['Ml']
        Ml_mean_error = mean_errors['Ml']
        Ml_wmean = wmeans['Ml']
        Ml_wmean_error = wmean_errors['Ml']
        Ml_perc = percentiles['Ml']
        Ml_perc_error = percentile_errors['Ml']
        Ml_comment_begin = ''
        Ml_comment_end = ''
    else:
        Ml_mean = Ml_wmean = Ml_perc = np.nan
        Ml_mean_error = Ml_wmean_error = Ml_perc_error = (np.nan, np.nan)
        Ml_comment_begin = '<!--'
        Ml_comment_end = '-->'
    replacements.update({
        '{ML_MEAN_AND_ERR}': _summary_value_and_err_text(
            Ml_mean, Ml_mean_error, '{:.2f}'),
        '{ML_WMEAN_AND_ERR}': _summary_value_and_err_text(
            Ml_wmean, Ml_wmean_error, '{:.2f}'),
        '{ML_PERC_AND_ERR}': _summary_value_and_err_text(
            Ml_perc, Ml_perc_error, '{:.2f}'),
        '{ML_COMMENT_BEGIN}': Ml_comment_begin,
        '{ML_COMMENT_END}': Ml_comment_end,
    })


def _add_box_plots_to_html(replacements):
    """
    Add box plots to HTML report.

    :param replacements: replacement dictionary
    :type replacements: dict
    """
    box_plots = ''
    with contextlib.suppress(KeyError, IndexError):
        box_plots = [
            b for b in config.figures['boxplots']
            if b.endswith(VALID_FIGURE_FORMATS)][0]
        box_plots = os.path.basename(box_plots)
    if box_plots:
        box_plots_comment_begin = box_plots_comment_end = ''
    else:
        box_plots_comment_begin = '<!--'
        box_plots_comment_end = '-->'
    replacements.update({
        '{BOX_PLOTS}': box_plots,
        '{BOX_PLOTS_COMMENT_BEGIN}': box_plots_comment_begin,
        '{BOX_PLOTS_COMMENT_END}': box_plots_comment_end,
    })


def _add_stacked_spectra_to_html(replacements):
    """
    Add stacked spectra to HTML report.

    :param replacements: replacement dictionary
    :type replacements: dict
    """
    stacked_spectra = ''
    with contextlib.suppress(KeyError, IndexError):
        stacked_spectra = [
            s for s in config.figures['stacked_spectra']
            if s.endswith(VALID_FIGURE_FORMATS)][0]
        stacked_spectra = os.path.basename(stacked_spectra)
    if stacked_spectra:
        stacked_spectra_comment_begin = stacked_spectra_comment_end = ''
    else:
        stacked_spectra_comment_begin = '<!--'
        stacked_spectra_comment_end = '-->'
    replacements.update({
        '{STACKED_SPECTRA}': stacked_spectra,
        '{STACKED_SPECTRA_COMMENT_BEGIN}': stacked_spectra_comment_begin,
        '{STACKED_SPECTRA_COMMENT_END}': stacked_spectra_comment_end,
    })


def _add_station_table_to_html(sspec_output, templates, replacements):
    """
    Add station table to HTML report.

    :param sspec_output: SourceSpec output
    :type sspec_output: :class:`~sourcespec.ssp_data_types.SourceSpecOutput`
    :param templates: template files
    :type templates: :class:`HTMLTemplates`
    :param replacements: replacement dictionary
    :type replacements: dict
    """
    with open(templates.station_table_row_html, encoding='utf-8') as fp:
        station_table_row = fp.read()
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
        ssd_text, ssd_err_text =\
            _station_value_and_err_text(par, 'ssd', '{:.3e}')
        ra_text, ra_err_text =\
            _station_value_and_err_text(par, 'radius', '{:.3f}')
        Er_text, _ = _station_value_and_err_text(par, 'Er', '{:.3e}')
        sigma_a_text, sigma_a_err_text =\
            _station_value_and_err_text(par, 'sigma_a', '{:.3e}')
        Ml_text, _ = _station_value_and_err_text(par, 'Ml', '{:.3f}')
        hyp_dist_text, _ =\
            _station_value_and_err_text(par, 'hypo_dist_in_km', '{:.3f}')
        az_text, _ = _station_value_and_err_text(par, 'azimuth', '{:.3f}')
        row_replacements = {
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
            '{STATION_SSD}': ssd_text,
            '{STATION_SSD_ERR}': ssd_err_text,
            '{STATION_RA}': ra_text,
            '{STATION_RA_ERR}': ra_err_text,
            '{STATION_SIGMA_A}': sigma_a_text,
            '{STATION_SIGMA_A_ERR}': sigma_a_err_text,
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
        row_replacements.update({
            '{ML_COMMENT_BEGIN}': Ml_comment_begin,
            '{ML_COMMENT_END}': Ml_comment_end
        })
        station_table_rows += _multireplace(
            station_table_row, row_replacements)
    replacements.update({
        '{STATION_TABLE_ROWS}': station_table_rows,
    })


def _add_misfit_plots_to_html(replacements):
    """
    Add misfit plots to HTML report.

    :param replacements: replacement dictionary
    :type replacements: dict
    """
    if 'misfit_1d' in config.figures:
        misfit_plot_comment_begin = ''
        misfit_plot_comment_end = ''
        _misfit_page()
    else:
        misfit_plot_comment_begin = '<!--'
        misfit_plot_comment_end = '-->'
    replacements.update({
        '{MISFIT_PLOT_COMMENT_BEGIN}': misfit_plot_comment_begin,
        '{MISFIT_PLOT_COMMENT_END}': misfit_plot_comment_end
    })


def _add_downloadable_files_to_html(templates, replacements):
    """
    Add links to downloadable files to HTML report.

    :param templates: template files
    :type templates: :class:`HTMLTemplates`
    :param replacements: replacement dictionary
    :type replacements: dict
    """
    # symlink to input files (not supported on Windows)
    input_files = '' if os.name == 'nt' else 'input_files'
    input_files_text = '' if os.name == 'nt'\
        else 'Click to navigate to input files'
    evid = config.event.event_id
    config_file = f'{evid}.ssp.conf'
    yaml_file = f'{evid}.ssp.yaml'
    log_file = f'{evid}.ssp.log'

    replacements.update({
        '{INPUT_FILES}': input_files,
        '{INPUT_FILES_TEXT}': input_files_text,
        '{CONF_FILE}': config_file,
        '{YAML_FILE}': yaml_file,
        '{LOG_FILE}': log_file,
    })

    # QuakeML file (if produced)
    if config.qml_file_out is not None:
        quakeml_file = os.path.basename(config.qml_file_out)
        with open(templates.quakeml_file_link_html, encoding='utf-8') as fp:
            quakeml_file_link = fp.read()
        quakeml_file_link = quakeml_file_link\
            .replace('{QUAKEML_FILE}', quakeml_file)
    else:
        quakeml_file_link = ''
    replacements.update({
        '{QUAKEML_FILE_LINK}': quakeml_file_link
    })

    suppl_file_list = [
        os.path.basename(fig) for fig in config.figures['traces_raw']
    ]
    suppl_file_list += [
        os.path.basename(fig) for fig in config.figures['spectra_weight']
    ]
    supplementary_files = ''.join(
        f'<a href="{suppl_file}">{suppl_file}</a></br>\n'
        for suppl_file in suppl_file_list
    )
    if supplementary_files:
        with open(
            templates.supplementary_file_links_html, encoding='utf-8'
        ) as fp:
            supplementary_file_links = fp.read()
        supplementary_file_links = supplementary_file_links\
            .replace('{SUPPLEMENTARY_FILES}', supplementary_files)
    else:
        supplementary_file_links = ''
    replacements.update({
        '{SUPPLEMENTARY_FILE_LINKS}': supplementary_file_links
    })


class HTMLtemplates:
    """Class to hold paths to HTML templates."""
    def __init__(self):
        template_dir = os.path.join(
            os.path.dirname(os.path.realpath(__file__)),
            'html_report_template'
        )
        self.style_css = os.path.join(template_dir, 'style.css')
        self.index_html = os.path.join(template_dir, 'index.html')
        self.traces_plot_html = os.path.join(
            template_dir, 'traces_plot.html')
        self.spectra_plot_html = os.path.join(
            template_dir, 'spectra_plot.html')
        self.station_table_row_html = os.path.join(
            template_dir, 'station_table_row.html')
        self.quakeml_file_link_html = os.path.join(
            template_dir, 'quakeml_file_link.html')
        self.supplementary_file_links_html = os.path.join(
            template_dir, 'supplementary_file_links.html')


def _cleanup_html(text):
    """
    Remove unnecessary comments and whitespace from HTML.

    :param text: HTML text
    :type text: str

    :returns: cleaned HTML text
    :rtype: str
    """
    # remove HTML-style comments
    text = re.sub(r'<!--.*?-->', '', text, flags=re.DOTALL)
    # strip spaces at the end of lines
    text = re.sub(r' +$', '', text, flags=re.MULTILINE)
    # replace multiple empty lines with a single empty line
    text = re.sub(r'\n\s*\n', '\n\n', text)
    return text


def html_report(sspec_output):
    """
    Generate an HTML report.

    :param sspec_output: Output from the SourceSpec inversion.
    :type sspec_output: :class:`sourcespec.ssp_data_types.SourceSpecOutput`
    """
    templates = HTMLtemplates()

    replacements = {}
    _add_run_info_to_html(replacements)
    _add_event_info_to_html(replacements)
    _add_maps_to_html(replacements)
    _add_traces_plots_to_html(templates, replacements)
    _add_spectra_plots_to_html(templates, replacements)
    _add_inversion_info_to_html(sspec_output, replacements)
    _add_summary_spectral_params_to_html(sspec_output, replacements)
    _add_box_plots_to_html(replacements)
    _add_stacked_spectra_to_html(replacements)
    _add_station_table_to_html(sspec_output, templates, replacements)
    _add_misfit_plots_to_html(replacements)
    _add_downloadable_files_to_html(templates, replacements)

    with open(templates.index_html, encoding='utf-8') as fp:
        index = fp.read()
    index = _multireplace(index, replacements)
    index = _cleanup_html(index)
    shutil.copy(templates.style_css, config.options.outdir)
    index_html_out = os.path.join(config.options.outdir, 'index.html')
    with open(index_html_out, 'w', encoding='utf-8') as fp:
        fp.write(index)
    logger.info(f'HTML report written to file: {index_html_out}')
