# -*- coding: utf-8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Output functions for source_spec.

:copyright:
    2012 Claudio Satriano <satriano@ipgp.fr>

    2013-2014 Claudio Satriano <satriano@ipgp.fr>,
              Emanuela Matrullo <matrullo@geologie.ens.fr>

    2015-2022 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
import os
import logging
import numpy as np
from collections.abc import Mapping
from datetime import datetime
from tzlocal import get_localzone
from sourcespec.ssp_data_types import sspec_out_comments
from sourcespec.ssp_qml_output import write_qml
from sourcespec.ssp_sqlite_output import write_sqlite
from sourcespec._version import get_versions
logger = logging.getLogger(__name__.split('.')[-1])


def _format_exponent(value, reference):
    """Format `value` to a string having the same exponent than `reference`."""
    # get the exponent of reference value
    xp = np.int(np.floor(np.log10(np.abs(reference))))
    # format value to print it with the same exponent of reference value
    return '{:5.3f}e{:+03d}'.format(value/10**xp, xp)


def _write_author_and_agency_to_parfile(config, parfile):
    author_str = empty_author_str = '\n*** Author:'
    if config.author_name is not None:
        author_str += ' {}'.format(config.author_name)
    if config.author_email is not None:
        if author_str != empty_author_str:
            author_str += ' <{}>'.format(config.author_email)
        else:
            author_str += ' {}'.format(config.author_email)
    if author_str != empty_author_str:
        parfile.write(author_str)
    agency_str = empty_agency_str = '\n*** Agency:'
    if config.agency_full_name is not None:
        agency_str += ' {}'.format(config.agency_full_name)
    if config.agency_short_name is not None:
        if agency_str != empty_agency_str:
            agency_str += ' ({})'.format(config.agency_short_name)
        else:
            agency_str += ' {}'.format(config.agency_short_name)
    if config.agency_url is not None:
        if agency_str != empty_agency_str:
            agency_str += ' -'
        agency_str += ' {}'.format(config.agency_url)
    if agency_str != empty_agency_str:
        parfile.write(agency_str)


def _write_parfile(config, sspec_output):
    """
    Write station source parameters to file.

    Note: this format is deprecated and will not evolve anymore
    (e.g., including new parameters or new way of computing statistics).
    """
    if not os.path.exists(config.options.outdir):
        os.makedirs(config.options.outdir)
    evid = config.hypo.evid
    parfilename = os.path.join(
        config.options.outdir, '{}.ssp.out'.format(evid))
    parfile = open(parfilename, 'w')
    parfile.write(
        '*** Note: this file is deprecated and might not contain all the '
        'output information. ***\n')
    parfile.write(
        '*** It will be removed in future versions of SourceSpec. ***\n')
    parfile.write(
        '*** Please look at the YAML file in this directory. ***\n')

    hypo = config.hypo
    parfile.write(
        '{} lon {:8.3f} lat {:7.3f} depth {:5.1f} km '
        'orig_time {}\n\n'.format(
            hypo.evid, hypo.longitude, hypo.latitude, hypo.depth,
            hypo.origin_time))
    parfile.write('*** Station source parameters ***\n')
    parfile.write(
        '*** Note: outliers are prepended by a star (*) symbol ***\n')
    parkeys = (
        'Mw', 'fc', 't_star', 'Qo', 'Mo',
        'bsd', 'ra', 'hyp_dist', 'az', 'Er'
    )
    formats = dict(
        Mo='{:.3e} ',
        Er='{:.3e} ',
        hyp_dist='{:7.3f} ',
        az='{:7.3f} ',
        Mw='{:6.3f} ',
        fc='{:6.3f} ',
        bsd='{:.3e} ',
        ra='{:8.3f} ',
        t_star='{:6.3f} ',
        Qo='{:7.1f} ',
        Ml='{:6.3f} '
    )
    formats_none = dict(
        Mo='{:>9} ',
        Er='{:>9} ',
        hyp_dist='{:>7} ',
        az='{:>7} ',
        Mw='{:>6} ',
        fc='{:>6} ',
        bsd='{:>9} ',
        ra='{:>8} ',
        t_star='{:>6} ',
        Qo='{:>7} ',
        Ml='{:>6} '
    )
    stationpar = sspec_output.station_parameters
    for statId in sorted(stationpar.keys()):
        par = stationpar[statId]
        parfile.write('{:>15} {:>6}\t'.format(statId, par.instrument_type))
        for key in parkeys:
            if key == 'ra':
                _pkey = 'radius'
            else:
                _pkey = key
            if key == 'hyp_dist':
                val = par['hypo_dist_in_km']
                outl = False
            elif key == 'az':
                val = par['azimuth']
                outl = False
            else:
                val = par[_pkey].value
                outl = par[_pkey].outlier
            if outl:
                space = ' *'
            else:
                space = '  '
            parfile.write('{}{} '.format(space, key))
            if val is not None and ~np.isnan(val):
                parfile.write(formats[key].format(val))
            else:
                parfile.write(formats_none[key].format('nan'))
        parfile.write('\n')
        parfile.write('{:>22}\t'.format('--- errmin'))
        for key in parkeys:
            if key == 'ra':
                _pkey = 'radius'
            else:
                _pkey = key
            if key in ['hyp_dist', 'az']:
                outl = False
                err = None
            else:
                outl = par[_pkey].outlier
                err = par[_pkey].lower_uncertainty
                if err is None:
                    err = par[_pkey].uncertainty
            if outl:
                space = ' *'
            else:
                space = '  '
            parfile.write('{}{} '.format(space, key))
            if err is not None:
                parfile.write(formats[key].format(err))
            else:
                parfile.write(formats_none[key].format('nan'))
        parfile.write('\n')
        parfile.write('{:>22}\t'.format('--- errmax'))
        for key in parkeys:
            if key == 'ra':
                _pkey = 'radius'
            else:
                _pkey = key
            if key in ['hyp_dist', 'az']:
                outl = False
                err = None
            else:
                outl = par[_pkey].outlier
                err = par[_pkey].upper_uncertainty
                if err is None:
                    err = par[_pkey].uncertainty
            if outl:
                space = ' *'
            else:
                space = '  '
            parfile.write('{}{} '.format(space, key))
            if err is not None:
                parfile.write(formats[key].format(err))
            else:
                parfile.write(formats_none[key].format('nan'))
        parfile.write('\n')

    means = sspec_output.mean_values()
    errors = sspec_output.mean_uncertainties()
    means_weight = sspec_output.weighted_mean_values()
    errors_weight = sspec_output.weighted_mean_uncertainties()

    parfile.write('\n*** Average source parameters ***\n')
    parfile.write('*** Note: averages computed after removing outliers ****\n')

    Mw_mean = means['Mw']
    Mw_error = errors['Mw']
    parfile.write('Mw: {:.2f} +/- {:.2f}\n'.format(Mw_mean, Mw_error))
    Mw_mean_weight = means_weight['Mw']
    Mw_error_weight = errors_weight['Mw']
    parfile.write('Mw (weighted): {:.2f} +/- {:.2f}\n'.format(
        Mw_mean_weight, Mw_error_weight))

    Mo_mean = means['Mo']
    Mo_minus, Mo_plus = errors['Mo']
    # format Mo_plus and Mo_minus to print it with the same exponent of Mo
    Mo_minus_str = _format_exponent(Mo_minus, Mo_mean)
    Mo_plus_str = _format_exponent(Mo_plus, Mo_mean)
    parfile.write('Mo: {:.3e} /- {} /+ {} N.m\n'.format(
        Mo_mean, Mo_minus_str, Mo_plus_str))
    Mo_mean_weight = means_weight['Mo']
    Mo_minus_weight, Mo_plus_weight = errors_weight['Mo']
    # format Mo_plus and Mo_minus to print it with the same exponent of Mo
    Mo_minus_str = _format_exponent(Mo_minus_weight, Mo_mean_weight)
    Mo_plus_str = _format_exponent(Mo_plus_weight, Mo_mean_weight)
    parfile.write('Mo (weighted): {:.3e} /- {} /+ {} N.m\n'.format(
        Mo_mean_weight, Mo_minus_str, Mo_plus_str))

    fc_mean = means['fc']
    fc_minus, fc_plus = errors['fc']
    parfile.write('fc: {:.3f} /- {:.3f} /+ {:.3f} Hz\n'.format(
        fc_mean, fc_minus, fc_plus))
    fc_mean_weight = means_weight['fc']
    fc_minus_weight, fc_plus_weight = errors_weight['fc']
    parfile.write('fc (weighted): {:.3f} /- {:.3f} /+ {:.3f} Hz\n'.format(
        fc_mean_weight, fc_minus_weight, fc_plus_weight))

    t_star_mean = means['t_star']
    t_star_error = errors['t_star']
    parfile.write('t_star: {:.3f} +/- {:.3f} s\n'.format(
        t_star_mean, t_star_error))
    t_star_mean_weight = means_weight['t_star']
    t_star_error_weight = errors_weight['t_star']
    parfile.write('t_star (weighted): {:.3f} +/- {:.3f} s\n'.format(
        t_star_mean_weight, t_star_error_weight))

    Qo_mean = means['Qo']
    Qo_error = errors['Qo']
    parfile.write('Qo: {:.1f} +/- {:.1f}\n'.format(Qo_mean, Qo_error))
    Qo_mean_weight = means_weight['Qo']
    Qo_error_weight = errors_weight['Qo']
    parfile.write('Qo (weighted): {:.1f} +/- {:.1f}\n'.format(
        Qo_mean_weight, Qo_error_weight))

    ra_mean = means['radius']
    ra_minus, ra_plus = errors['radius']
    parfile.write('Source radius: {:.3f} /- {:.3f} /+ {:.3f} m\n'.format(
        ra_mean, ra_minus, ra_plus))
    ra_mean_weight = means_weight['radius']
    ra_minus_weight, ra_plus_weight = errors_weight['radius']
    parfile.write(
        'Source radius (weighted): {:.3f} /- {:.3f} /+ {:.3f} m\n'.format(
            ra_mean_weight, ra_minus_weight, ra_plus_weight))

    bsd_mean = means['bsd']
    bsd_minus, bsd_plus = errors['bsd']
    bsd_minus_str = _format_exponent(bsd_minus, bsd_mean)
    bsd_plus_str = _format_exponent(bsd_plus, bsd_mean)
    parfile.write('Brune stress drop: {:.3e} /- {} /+ {} MPa\n'.format(
        bsd_mean, bsd_minus_str, bsd_plus_str))
    bsd_mean_weight = means_weight['bsd']
    bsd_minus_weight, bsd_plus_weight = errors_weight['bsd']
    bsd_minus_str = _format_exponent(bsd_minus_weight, bsd_mean_weight)
    bsd_plus_str = _format_exponent(bsd_plus_weight, bsd_mean_weight)
    parfile.write(
        'Brune stress drop (weighted): {:.3e} /- {} /+ {} MPa\n'.format(
            bsd_mean_weight, bsd_minus_str, bsd_plus_str))

    Ml_mean = means.get('Ml', None)
    Ml_error = errors.get('Ml', None)
    if Ml_mean is not None:
        parfile.write('Ml: {:.3f} +/- {:.3f}\n'.format(Ml_mean, Ml_error))

    Er_mean = means['Er']
    Er_minus, Er_plus = errors['Er']
    # format Er_plus and Er_minus to print it with the same exponent of Er
    Er_minus_str = _format_exponent(Er_minus, Er_mean)
    Er_plus_str = _format_exponent(Er_plus, Er_mean)
    parfile.write('Er: {:.3e} /- {} /+ {} N.m\n'.format(
        Er_mean, Er_minus_str, Er_plus_str))

    parfile.write('\n*** SourceSpec: {}'.format(get_versions()['version']))
    parfile.write('\n*** Run completed on: {} {}'.format(
        config.end_of_run, config.end_of_run_tz))
    if config.options.run_id:
        parfile.write('\n*** Run ID: {}'.format(config.options.run_id))
    _write_author_and_agency_to_parfile(config, parfile)

    parfile.close()

    logger.info('Output written to file: ' + parfilename)


def _dict2yaml(dict_like, level=0, comments={}):
    """Serialize a dict-like object into YAML format."""
    if not isinstance(dict_like, Mapping):
        return
    fmt = dict_like.get('_format', None)
    value_uncertainty_keys = (
        'value', 'uncertainty', 'lower_uncertainty', 'upper_uncertainty')
    target_dict = {
        key: (fmt.format(value)
              if fmt is not None and key in value_uncertainty_keys
              else value)
        for key, value in dict_like.items()
        if not key.startswith('_') and value is not None
    }
    indent = ' '*2*level
    # use oneliners for dict-like objects containing value and uncertainty keys
    if set(target_dict.keys()).intersection(set(value_uncertainty_keys)):
        oneliner = str(target_dict).replace("'", "")
        lines = indent + oneliner + '\n'
        return lines
    lines = ''
    for key, value in target_dict.items():
        if key.startswith('_'):
            continue
        if value is None:
            continue
        if isinstance(value, Mapping):
            if level == 0:
                lines += '\n'
            try:
                comment = comments[key]
                for line in comment.split('\n'):
                    lines += '# {}\n'.format(line)
            except KeyError:
                pass
            lines += '{}{}:\n'.format(indent, key)
            lines += _dict2yaml(value, level+1, comments)
        else:
            lines += '{}{}: {}\n'.format(indent, key, value)
    return lines


def _write_yaml(config, sspec_output):
    """Write sspec output in a YAML file."""
    if not os.path.exists(config.options.outdir):
        os.makedirs(config.options.outdir)
    evid = config.hypo.evid
    yamlfilename = os.path.join(
        config.options.outdir, '{}.ssp.yaml'.format(evid))
    lines = _dict2yaml(sspec_output, comments=sspec_out_comments)
    with open(yamlfilename, 'w') as fp:
        comment = sspec_out_comments['begin']
        for line in comment.split('\n'):
            fp.write('# {}\n'.format(line))
        fp.write(lines)
    logger.info('Output written to file: ' + yamlfilename)


def _write_hypo(config, sspec_output):
    if not config.options.hypo_file:
        return
    with open(config.options.hypo_file, 'r') as fp:
        line = fp.readline()
        # Check if first 10 digits of the line contain characters
        if any(c.isalpha() for c in line[0:10]):
            line1 = line
            line = fp.readline()
        line = list(line)

    summary_values = sspec_output.reference_values()
    mw_str = '{:03.2f}'.format(summary_values['Mw'])
    Ml = summary_values.get('Ml', None)
    if Ml is not None and ~np.isnan(Ml):
        ml_str = '{:03.2f}'.format(Ml)
    else:
        ml_str = ' '*4
    for i in range(0, 4):
        line[49+i] = mw_str[0+i]
        # line[45+i] = mw_str[0+i]
        line[69+i] = ml_str[0+i]
    outline = ''.join(line)
    evid = config.hypo.evid
    hypo_file_out = os.path.join(
        config.options.outdir, '{}.ssp.h'.format(evid))
    with open(hypo_file_out, 'w') as fp:
        try:
            fp.write(line1)
        except Exception:
            pass
        fp.write(outline)
    logger.info('Hypo file written to: ' + hypo_file_out)


def write_output(config, sspec_output):
    """Write results into different formats."""
    # Add run info to output object
    run_info = sspec_output.run_info
    run_info.SourceSpec_version = get_versions()['version']
    config.end_of_run = datetime.now()
    tz = get_localzone()
    config.end_of_run_tz = tz.tzname(config.end_of_run)
    run_info.run_completed = '{} {}'.format(
        config.end_of_run, config.end_of_run_tz)
    if config.options.run_id:
        run_info.run_id: config.options.run_id
    run_info.author_name = config.author_name
    run_info.author_email = config.author_email
    run_info.agency_full_name = config.agency_full_name
    run_info.agency_short_name = config.agency_short_name
    run_info.agency_url = config.agency_url
    # Write to parfile (deprecated)
    _write_parfile(config, sspec_output)
    # Write to YAML file
    _write_yaml(config, sspec_output)
    # Write to SQLite database, if requested
    write_sqlite(config, sspec_output)
    # Write to hypo file, if requested
    _write_hypo(config, sspec_output)
    # Write to quakeml file, if requested
    write_qml(config, sspec_output)
