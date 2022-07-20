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
import sqlite3
import numpy as np
from datetime import datetime
from tzlocal import get_localzone
from sourcespec.ssp_setup import ssp_exit
from sourcespec.ssp_qml_output import write_qml
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


def _write_parfile(config, sourcepar):
    """Write station source parameters to file."""
    if not os.path.exists(config.options.outdir):
        os.makedirs(config.options.outdir)
    evid = config.hypo.evid
    parfilename = os.path.join(
        config.options.outdir, '{}.ssp.out'.format(evid))
    parfile = open(parfilename, 'w')

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
    stationpar = sourcepar.station_parameters
    for statId in sorted(stationpar.keys()):
        par = stationpar[statId]
        parfile.write('{:>15} {:>6}\t'.format(*statId.split()))
        for key in parkeys:
            val = par[key]
            outl = par.get(key + '_outlier', False)
            if outl:
                space = ' *'
            else:
                space = '  '
            parfile.write('{}{} '.format(space, key))
            if val is not None:
                parfile.write(formats[key].format(val))
            else:
                parfile.write(formats_none[key].format('nan'))
        parfile.write('\n')
        parfile.write('{:>22}\t'.format('--- errmin'))
        for key in parkeys:
            outl = par.get(key + '_outlier', False)
            if outl:
                space = ' *'
            else:
                space = '  '
            parfile.write('{}{} '.format(space, key))
            try:
                err = par[key + '_err'][0]
                parfile.write(formats[key].format(err))
            except KeyError:
                parfile.write(formats_none[key].format('nan'))
        parfile.write('\n')
        parfile.write('{:>22}\t'.format('--- errmax'))
        for key in parkeys:
            outl = par.get(key + '_outlier', False)
            if outl:
                space = ' *'
            else:
                space = '  '
            parfile.write('{}{} '.format(space, key))
            try:
                err = par[key + '_err'][1]
                parfile.write(formats[key].format(err))
            except KeyError:
                parfile.write(formats_none[key].format('nan'))
        parfile.write('\n')

    means = sourcepar.means
    errors = sourcepar.errors
    means_weight = sourcepar.means_weight
    errors_weight = sourcepar.errors_weight

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

    ra_mean = means['ra']
    ra_minus, ra_plus = errors['ra']
    parfile.write('Source radius: {:.3f} /- {:.3f} /+ {:.3f} m\n'.format(
        ra_mean, ra_minus, ra_plus))
    ra_mean_weight = means_weight['ra']
    ra_minus_weight, ra_plus_weight = errors_weight['ra']
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

    if means['Ml'] is not None:
        Ml_mean = means['Ml']
        Ml_error = errors['Ml']
        parfile.write('Ml: {:.3f} +/- {:.3f} \n'.format(Ml_mean, Ml_error))

    Er_mean = means['Er']
    Er_minus, Er_plus = errors['Er']
    # format Er_plus and Er_minus to print it with the same exponent of Er
    Er_minus_str = _format_exponent(Er_minus, Er_mean)
    Er_plus_str = _format_exponent(Er_plus, Er_mean)
    parfile.write('Er: {:.3e} /- {} /+ {} N.m\n'.format(
        Er_mean, Er_minus_str, Er_plus_str))

    parfile.write('\n*** SourceSpec: {}'.format(get_versions()['version']))
    now = datetime.now()
    tz = get_localzone()
    timezone = tz.tzname(now)
    parfile.write('\n*** Run completed on: {} {}'.format(now, timezone))
    config.end_of_run = now
    config.end_of_run_tz = timezone
    if config.options.run_id:
        parfile.write('\n*** Run ID: {}'.format(config.options.run_id))
    _write_author_and_agency_to_parfile(config, parfile)

    parfile.close()

    logger.info('Output written to file: ' + parfilename)


def _log_db_write_error(db_err, db_file):
    msg = 'Unable to insert values: {}'.format(db_err)
    logger.error(msg)
    msg = 'Maybe your sqlite database has an old format.'
    logger.info(msg)
    msg = 'Try to remove or rename your database file.'
    logger.info(msg)
    msg = '(Current database file: {})'.format(db_file)
    logger.info(msg)
    ssp_exit(1)


def _write_db(config, sourcepar):
    try:
        database_file = config.database_file
    except KeyError:
        database_file = None
    if not database_file:
        return

    hypo = config.hypo
    evid = hypo.evid
    runid = config.options.run_id

    # Open SQLite database
    conn = sqlite3.connect(database_file, timeout=60)
    c = conn.cursor()

    # Init Station table
    c.execute(
        'create table if not exists Stations '
        '(stid, evid, runid,'
        'Mo, Mo_err_minus, Mo_err_plus,'
        'Mw, Mw_err_minus, Mw_err_plus,'
        'fc, fc_err_minus, fc_err_plus,'
        't_star, t_star_err_minus, t_star_err_plus,'
        'Qo, Qo_err_minus, Qo_err_plus,'
        'bsd, bsd_err_minus, bsd_err_plus,'
        'ra, ra_err_minus, ra_err_plus,'
        'dist, azimuth, Er);')
    # Write station source parameters to database
    nobs = 0
    stationpar = sourcepar.station_parameters
    for statId in sorted(stationpar.keys()):
        nobs += 1
        par = stationpar[statId]
        # Remove existing line, if present
        t = (statId, evid, runid)
        try:
            c.execute(
                'delete from Stations where stid=? and evid=? and runid=?;', t)
        except Exception as msg:
            _log_db_write_error(msg, database_file)
            ssp_exit(1)
        # Insert new line
        t = (
            statId, evid, runid,
            par['Mo'], *par['Mo_err'],
            par['Mw'], *par['Mw_err'],
            par['fc'], *par['fc_err'],
            par['t_star'], *par['t_star_err'],
            par['Qo'], *par['Qo_err'],
            par['bsd'], *par['bsd_err'],
            par['ra'], *par['ra_err'],
            par['hyp_dist'], par['az'], par['Er']
        )
        # Create a string like ?,?,?,?
        values = ','.join('?'*len(t))
        try:
            c.execute('insert into Stations values(' + values + ');', t)
        except Exception as msg:
            _log_db_write_error(msg, database_file)
            ssp_exit(1)
    # Commit changes
    conn.commit()

    # Init Event table
    c.execute(
        'create table if not exists Events '
        '(evid, runid, orig_time, lon, lat, depth, nobs,'
        'Mo, Mo_err_minus, Mo_err_plus,'
        'Mo_wavg, Mo_wavg_err_minus, Mo_wavg_err_plus,'
        'Mw, Mw_err, Mw_wavg, Mw_wavg_err,'
        'fc, fc_err_minus, fc_err_plus,'
        'fc_wavg, fc_wavg_err_minus, fc_wavg_err_plus,'
        't_star, t_star_err,'
        't_star_wavg, t_star_wavg_err,'
        'Qo, Qo_err,'
        'Qo_wavg, Qo_wavg_err,'
        'ra, ra_err_minus, ra_err_plus,'
        'ra_wavg, ra_wavg_err_minus, ra_wavg_err_plus,'
        'bsd, bsd_err_minus, bsd_err_plus,'
        'bsd_wavg, bsd_wavg_err_minus, bsd_wavg_err_plus,'
        'Er, Er_err_minus, Er_err_plus,'
        'Ml, Ml_err,'
        'run_completed, sourcespec_version,'
        'author_name, author_email,'
        'agency_full_name, agency_short_name, agency_url);')
    means = sourcepar.means
    means_weight = sourcepar.means_weight
    errors = sourcepar.errors
    errors_weight = sourcepar.errors_weight
    run_completed = '{} {}'.format(config.end_of_run, config.end_of_run_tz)
    ssp_version = get_versions()['version']
    # Remove event from Event table, if present
    t = (evid, runid)
    try:
        c.execute('delete from Events where evid=? and runid=?;', t)
    except Exception as msg:
        _log_db_write_error(msg, database_file)
        ssp_exit(1)
    t = (
        evid, runid,
        str(hypo.origin_time),
        float(hypo.longitude), float(hypo.latitude), float(hypo.depth),
        nobs,
        means['Mo'], *errors['Mo'],
        means_weight['Mo'], *errors_weight['Mo'],
        means['Mw'], errors['Mw'],
        means_weight['Mw'], errors_weight['Mw'],
        means['fc'], *errors['fc'],
        means_weight['fc'], *errors_weight['fc'],
        means['t_star'], errors['t_star'],
        means_weight['t_star'], errors_weight['t_star'],
        means['Qo'], errors['Qo'],
        means_weight['Qo'], errors_weight['Qo'],
        means['ra'], *errors['ra'],
        means_weight['ra'], *errors_weight['ra'],
        means['bsd'], *errors['bsd'],
        means_weight['bsd'], *errors_weight['bsd'],
        means['Er'], *errors['Er'],
        means['Ml'], errors['Ml'],
        run_completed, ssp_version,
        config.author_name, config.author_email,
        config.agency_full_name, config.agency_short_name,
        config.agency_url
    )
    # Create a string like ?,?,?,?
    values = ','.join('?'*len(t))
    try:
        c.execute(
            'insert into Events values(' + values + ');', t)
    except Exception as msg:
        _log_db_write_error(msg, database_file)
        ssp_exit(1)
    # Commit changes and close database
    conn.commit()
    conn.close()
    logger.info('Output written to database: ' + database_file)


def _write_hypo(config, sourcepar):
    if not config.options.hypo_file:
        return
    with open(config.options.hypo_file, 'r') as fp:
        line = fp.readline()
        # Check if first 10 digits of the line contain characters
        if any(c.isalpha() for c in line[0:10]):
            line1 = line
            line = fp.readline()
        line = list(line)

    means = sourcepar.means
    mw_str = '{:03.2f}'.format(means['Mw'])
    if means['Ml'] is not None and ~np.isnan(means['Ml']):
        ml_str = '{:03.2f}'.format(means['Ml'])
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


def write_output(config, sourcepar):
    """Write results to a plain text file and/or to a SQLite database file."""
    # Write to parfile
    _write_parfile(config, sourcepar)
    # Write to database, if requested
    _write_db(config, sourcepar)
    # Write to hypo file, if requested
    _write_hypo(config, sourcepar)
    # Write to quakeml file, if requested
    write_qml(config, sourcepar)
