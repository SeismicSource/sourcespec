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
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import os
import logging
import sqlite3
import numpy as np
from datetime import datetime
from tzlocal import get_localzone
from sourcespec.ssp_setup import ssp_exit
from sourcespec.ssp_util import mag_to_moment
logger = logging.getLogger(__name__.split('.')[-1])


def _avg_and_std(values, errors=None, logarithmic=False, std_cutoff=True):
    """
    Return the average and standard deviation.

    Optionally:
    - errors can be specfied for weighted statistics
    - logarithmic average and standard deviation
    - cutoff outlier values that are larger than 2*std
    """
    if len(values) == 0:
        return np.nan, np.nan
    if errors is not None and len(errors) == 0:
        return np.nan, np.nan
    if logarithmic:
        values = np.log10(values)
    if errors is None:
        weights = None
    else:
        # compute a symmetric error
        errors = (errors[:, 0] + errors[:, 1])/2
        # fix for infinite weight (zero error)
        errors[errors == 0] = np.nanmin(errors[errors > 0])
        if logarithmic:
            # use relative errors in logarithmic statistics
            errors /= values
        weights = 1./(errors**2.)
    average = np.average(values, weights=weights)
    variance = np.average((values-average)**2, weights=weights)
    std = np.sqrt(variance)
    if std_cutoff and std > 0:
        nstd = 2
        for _ in range(10):
            condition = np.abs(values-average) < nstd*std
            # break if condition is always true: no reason to go on
            if np.all(condition):
                break
            values = values[condition]
            # break if less than three values: no std can be computed
            if len(values) < 3:
                break
            if weights is not None:
                weights = weights[condition]
            average = np.average(values, weights=weights)
            variance = np.average((values-average)**2, weights=weights)
            std = np.sqrt(variance)
            # raise nstd for subsequent iterations
            nstd = 2.5
    if logarithmic:
        log_average = 10.**average
        minus = log_average - 10.**(average-std)
        plus = 10.**(average+std) - log_average
        return log_average, (minus, plus)
    else:
        return average, std


def _M0_avg_and_std(Mw_mean, Mw_error):
    """Compute average and standard deviation for Mo starting from Mw."""
    Mo_mean = mag_to_moment(Mw_mean)
    Mo_min = mag_to_moment(Mw_mean - Mw_error)
    Mo_max = mag_to_moment(Mw_mean + Mw_error)
    return Mo_mean, (Mo_mean - Mo_min, Mo_max - Mo_mean)


def _format_exponent(value, reference):
    """Format `value` to a string having the same exponent than `reference`."""
    # get the exponent of reference value
    xp = np.int(np.floor(np.log10(np.abs(reference))))
    # format value to print it with the same exponent of reference value
    return '{:5.3f}e{:+03d}'.format(value/10**xp, xp)


def _write_parfile(config, sourcepar, sourcepar_err):
    """Write station source parameters to file."""
    if not os.path.exists(config.options.outdir):
        os.makedirs(config.options.outdir)
    evid = config.hypo.evid
    parfilename = os.path.join(
        config.options.outdir, '{}.ssp.out'.format(evid))

    with open(parfilename, 'w') as parfile:
        hypo = config.hypo
        parfile.write(
            '{} lon {:8.3f} lat {:7.3f} depth {:5.1f} km '
            'orig_time {}\n\n'.format(
                hypo.evid, hypo.longitude, hypo.latitude, hypo.depth,
                hypo.origin_time))
        parfile.write('*** Station source parameters ***\n')
        parkeys = (
            'Mw', 'fc', 't_star', 'Qo', 'Mo',
            'bsd', 'ra', 'hyp_dist', 'az', 'Er'
        )
        formats = dict(
            Mo='  {} {:.3e} ',
            Er='  {} {:.3e} ',
            hyp_dist='  {} {:7.3f} ',
            az='  {} {:7.3f} ',
            Mw='  {} {:6.3f} ',
            fc='  {} {:6.3f} ',
            bsd='  {} {:.3e} ',
            ra='  {} {:7.3f} ',
            t_star='  {} {:6.3f} ',
            Qo='  {} {:5.1f} ',
            Ml='  {} {:6.3f} '
        )
        formats_none = dict(
            Mo='  {} {:>9} ',
            Er='  {} {:>9} ',
            hyp_dist='  {} {:>7} ',
            az='  {} {:>7} ',
            Mw='  {} {:>6} ',
            fc='  {} {:>6} ',
            bsd='  {} {:>9} ',
            ra='  {} {:>7} ',
            t_star='  {} {:>6} ',
            Qo='  {} {:>5} ',
            Ml='  {} {:>6} '
        )
        for statId in sorted(sourcepar.keys()):
            if statId in ['means', 'errors', 'means_weight', 'errors_weight']:
                continue
            par = sourcepar[statId]
            parfile.write('{:>14} {:>6}\t'.format(*statId.split()))
            for key in parkeys:
                val = par[key]
                if val is not None:
                    parfile.write(formats[key].format(key, val))
                else:
                    parfile.write(formats_none[key].format(key, 'nan'))
            parfile.write('\n')
            err = sourcepar_err[statId]
            parfile.write('{:>21}\t'.format('--- errmin'))
            for key in parkeys:
                try:
                    val = err[key][0]
                    parfile.write(formats[key].format(key, val))
                except KeyError:
                    parfile.write(formats_none[key].format(key, 'nan'))
            parfile.write('\n')
            parfile.write('{:>21}\t'.format('--- errmax'))
            for key in parkeys:
                try:
                    val = err[key][1]
                    parfile.write(formats[key].format(key, val))
                except KeyError:
                    parfile.write(formats_none[key].format(key, 'nan'))
            parfile.write('\n')

        means = sourcepar['means']
        errors = sourcepar['errors']
        means_weight = sourcepar['means_weight']
        errors_weight = sourcepar['errors_weight']

        parfile.write('\n*** Average source parameters ***\n')

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

        now = datetime.now()
        tz = get_localzone()
        timezone = tz.tzname(now)
        parfile.write('\n*** Run completed on: {} {}\n'.format(now, timezone))
        config.end_of_run = now
        config.end_of_run_tz = timezone

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


def _write_db(config, sourcepar, sourcepar_err):
    try:
        database_file = config.database_file
    except KeyError:
        database_file = None
    if not database_file:
        return

    hypo = config.hypo
    evid = hypo.evid

    # Open SQLite database
    conn = sqlite3.connect(database_file, timeout=60)
    c = conn.cursor()

    # Init Station table
    c.execute(
        'create table if not exists Stations '
        '(stid, evid,'
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
    for statId in sorted(sourcepar.keys()):
        if statId in ['means', 'errors', 'means_weight', 'errors_weight']:
            continue
        nobs += 1
        par = sourcepar[statId]
        par_err = sourcepar_err[statId]
        # Remove existing line, if present
        t = (statId, evid)
        c.execute('delete from Stations where stid=? and evid=?;', t)
        # Insert new line
        t = (
            statId, evid,
            par['Mo'], *par_err['Mo'],
            par['Mw'], *par_err['Mw'],
            par['fc'], *par_err['fc'],
            par['t_star'], *par_err['t_star'],
            par['Qo'], *par_err['Qo'],
            par['bsd'], *par_err['bsd'],
            par['ra'], *par_err['ra'],
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
        '(evid, orig_time, lon, lat, depth, nobs,'
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
        'Ml, Ml_err);')
    means = sourcepar['means']
    means_weight = sourcepar['means_weight']
    errors = sourcepar['errors']
    errors_weight = sourcepar['errors_weight']
    # Remove event from Event table, if present
    t = (evid, )
    c.execute('delete from Events where evid=?;', t)
    t = (
        evid,
        str(hypo.origin_time), hypo.longitude, hypo.latitude, hypo.depth,
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
        means['Ml'], errors['Ml']
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

    means = sourcepar['means']
    mw_str = '{:03.2f}'.format(means['Mw'])
    if means['Ml'] is not None:
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


def write_output(config, sourcepar, sourcepar_err):
    """Write results to a plain text file and/or to a SQLite database file."""
    if len(sourcepar) == 0:
        logger.info('No source parameter calculated')
        ssp_exit()

    means = dict()
    errors = dict()
    means_weight = dict()
    errors_weight = dict()

    # Compute average source parameters
    # Mw
    Mw_values = np.array([x['Mw'] for x in sourcepar.values()])
    Mw_err = np.array([x['Mw'] for x in sourcepar_err.values()])
    means['Mw'], errors['Mw'] = _avg_and_std(Mw_values)
    means_weight['Mw'], errors_weight['Mw'] = \
        _avg_and_std(Mw_values, errors=Mw_err)

    # Mo (N.m)
    means['Mo'], errors['Mo'] = _M0_avg_and_std(means['Mw'], errors['Mw'])
    means_weight['Mo'], errors_weight['Mo'] = \
        _M0_avg_and_std(means_weight['Mw'], errors_weight['Mw'])

    # fc , hertz
    fc_values = np.array([x['fc'] for x in sourcepar.values()])
    fc_err = np.array([x['fc'] for x in sourcepar_err.values()])
    means['fc'], errors['fc'] = _avg_and_std(fc_values, logarithmic=True)
    means_weight['fc'], errors_weight['fc'] = \
        _avg_and_std(fc_values, errors=fc_err, logarithmic=True)

    # t_star
    t_star_values = np.array([x['t_star'] for x in sourcepar.values()])
    t_star_err = np.array([x['t_star'] for x in sourcepar_err.values()])
    means['t_star'], errors['t_star'] = _avg_and_std(t_star_values)
    means_weight['t_star'], errors_weight['t_star'] = \
        _avg_and_std(t_star_values, errors=t_star_err)

    # ra, radius (meters)
    ra_values = np.array([x['ra'] for x in sourcepar.values()])
    ra_err = np.array([x['ra'] for x in sourcepar_err.values()])
    means['ra'], errors['ra'] = _avg_and_std(ra_values, logarithmic=True)
    means_weight['ra'], errors_weight['ra'] = \
        _avg_and_std(ra_values, errors=ra_err, logarithmic=True)

    # bsd, Brune stress drop (MPa)
    bsd_values = np.array([x['bsd'] for x in sourcepar.values()])
    bsd_err = np.array([x['bsd'] for x in sourcepar_err.values()])
    means['bsd'], errors['bsd'] = _avg_and_std(bsd_values, logarithmic=True)
    means_weight['bsd'], errors_weight['bsd'] = \
        _avg_and_std(bsd_values, errors=bsd_err, logarithmic=True)

    # Quality factor
    Qo_values = np.array([x['Qo'] for x in sourcepar.values()])
    Qo_err = np.array([x['Qo'] for x in sourcepar_err.values()])
    cond = ~np.isinf(Qo_values)
    means['Qo'], errors['Qo'] = _avg_and_std(Qo_values[cond])
    cond_err = np.logical_and(~np.isinf(Qo_err[:, 0]), ~np.isinf(Qo_err[:, 1]))
    cond = np.logical_and(cond, cond_err)
    means_weight['Qo'], errors_weight['Qo'] = \
        _avg_and_std(Qo_values[cond], errors=Qo_err[cond])

    # Ml
    # build Ml_values: use np.nan for missing values
    Ml_values = np.array([x.get('Ml', np.nan) for x in sourcepar.values()])
    Ml_values = Ml_values[~np.isnan(Ml_values)]
    if Ml_values.size:
        means['Ml'], errors['Ml'] = _avg_and_std(Ml_values)
    else:
        means['Ml'] = None
        errors['Ml'] = None

    # Er
    Er_values = np.array([x['Er'] for x in sourcepar.values()])
    means['Er'], errors['Er'] = _avg_and_std(Er_values, logarithmic=True)

    sourcepar['means'] = means
    sourcepar['errors'] = errors
    sourcepar['means_weight'] = means_weight
    sourcepar['errors_weight'] = errors_weight

    # Write to parfile
    _write_parfile(config, sourcepar, sourcepar_err)

    # Write to database, if requested
    _write_db(config, sourcepar, sourcepar_err)

    # Write to hypo file, if requested
    _write_hypo(config, sourcepar)

    params_name = ('Mw', 'fc', 't_star')
    sourcepar_mean = dict(
        zip(params_name, [means['Mw'], means['fc'], means['t_star']]))
    logger.info('params_mean: {}'.format(sourcepar_mean))
    sourcepar_mean_weight = dict(
        zip(params_name,
            [means_weight['Mw'], means_weight['fc'], means_weight['t_star']]))
    logger.info('params_mean_weighted: {}'.format(sourcepar_mean_weight))

    return sourcepar_mean
