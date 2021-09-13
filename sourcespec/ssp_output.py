# -*- coding: utf-8 -*-
"""
Output functions for source_spec.

:copyright:
    2012 Claudio Satriano <satriano@ipgp.fr>

    2013-2014 Claudio Satriano <satriano@ipgp.fr>,
              Emanuela Matrullo <matrullo@geologie.ens.fr>

    2015-2021 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement, Version 2.1
    (http://www.cecill.info/index.en.html)
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import os
import logging
import sqlite3
import numpy as np
from datetime import datetime
from pytz import reference
from sourcespec.ssp_setup import ssp_exit
from sourcespec.ssp_util import mag_to_moment
logger = logging.getLogger(__name__.split('.')[-1])


def _avg_and_std(values, errors=None, logarithmic=False):
    """
    Return the average and standard deviation.

    Optionally:
    - errors can be specfied for weighted statistics
    - logarithmic average and standard deviation
    """
    if logarithmic:
        values = np.log10(values)
    if errors is None:
        weights = None
    else:
        weights = 1./(errors**2.)
    average = np.average(values, weights=weights)
    variance = np.average((values-average)**2, weights=weights)
    std = np.sqrt(variance)
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
            'Mw', 'fc', 't_star', 'Mo', 'hyp_dist', 'az', 'Er'
        )
        formats = dict(
            Mo='  {} {:.3e} ',
            Er='  {} {:.3e} ',
            hyp_dist='  {} {:7.3f} ',
            az='  {} {:7.3f} ',
            Mw='  {} {:6.3f} ',
            fc='  {} {:6.3f} ',
            t_star='  {} {:6.3f} ',
            Ml='  {} {:6.3f} '
        )
        formats_none = dict(
            Mo='  {} {:>9} ',
            Er='  {} {:>9} ',
            hyp_dist='  {} {:>7} ',
            az='  {} {:>7} ',
            Mw='  {} {:>6} ',
            fc='  {} {:>6} ',
            t_star='  {} {:>6} ',
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
            parfile.write('{:>21}\t'.format('--- errors'))
            for key in parkeys:
                try:
                    val = err[key]
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

        ra_mean = means['ra']
        ra_minus, ra_plus = errors['ra']
        parfile.write('Source radius: {:.3f} /- {:.3f} /+ {:.3f} m\n'.format(
            ra_mean, ra_minus, ra_plus))

        bsd_mean = means['bsd']
        bsd_minus, bsd_plus = errors['bsd']
        bsd_minus_str = _format_exponent(bsd_minus, bsd_mean)
        bsd_plus_str = _format_exponent(bsd_plus, bsd_mean)
        parfile.write('Brune stress drop: {:.3e} /- {} /+ {} MPa\n'.format(
            bsd_mean, bsd_minus_str, bsd_plus_str))

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
        localtime = reference.LocalTimezone()
        timezone = localtime.tzname(now)
        parfile.write('\n*** Run completed on: {} {}\n'.format(now, timezone))
        config.end_of_run = now
        config.end_of_run_tz = timezone

    logger.info('Output written to file: ' + parfilename)


def _write_db(config, sourcepar):
    try:
        database_file = config.database_file
    except KeyError:
        database_file = None
    if not database_file:
        return

    evid = config.hypo.evid

    # Open SQLite database
    conn = sqlite3.connect(database_file)
    c = conn.cursor()

    # Init database schema
    c.execute('create table if not exists Stations '
              '(stid, evid, Mo, Mw, fc, t_star, dist, azimuth);')

    # Write station source parameters to database
    for statId in sorted(sourcepar.keys()):
        if statId in ['means', 'errors', 'means_weight', 'errors_weight']:
            continue
        par = sourcepar[statId]

        # Remove existing line, if present
        t = (statId, evid)
        c.execute('delete from Stations where stid=? and evid=?;', t)

        # Insert new line
        t = (statId, evid, par['Mo'], par['Mw'], par['fc'], par['t_star'],
             par['hyp_dist'], par['az'])
        c.execute('insert into Stations values(?, ?, ?, ?, ?, ?, ?, ?);', t)

    # Commit changes
    conn.commit()

    means = sourcepar['means']

    c.execute('create table if not exists Events '
              '(evid, Mo_mean, Mw_mean, fc_mean, t_star_mean, '
              'ra_mean, bsd_mean, Ml_mean);')
    t = (evid, means['Mw'])
    c.execute('delete from Events where evid=? and Mw_mean=?;', t)
    t = (evid, means['Mo'], means['Mw'], means['fc'], means['t_star'],
         means['ra'], means['bsd'], means['Ml'])
    c.execute('insert into Events values(?, ?, ?, ?, ?, ?, ?, ?);', t)
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
    vs_m = config.hypo.vs*1e3
    # Madariaga (2009), doi:10.1007/978-1-4419-7695-6_22, eq. 31:
    ra_values = 0.3724 * vs_m / fc_values
    means['ra'], errors['ra'] = _avg_and_std(ra_values, logarithmic=True)

    # bsd, Brune stress drop (MPa)
    Mo_values = mag_to_moment(Mw_values)
    # Madariaga (2009), doi:10.1007/978-1-4419-7695-6_22, eq. 27:
    bsd_values = 7./16 * Mo_values / np.power(ra_values, 3) * 1e-6
    means['bsd'], errors['bsd'] = _avg_and_std(bsd_values, logarithmic=True)

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
    _write_db(config, sourcepar)

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
