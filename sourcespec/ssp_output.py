# -*- coding: utf-8 -*-
"""
Output functions for source_spec.

:copyright:
    2012 Claudio Satriano <satriano@ipgp.fr>

    2013-2014 Claudio Satriano <satriano@ipgp.fr>,
              Emanuela Matrullo <matrullo@geologie.ens.fr>

    2015-2016 Claudio Satriano <satriano@ipgp.fr>
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
from sourcespec.ssp_setup import ssp_exit
from sourcespec.ssp_util import mag_to_moment


def _logmean(array):
    """Logarithmic mean of a numpy array."""
    return 10**(np.log10(array).mean())


def _log_error(array):
    """Asymmetric logarithmic error bars."""
    log_array = np.log10(array)
    mean = log_array.mean()
    std = log_array.std()
    minus = _logmean(array) - 10**(mean-std)
    plus = 10**(mean+std) - _logmean(array)
    return minus, plus


def _format_exponent(value, reference):
    """Format `value` to a string having the same exponent than `reference`."""
    # get the exponent of reference value
    xp = int(np.log10(reference))
    # format value to print it with the same exponent of reference value
    return '%5.3fe%+03d' % (value/10**xp, xp)


def _write_parfile(config, evid, sourcepar):
    """Write station source parameters to file."""
    if not os.path.exists(config.options.outdir):
        os.makedirs(config.options.outdir)
    parfilename = os.path.join(config.options.outdir, '%s.ssp.out' % evid)

    with open(parfilename, 'w') as parfile:
        parfile.write('*** Station source parameters ***\n')
        for statId in sorted(sourcepar.keys()):
            if statId in ['means', 'errors']:
                continue
            par = sourcepar[statId]
            parfile.write('%s\t' % statId)
            for key in par:
                if key in ['Mo', 'Er']:
                    parfile.write('  %s %.3e ' % (key, par[key]))
                else:
                    if par[key] is not None:
                        parfile.write('  %s %6.3f ' % (key, par[key]))
                    else:
                        parfile.write('  %s %s ' % (key, par[key]))
            parfile.write('\n')

        means = sourcepar['means']
        errors = sourcepar['errors']

        parfile.write('\n*** Average source parameters ***\n')

        Mw_mean = means['Mw']
        Mw_error = errors['Mw']
        parfile.write('Mw: %.2f +/- %.2f\n' % (Mw_mean, Mw_error))

        Mo_mean = means['Mo']
        Mo_minus, Mo_plus = errors['Mo']
        # format Mo_plus and Mo_minus to print it with the same exponent of Mo
        Mo_minus_str = _format_exponent(Mo_minus, Mo_mean)
        Mo_plus_str = _format_exponent(Mo_plus, Mo_mean)
        parfile.write('Mo: %.3e /- %s /+ %s N.m\n' %
                      (Mo_mean, Mo_minus_str, Mo_plus_str))

        fc_mean = means['fc']
        fc_minus, fc_plus = errors['fc']
        parfile.write('fc: %.3f /- %.3f /+ %.3f Hz\n' %
                      (fc_mean, fc_minus, fc_plus))

        t_star_mean = means['t_star']
        t_star_error = errors['t_star']
        parfile.write('t_star: %.3f +/- %.3f Hz\n' %
                      (t_star_mean, t_star_error))

        ra_mean = means['ra']
        ra_minus, ra_plus = errors['ra']
        parfile.write('Source radius: %.3f /- %.3f /+ %.3f m\n' %
                      (ra_mean, ra_minus, ra_plus))

        bsd_mean = means['bsd']
        bsd_error = errors['bsd']
        # format Mo_std to print it with the same exponent of Mo
        bsd_error_str = _format_exponent(bsd_error, bsd_mean)
        parfile.write('Brune stress drop: %.3e +/- %s MPa\n' %
                      (bsd_mean, bsd_error_str))

        if means['Ml'] is not None:
            Ml_mean = means['Ml']
            Ml_error = errors['Ml']
            parfile.write('Ml: %.3f +/- %.3f \n' % (Ml_mean, Ml_error))

        Er_mean = means['Er']
        Er_minus, Er_plus = errors['Er']
        # format Er_plus and Er_minus to print it with the same exponent of Er
        Er_minus_str = _format_exponent(Er_minus, Er_mean)
        Er_plus_str = _format_exponent(Er_plus, Er_mean)
        parfile.write('Er: %.3e /- %s /+ %s N.m\n' %
                      (Er_mean, Er_minus_str, Er_plus_str))

    logging.info('Output written to file: ' + parfilename)


def _write_db(config, evid, sourcepar):
    try:
        database_file = config.database_file
    except KeyError:
        database_file = None
    if not database_file:
        return

    # Open SQLite database
    conn = sqlite3.connect(database_file)
    c = conn.cursor()

    # Init database schema
    c.execute('create table if not exists Stations '
              '(stid, evid, Mo, Mw, fc, t_star, dist, azimuth);')

    # Write station source parameters to database
    for statId in sorted(sourcepar.keys()):
        if statId in ['means', 'errors']:
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
    logging.info('Output written to database: ' + database_file)


def _write_hypo(config, evid, sourcepar):
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
    mw_str = '%03.2f' % means['Mw']
    if means['Ml'] is not None:
        ml_str = '%03.2f' % means['Ml']
    else:
        ml_str = ' '*4
    for i in range(0, 4):
        line[49+i] = mw_str[0+i]
        #line[45+i] = mw_str[0+i]
        line[69+i] = ml_str[0+i]
    outline = ''.join(line)
    hypo_file_out = os.path.join(config.options.outdir, '%s.ssp.h' % evid)
    with open(hypo_file_out, 'w') as fp:
        try:
            fp.write(line1)
        except:
            pass
        fp.write(outline)
    logging.info('Hypo file written to: ' + hypo_file_out)


def write_output(config, evid, sourcepar, sourcepar_err):
    """Write results to a plain text file and/or to a SQLite database file."""
    if len(sourcepar) == 0:
        logging.info('No source parameter calculated')
        ssp_exit()

    means = dict()
    errors = dict()

    # Compute average source parameters
    # Mw
    Mw_array = np.array([x['Mw'] for x in sourcepar.values()])
    means['Mw'] = Mw_array.mean()
    errors['Mw'] = Mw_array.std()

    # Mo (N.m)
    Mo_array = mag_to_moment(Mw_array)
    means['Mo'] = _logmean(Mo_array)
    errors['Mo'] = _log_error(Mo_array)

    # fc , hertz
    fc_array = np.array([x['fc'] for x in sourcepar.values()])
    means['fc'] = _logmean(fc_array)
    errors['fc'] = _log_error(fc_array)

    # t_star
    t_star_array = np.array([x['t_star'] for x in sourcepar.values()])
    means['t_star'] = t_star_array.mean()
    errors['t_star'] = t_star_array.std()

    # ra, radius (meters)
    vs_m = config.hypo.vs*1000
    ra_array = 0.37 * vs_m / fc_array
    means['ra'] = _logmean(ra_array)
    errors['ra'] = _log_error(ra_array)

    # bsd, Brune stress drop (MPa)
    bsd_array = 0.4375 * Mo_array / np.power(ra_array, 3) * 1e-6
    means['bsd'] = bsd_array.mean()
    errors['bsd'] = bsd_array.std()

    # Ml
    Ml_array = np.array([x['Ml'] for x in sourcepar.values()
                         if x['Ml'] is not None])
    if Ml_array.size:
        means['Ml'] = Ml_array.mean()
        errors['Ml'] = Ml_array.std()
    else:
        means['Ml'] = None
        errors['Ml'] = None

    # Er
    Er_array = np.array([x['Er'] for x in sourcepar.values()])
    means['Er'] = _logmean(Er_array)
    errors['Er'] = _log_error(Er_array)

    sourcepar['means'] = means
    sourcepar['errors'] = errors

    # Write to parfile
    _write_parfile(config, evid, sourcepar)

    # Write to database, if requested
    _write_db(config, evid, sourcepar)

    # Write to hypo file, if requested
    _write_hypo(config, evid, sourcepar)

    params_name = ('Mw', 'fc', 't_star')
    sourcepar_mean = dict(zip(params_name,
                              [means['Mw'], means['fc'], means['t_star']]))
    logging.info('params_mean: %s' % (sourcepar_mean))

    return sourcepar_mean
