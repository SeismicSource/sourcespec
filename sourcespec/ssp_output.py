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


def logmean(array):
    """Logarithmic mean of a numpy array."""
    return 10**(np.log10(array).mean())


def log_error(array):
    """Asymmetric logarithmic error bars."""
    log_array = np.log10(array)
    mean = log_array.mean()
    std = log_array.std()
    minus = logmean(array) - 10**(mean-std)
    plus = 10**(mean+std) - logmean(array)
    return minus, plus


def format_exponent(value, reference):
    """Format `value` to a string having the same exponent than `reference`."""
    # get the exponent of reference value
    xp = int(np.log10(reference))
    # format value to print it with the same exponent of reference value
    return '%5.3fe%+03d' % (value/10**xp, xp)


def write_output(config, evid, sourcepar):
    """Write results to a plain text file and/or to a SQLite database file."""
    if len(sourcepar) == 0:
        logging.info('No source parameter calculated')
        ssp_exit()

    # Optional: write station source parameters to SQLite database
    try:
        database_file = config.database_file
    except KeyError:
        database_file = None

    if database_file:
        # Open SQLite database
        conn = sqlite3.connect(database_file)
        c = conn.cursor()

        # Init database schema
        c.execute('create table if not exists Stations '
                  '(stid, evid, Mo, Mw, fc, t_star, dist, azimuth);')

        # Write station source parameters to database
        for statId in sorted(sourcepar.keys()):
            par = sourcepar[statId]

            # Remove existing line, if present
            t = (statId, evid)
            c.execute('delete from Stations where stid=? and evid=?;', t)

            # Insert new line
            t = (statId, evid, par['Mo'], par['Mw'], par['fc'], par['t_star'],
                 par['hyp_dist'], par['az'])
            c.execute('insert into Stations values(?, ?, ?, ?, ?, ?, ?, ?);',
                      t)

        # Commit changes
        conn.commit()

    # Write station source parameters to file
    if not os.path.exists(config.options.outdir):
        os.makedirs(config.options.outdir)
    parfilename = os.path.join(config.options.outdir, '%s.ssp.out' % evid)

    with open(parfilename, 'w') as parfile:
        parfile.write('*** Station source parameters ***\n')
        for statId in sorted(sourcepar.keys()):
            par = sourcepar[statId]
            parfile.write('%s\t' % statId)
            for key in par:
                if key == 'Mo':
                    parfile.write('  %s %.3e ' % (key, par[key]))
                else:
                    parfile.write('  %s %6.3f ' % (key, par[key]))
            parfile.write('\n')
        parfile.write('\n')

        # Compute and write average source parameters
        parfile.write('*** Average source parameters ***\n')
        # Mw
        Mw_array = np.array([x['Mw'] for x in sourcepar.values()])
        Mw_mean = Mw_array.mean()
        Mw_std = Mw_array.std()
        parfile.write('Mw: %.2f +/- %.2f\n' % (Mw_mean, Mw_std))

        # Mo (N.m)
        Mo_array = mag_to_moment(Mw_array)
        Mo_mean = logmean(Mo_array)
        Mo_minus, Mo_plus = log_error(Mo_array)
        # format Mo_plus and Mo_minus to print it with the same exponent of Mo
        Mo_plus_str = format_exponent(Mo_plus, Mo_mean)
        Mo_minus_str = format_exponent(Mo_minus, Mo_mean)
        parfile.write('Mo: %.3e /- %s /+ %s N.m\n' %
                      (Mo_mean, Mo_minus_str, Mo_plus_str))

        # fc , hertz
        fc_array = np.array([x['fc'] for x in sourcepar.values()])
        fc_mean = logmean(fc_array)
        fc_minus, fc_plus = log_error(fc_array)
        parfile.write('fc: %.3f /- %.3f /+ %.3f Hz\n' %
                      (fc_mean, fc_minus, fc_plus))

        # t_star
        t_star_array = np.array([x['t_star'] for x in sourcepar.values()])
        t_star_mean = t_star_array.mean()
        t_star_std = t_star_array.std()
        parfile.write('t_star: %.3f +/- %.3f Hz\n' % (t_star_mean, t_star_std))

        # ra, radius (meters)
        vs_m = config.vs*1000
        ra_array = 0.37 * vs_m / fc_array
        ra_mean = logmean(ra_array)
        ra_minus, ra_plus = log_error(ra_array)
        parfile.write('Source radius: %.3f /- %.3f /+ %.3f m\n' %
                      (ra_mean, ra_minus, ra_plus))

        # bsd, Brune stress drop (MPa)
        bsd_array = 0.4375 * Mo_array / np.power(ra_array, 3) * 1e-6
        bsd_mean = bsd_array.mean()
        bsd_std = bsd_array.std()
        # format Mo_std to print it with the same exponent of Mo
        bsd_std_str = format_exponent(bsd_std, bsd_mean)
        parfile.write('Brune stress drop: %.3e +/- %s MPa\n' %
                      (bsd_mean, bsd_std_str))

        # Ml
        Ml_array = np.array([x['Ml'] for x in sourcepar.values()])
        Ml_mean = Ml_array.mean()
        Ml_std = Ml_array.std()
        parfile.write('Ml: %.3f +/- %.3f \n' % (Ml_mean, Ml_std))

    logging.info('Output written to file: ' + parfilename)

    # Write average source parameters to database
    if database_file:
        c.execute('create table if not exists Events '
                  '(evid, Mw_mean, Mo_mean, fc_mean, t_star_mean, '
                  'ra_mean, bsd_mean, Ml_mean);')
        t = (evid, Mw_mean)
        c.execute('delete from Events where evid=? and Mw_mean=?;', t)
        t = (evid, Mw_mean, Mo_mean, fc_mean, t_star_mean, ra_mean,
             bsd_mean, Ml_mean)
        c.execute('insert into Events values(?, ?, ?, ?, ?, ?, ?, ?);', t)
        # Commit changes and close database
        conn.commit()
        conn.close()
        logging.info('Output written to database: ' + database_file)

    if config.options.hypo_file:
        with open(config.options.hypo_file, 'r') as fp:
            line = fp.readline()
            # Check if first 10 digits of the line contain characters
            if any(c.isalpha() for c in line[0:10]):
                line1 = line
                line = fp.readline()
            line = list(line)
        mw_str = '%03.2f' % Mw_mean
        ml_str = '%03.2f' % Ml_mean
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

    params_name = ('Mw', 'fc', 't_star')
    sourcepar_mean = dict(zip(params_name, [Mw_mean, fc_mean, t_star_mean]))
    logging.info('params_mean: %s' % (sourcepar_mean))

    return sourcepar_mean
