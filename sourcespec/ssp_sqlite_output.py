# -*- coding: utf-8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
SQLite output for source_spec.

:copyright:
    2013-2022 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
import os.path
import logging
import sqlite3
from sourcespec.ssp_setup import ssp_exit
from sourcespec._version import get_versions
logger = logging.getLogger(__name__.split('.')[-1])


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


def write_sqlite(config, sspec_output):
    try:
        database_file = config.database_file
    except KeyError:
        database_file = None
    if not database_file:
        return

    hypo = config.hypo
    evid = hypo.evid
    runid = config.options.run_id

    # Current supported DB version
    DB_VERSION = 1

    # only check version if database_file exsists
    check_version = os.path.isfile(database_file)

    # Open SQLite database
    try:
        conn = sqlite3.connect(database_file, timeout=60)
    except Exception as msg:
        logger.error(msg)
        logger.info(
            'Please check whether "{}" is a valid SQLite file.'.format(
                database_file)
        )
        ssp_exit(1)
    c = conn.cursor()

    if check_version:
        # Get current DB version
        db_version = c.execute('PRAGMA user_version').fetchone()[0]
        if db_version < DB_VERSION:
            msg = '"{}" has an old database version: "{}". '.format(
                database_file, db_version)
            msg += 'Current supported version is "{}".'.format(DB_VERSION)
            logger.error(msg)
            msg = 'Remove or rename your old database file, '
            msg += 'so that a new one can be created.'
            logger.info(msg)
            exit(1)
    else:
        # Set the DB version
        c.execute('PRAGMA user_version = {v:d}'.format(v=DB_VERSION))

    # Create Station table
    sql_create_stations_table = """CREATE TABLE IF NOT EXISTS Stations (
        stid TEXT PRIMARY KEY,
        evid TEXT,
        runid TEXT,
        Mo REAL,
        Mo_err_minus REAL,
        Mo_err_plus REAL,
        Mw REAL,
        Mw_err_minus REAL,
        Mw_err_plus REAL,
        fc REAL,
        fc_err_minus REAL,
        fc_err_plus REAL,
        t_star REAL,
        t_star_err_minus REAL,
        t_star_err_plus REAL,
        Qo REAL,
        Qo_err_minus REAL,
        Qo_err_plus REAL,
        bsd REAL,
        bsd_err_minus REAL,
        bsd_err_plus REAL,
        ra REAL,
        ra_err_minus REAL,
        ra_err_plus REAL,
        dist REAL,
        azimuth REAL,
        Er REAL
    );"""
    c.execute(sql_create_stations_table)
    # Write station source parameters to database
    nobs = 0
    stationpar = sspec_output.station_parameters
    for statId in sorted(stationpar.keys()):
        nobs += 1
        par = stationpar[statId]
        # Remove existing line, if present
        t = (statId, evid, runid)
        sql_delete_from_sations =\
            'DELETE FROM Stations WHERE stid=? AND evid=? AND runid=?;'
        try:
            c.execute(sql_delete_from_sations, t)
        except Exception as msg:
            _log_db_write_error(msg, database_file)
            ssp_exit(1)
        # Insert new line
        t = (
            statId, evid, runid,
            *par.Mo.value_uncertainty(),
            *par.Mw.value_uncertainty(),
            *par.fc.value_uncertainty(),
            *par.t_star.value_uncertainty(),
            *par.Qo.value_uncertainty(),
            *par.bsd.value_uncertainty(),
            *par.radius.value_uncertainty(),
            par.hypo_dist_in_km,
            par.azimuth,
            par.Er.value
        )
        # Create a string like ?,?,?,?
        values = ','.join('?'*len(t))
        sql_insert_into_stations =\
            'INSERT INTO Stations VALUES({});'.format(values)
        try:
            c.execute(sql_insert_into_stations, t)
        except Exception as msg:
            _log_db_write_error(msg, database_file)
            ssp_exit(1)
    # Commit changes
    conn.commit()

    # Create Event table
    sql_create_events_table = """CREATE TABLE IF NOT EXISTS Events (
        /* Event info */
            evid TEXT PRIMARY KEY,
            runid TEXT,
            orig_time REAL,
            lon REAL,
            lat REAL,
            depth REAL,
        /* Statistical info */
            nobs INTEGER,
            nsigma REAL,
            mid_pct REAL,
            lower_pct REAL,
            upper_pct REAL,
        /* Seismic moment */
            Mo_mean REAL,
            Mo_mean_err_minus REAL,
            Mo_mean_err_plus REAL,
            Mo_wmean REAL,
            Mo_wmean_err_minus REAL,
            Mo_wmean_err_plus REAL,
            Mo_pctl REAL,
            Mo_pctl_err_minus REAL,
            Mo_pctl_err_plus REAL,
        /* Moment magnitude */
            Mw_mean REAL,
            Mw_mean_err_minus REAL,
            Mw_mean_err_plus REAL,
            Mw_wmean REAL,
            Mw_wmean_err_minus REAL,
            Mw_wmean_err_plus REAL,
            Mw_pctl REAL,
            Mw_pctl_err_minus REAL,
            Mw_pctl_err_plus REAL,
        /* Corner frequency */
            fc_mean REAL,
            fc_mean_err_minus REAL,
            fc_mean_err_plus REAL,
            fc_wmean REAL,
            fc_wmean_err_minus REAL,
            fc_wmean_err_plus REAL,
            fc_pctl REAL,
            fc_pctl_err_minus REAL,
            fc_pctl_err_plus REAL,
        /* t-star */
            t_star_mean REAL,
            t_star_mean_err_minus REAL,
            t_star_mean_err_plus REAL,
            t_star_wmean REAL,
            t_star_wmean_err_minus REAL,
            t_star_wmean_err_plus REAL,
            t_star_pctl REAL,
            t_star_pctl_err_minus REAL,
            t_star_pctl_err_plus REAL,
        /* Qo */
            Qo_mean REAL,
            Qo_mean_err_minus REAL,
            Qo_mean_err_plus REAL,
            Qo_wmean REAL,
            Qo_wmean_err_minus REAL,
            Qo_wmean_err_plus REAL,
            Qo_pctl REAL,
            Qo_pctl_err_minus REAL,
            Qo_pctl_err_plus REAL,
        /* Source radius */
            ra_mean REAL,
            ra_mean_err_minus REAL,
            ra_mean_err_plus REAL,
            ra_wmean REAL,
            ra_wmean_err_minus REAL,
            ra_wmean_err_plus REAL,
            ra_pctl REAL,
            ra_pctl_err_minus REAL,
            ra_pctl_err_plus REAL,
        /* Brune stress drop */
            bsd_mean REAL,
            bsd_mean_err_minus REAL,
            bsd_mean_err_plus REAL,
            bsd_wmean REAL,
            bsd_wmean_err_minus REAL,
            bsd_wmean_err_plus REAL,
            bsd_pctl REAL,
            bsd_pctl_err_minus REAL,
            bsd_pctl_err_plus REAL,
        /* Radiated energy */
            Er_mean REAL,
            Er_mean_err_minus REAL,
            Er_mean_err_plus REAL,
            Er_pctl REAL,
            Er_pctl_err_minus REAL,
            Er_pctl_err_plus REAL,
        /* Local magnitude */
            Ml_mean REAL,
            Ml_mean_err_minus REAL,
            Ml_mean_err_plus REAL,
            Ml_pctl REAL,
            Ml_pctl_err_minus REAL,
            Ml_pctl_err_plus REAL,
        /* Run info */
            run_completed TEXT,
            sourcespec_version TEXT,
            author_name TEXT,
            author_email TEXT,
            agency_full_name TEXT,
            agency_short_name TEXT,
            agency_url TEXT
    );"""
    c.execute(sql_create_events_table)
    means = sspec_output.mean_values()
    mean_errors = sspec_output.mean_uncertainties()
    wmeans = sspec_output.weighted_mean_values()
    wmean_errors = sspec_output.weighted_mean_uncertainties()
    percentiles = sspec_output.percentiles_values()
    percentile_errors = sspec_output.percentiles_uncertainties()
    run_completed = '{} {}'.format(config.end_of_run, config.end_of_run_tz)
    ssp_version = get_versions()['version']
    # Remove event from Event table, if present
    t = (evid, runid)
    sql_delete_from_events = 'DELETE FROM Events WHERE evid=? AND runid=?;'
    try:
        c.execute(sql_delete_from_events, t)
    except Exception as msg:
        _log_db_write_error(msg, database_file)
        ssp_exit(1)
    t = (
        # Event info
        evid,
        runid,
        str(hypo.origin_time),
        float(hypo.longitude),
        float(hypo.latitude),
        float(hypo.depth),
        # Statistical info
        nobs,
        config.n_sigma,
        config.lower_percentage,
        config.mid_percentage,
        config.upper_percentage,
        # Seismic moment
        means['Mo'],
        *mean_errors['Mo'],
        wmeans['Mo'],
        *wmean_errors['Mo'],
        percentiles['Mo'],
        *percentile_errors['Mo'],
        # Moment magnitude
        means['Mw'],
        *mean_errors['Mw'],
        wmeans['Mw'],
        *wmean_errors['Mw'],
        percentiles['Mw'],
        *percentile_errors['Mw'],
        # Corner frequency
        means['fc'],
        *mean_errors['fc'],
        wmeans['fc'],
        *wmean_errors['fc'],
        percentiles['fc'],
        *percentile_errors['fc'],
        # t-star
        means['t_star'],
        *mean_errors['t_star'],
        wmeans['t_star'],
        *wmean_errors['t_star'],
        percentiles['t_star'],
        *percentile_errors['t_star'],
        # Qo
        means['Qo'],
        *mean_errors['Qo'],
        wmeans['Qo'],
        *wmean_errors['Qo'],
        percentiles['Qo'],
        *percentile_errors['Qo'],
        # Source radius
        means['radius'],
        *mean_errors['radius'],
        wmeans['radius'],
        *wmean_errors['radius'],
        percentiles['radius'],
        *percentile_errors['radius'],
        # Brune stress drop
        means['bsd'],
        *mean_errors['bsd'],
        wmeans['bsd'],
        *wmean_errors['bsd'],
        percentiles['bsd'],
        *percentile_errors['bsd'],
        # Radiated energy
        means['Er'],
        *mean_errors['Er'],
        percentiles['Er'],
        *percentile_errors['Er'],
        # Local magnitude
        means.get('Ml', None),
        *mean_errors.get('Ml', (None, None)),
        percentiles.get('Ml', None),
        *percentile_errors.get('Ml', (None, None)),
        # Run info
        run_completed,
        ssp_version,
        config.author_name,
        config.author_email,
        config.agency_full_name,
        config.agency_short_name,
        config.agency_url
    )
    # Create a string like ?,?,?,?
    values = ','.join('?'*len(t))
    sql_insert_into_events = 'INSERT INTO Events VALUES({});'.format(values)
    try:
        c.execute(sql_insert_into_events, t)
    except Exception as msg:
        _log_db_write_error(msg, database_file)
        ssp_exit(1)
    # Commit changes and close database
    conn.commit()
    conn.close()
    logger.info('Output written to SQLite database: ' + database_file)
