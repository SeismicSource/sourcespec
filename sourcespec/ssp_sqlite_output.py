# -*- coding: utf-8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
SQLite output for source_spec.

:copyright:
    2013-2023 Claudio Satriano <satriano@ipgp.fr>
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
    logger.error(f'Unable to insert values: {db_err}')
    logger.info('Maybe your sqlite database has an old format.')
    logger.info('Try to remove or rename your database file.')
    logger.info(f'(Current database file: {db_file})')
    ssp_exit(1)


def write_sqlite(config, sspec_output):
    try:
        database_file = config.database_file
    except KeyError:
        database_file = None
    if not database_file:
        return

    event = config.event
    evid = event.event_id
    runid = config.options.run_id

    # Current supported DB version
    DB_VERSION = 1

    # only check version if database_file exists
    check_version = os.path.isfile(database_file)

    # Open SQLite database
    try:
        conn = sqlite3.connect(database_file, timeout=60)
    except Exception as msg:
        logger.error(msg)
        logger.info(
            f'Please check whether "{database_file}" is a valid SQLite file.')
        ssp_exit(1)
    c = conn.cursor()

    if check_version:
        # Get current DB version
        db_version = c.execute('PRAGMA user_version').fetchone()[0]
        if db_version < DB_VERSION:
            logger.error(
                f'"{database_file}" has an old database version: '
                f'"{db_version}" Current supported version is "{DB_VERSION}".'
            )
            logger.info(
                'Remove or rename your old database file, '
                'so that a new one can be created.'
            )
            ssp_exit(1)
    else:
        # Set the DB version
        c.execute('PRAGMA user_version = {v:d}'.format(v=DB_VERSION))

    # Create Station table
    sql_create_stations_table = """CREATE TABLE IF NOT EXISTS Stations (
        stid TEXT,
        evid TEXT,
        runid TEXT,
        Mo REAL,
        Mo_err_minus REAL,
        Mo_err_plus REAL,
        Mo_is_outlier INT,
        Mw REAL,
        Mw_err_minus REAL,
        Mw_err_plus REAL,
        Mw_is_outlier INT,
        fc REAL,
        fc_err_minus REAL,
        fc_err_plus REAL,
        fc_is_outlier INT,
        t_star REAL,
        t_star_err_minus REAL,
        t_star_err_plus REAL,
        t_star_is_outlier INT,
        Qo REAL,
        Qo_err_minus REAL,
        Qo_err_plus REAL,
        Qo_is_outlier INT,
        bsd REAL,
        bsd_err_minus REAL,
        bsd_err_plus REAL,
        bsd_is_outlier INT,
        ra REAL,
        ra_err_minus REAL,
        ra_err_plus REAL,
        ra_is_outlier INT,
        Er REAL,
        Er_is_outlier INT,
        dist REAL,
        azimuth REAL
    );"""
    c.execute(sql_create_stations_table)
    # Write station source parameters to database
    nobs = 0
    stationpar = sspec_output.station_parameters
    sql_delete_from_stations =\
        'DELETE FROM Stations WHERE stid=? AND evid=? AND runid=?;'
    for statId in sorted(stationpar.keys()):
        nobs += 1
        par = stationpar[statId]
        # Remove existing line, if present
        t = (statId, evid, runid)
        try:
            c.execute(sql_delete_from_stations, t)
        except Exception as msg:
            _log_db_write_error(msg, database_file)
            ssp_exit(1)
        # Insert new line
        t = (
            statId, evid, runid,
            *par.Mo.value_uncertainty(),
            int(par.Mo.outlier),
            *par.Mw.value_uncertainty(),
            int(par.Mw.outlier),
            *par.fc.value_uncertainty(),
            int(par.fc.outlier),
            *par.t_star.value_uncertainty(),
            int(par.t_star.outlier),
            *par.Qo.value_uncertainty(),
            int(par.Qo.outlier),
            *par.bsd.value_uncertainty(),
            int(par.bsd.outlier),
            *par.radius.value_uncertainty(),
            int(par.radius.outlier),
            par.Er.value,
            int(par.Er.outlier),
            par.hypo_dist_in_km,
            par.azimuth
        )
        # Create a string like ?,?,?,?
        values = ','.join('?'*len(t))
        sql_insert_into_stations = f'INSERT INTO Stations VALUES({values});'
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
            evid TEXT,
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
    run_completed = f'{config.end_of_run} {config.end_of_run_tz}'
    ssp_version = get_versions()['version']
    # Remove event from Event table, if present
    t = (evid, runid)
    sql_delete_from_events = 'DELETE FROM Events WHERE evid=? AND runid=?;'
    try:
        c.execute(sql_delete_from_events, t)
    except Exception as msg:
        _log_db_write_error(msg, database_file)
        ssp_exit(1)
    ev_lon = event.hypocenter.longitude.value_in_deg
    ev_lat = event.hypocenter.latitude.value_in_deg
    ev_depth = event.hypocenter.depth.value_in_km
    ev_origin_time = event.hypocenter.origin_time
    t = (
        # Event info
        evid,
        runid,
        str(ev_origin_time),
        float(ev_lon),
        float(ev_lat),
        float(ev_depth),
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
    sql_insert_into_events = f'INSERT INTO Events VALUES({values});'
    try:
        c.execute(sql_insert_into_events, t)
    except Exception as msg:
        _log_db_write_error(msg, database_file)
        ssp_exit(1)
    # Commit changes and close database
    conn.commit()
    conn.close()
    logger.info(f'Output written to SQLite database: {database_file}')
