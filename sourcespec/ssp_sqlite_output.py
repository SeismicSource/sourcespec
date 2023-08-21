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
import numpy as np
from sourcespec.ssp_setup import ssp_exit
from sourcespec.ssp_db_definitions import (
    DB_VERSION,
    STATIONS_TABLE, STATIONS_PRIMARY_KEYS, EVENTS_TABLE, EVENTS_PRIMARY_KEYS)
from sourcespec._version import get_versions
logger = logging.getLogger(__name__.split('.')[-1])


def _db_file_exists(db_file):
    """
    Check if SQLite database file exists.

    :param db_file: SQLite database file
    :type db_file: str
    :return: True if file exists, False otherwise
    :rtype: bool
    """
    return os.path.isfile(db_file)


def _open_sqlite_db(db_file):
    """
    Open SQLite database.

    :param db_file: SQLite database file
    :type db_file: str
    :return: SQLite connection and cursor
    :rtype: tuple
    """
    try:
        conn = sqlite3.connect(db_file, timeout=60)
    except Exception as msg:
        logger.error(msg)
        logger.info(
            f'Please check whether "{db_file}" is a valid SQLite file.')
        ssp_exit(1)
    return conn, conn.cursor()


def _check_db_version(cursor, db_file):
    """
    Check database version.

    :param cursor: SQLite cursor
    :type cursor: sqlite3.Cursor
    :param db_file: SQLite database file
    :type db_file: str
    """
    db_version = cursor.execute('PRAGMA user_version').fetchone()[0]
    if db_version < DB_VERSION:
        logger.error(
            f'"{db_file}" has an old database version: '
            f'"{db_version}" Current supported version is "{DB_VERSION}".'
        )
        logger.info(
            'Use "source_spec --updatedb" to update your database. '
            'The current database will be backed up.'
        )
        ssp_exit(1)


def _set_db_version(cursor):
    """
    Set database version.

    :param cursor: SQLite cursor
    :type cursor: sqlite3.Cursor
    """
    cursor.execute('PRAGMA user_version = {v:d}'.format(v=DB_VERSION))


def _log_db_write_error(db_err, db_file):
    """
    Log database write error.

    :param db_err: database error
    :type db_err: Exception
    :param db_file: SQLite database file
    :type db_file: str
    """
    logger.error(f'Unable to insert values: {db_err}')
    logger.info('Maybe your sqlite database has an old format.')
    logger.info(
        'Use "source_spec --updatedb" to update your database. '
        'The current database will be backed up.'
    )
    logger.info(f'(Current database file: {db_file})')
    ssp_exit(1)


def _create_stations_table(cursor, db_file):
    """
    Create Stations table.

    :param cursor: SQLite cursor
    :type cursor: sqlite3.Cursor
    """
    sql_create_stations_table = (
        'CREATE TABLE IF NOT EXISTS Stations ('
        + '\n'.join(
            [f'{key} {value},' for key, value in STATIONS_TABLE.items()]
        )
        + 'PRIMARY KEY (' + ', '.join(STATIONS_PRIMARY_KEYS) + ')'
        + ');'
    )
    try:
        cursor.execute(sql_create_stations_table)
    except Exception as db_err:
        _log_db_write_error(db_err, db_file)


def _write_stations_table(cursor, db_file, sspec_output, config):
    """
    Write station source parameters to database.

    :param cursor: SQLite cursor
    :type cursor: sqlite3.Cursor
    :param db_file: SQLite database file
    :type db_file: str
    :param sspec_output: sspec output object
    :type sspec_output: ssp_data_types.SourceSpecOutput
    :param config: sspec configuration object
    :type config: config.Config
    """
    event = config.event
    evid = event.event_id
    runid = config.options.run_id
    stationpar = sspec_output.station_parameters
    nobs = 0
    for statId in sorted(stationpar.keys()):
        nobs += 1
        par = stationpar[statId]
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
            # Er_err_minus and Er_err_plus are not defined
            np.nan, np.nan,
            int(par.Er.outlier),
            par.hypo_dist_in_km,
            par.azimuth
        )
        # Create a string like ?,?,?,?
        values = ','.join('?' * len(t))
        sql_insert_into_stations =\
            f'INSERT OR REPLACE INTO Stations VALUES({values});'
        try:
            cursor.execute(sql_insert_into_stations, t)
        except Exception as msg:
            _log_db_write_error(msg, db_file)
            ssp_exit(1)
    return nobs


def _create_events_table(cursor, db_file):
    """
    Create Events table.

    :param cursor: SQLite cursor
    :type cursor: sqlite3.Cursor
    """
    sql_create_events_table = (
        'CREATE TABLE IF NOT EXISTS Events ('
        + '\n'.join(
            [f'{key} {value},' for key, value in EVENTS_TABLE.items()]
        )
        + 'PRIMARY KEY (' + ', '.join(EVENTS_PRIMARY_KEYS) + ')'
        + ');'
    )
    try:
        cursor.execute(sql_create_events_table)
    except Exception as db_err:
        _log_db_write_error(db_err, db_file)


def _write_events_table(cursor, db_file, sspec_output, config, nobs):
    """
    Write Events table.

    :param cursor: SQLite cursor
    :type cursor: sqlite3.Cursor
    :param db_file: SQLite database file
    :type db_file: str
    :param sspec_output: SSP output object
    :type sspec_output: ssp_data_types.SourceSpecOutput
    :param config: SSP configuration object
    :type config: config.Config
    :param nobs: Number of observations
    :type nobs: int
    """
    event = config.event
    evid = event.event_id
    runid = config.options.run_id
    means = sspec_output.mean_values()
    mean_errors = sspec_output.mean_uncertainties()
    wmeans = sspec_output.weighted_mean_values()
    wmean_errors = sspec_output.weighted_mean_uncertainties()
    percentiles = sspec_output.percentiles_values()
    percentile_errors = sspec_output.percentiles_uncertainties()
    run_completed = f'{config.end_of_run} {config.end_of_run_tz}'
    ssp_version = get_versions()['version']
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
        wmeans['Er'],
        *wmean_errors['Er'],
        percentiles['Er'],
        *percentile_errors['Er'],
        # Local magnitude
        means.get('Ml', None),
        *mean_errors.get('Ml', (None, None)),
        wmeans.get('Ml', None),
        *wmean_errors.get('Ml', (None, None)),
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
    values = ','.join('?' * len(t))
    sql_insert_into_events = f'INSERT OR REPLACE INTO Events VALUES({values});'
    try:
        cursor.execute(sql_insert_into_events, t)
    except Exception as msg:
        _log_db_write_error(msg, db_file)
        ssp_exit(1)


def write_sqlite(config, sspec_output):
    """
    Write SSP output to SQLite database.

    :param config: SSP configuration object
    :type config: config.Config
    :param sspec_output: SSP output object
    :type sspec_output: ssp_data_types.SourceSpecOutput
    """
    db_file = config.get('database_file', None)
    if not db_file:
        return

    db_file_exists = _db_file_exists(db_file)
    conn, cursor = _open_sqlite_db(db_file)
    if db_file_exists:
        _check_db_version(cursor, db_file)
    else:
        _set_db_version(cursor)

    # Create Stations table
    _create_stations_table(cursor, db_file)
    # Write station source parameters to database
    nobs = _write_stations_table(cursor, db_file, sspec_output, config)
    # Commit changes
    conn.commit()
    # Create Events table
    _create_events_table(cursor, db_file)
    # Write event source parameters to database
    _write_events_table(cursor, db_file, sspec_output, config, nobs)
    # Commit changes and close database
    conn.commit()
    conn.close()
    logger.info(f'Output written to SQLite database: {db_file}')
