
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Update an existing SourceSpec database from a previous version.

:copyright:
    2013-2023 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
import os
import sys
import shutil
import sqlite3
from sourcespec.ssp_db_definitions import (
    DB_VERSION,
    STATIONS_TABLE, STATIONS_PRIMARY_KEYS, EVENTS_TABLE, EVENTS_PRIMARY_KEYS)


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
        sys.stderr.write(f'{msg}\n')
        sys.stderr.write(
            f'Please check whether "{db_file}" is a valid SQLite file.\n')
        sys.exit(1)
    return conn, conn.cursor()


def _get_db_version(cursor, db_file):
    """
    Get database version.

    :param cursor: SQLite cursor
    :type cursor: sqlite3.Cursor
    :param db_file: SQLite database file
    :type db_file: str
    """
    try:
        return cursor.execute('PRAGMA user_version').fetchone()[0]
    except Exception as msg:
        sys.stderr.write(f'{msg}\n')
        sys.stderr.write(
            f'Please check whether "{db_file}" is a valid SQLite file.\n')
        sys.exit(1)


def _version_1_to_2(cursor):
    """
    Update a version 1 database to version 2.

    :param cursor: SQLite cursor
    :type cursor: sqlite3.Cursor
    """
    # Stations table:
    # New in version 2:
    #   - primary keys: stid, evid, runid
    #   - renamed keys:
    #     bsd -> ssd,
    #     bsd_err_minus -> ssd_err_minus,
    #     bsd_err_plus -> ssd_err_plus
    #   - new keys:
    #     Mo_is_outlier, Mw_is_outlier, fc_is_outlier,
    #     t_star_is_outlier, Qo_is_outlier, ssd_is_outlier, ra_is_outlier,
    #     Er_err_minus, Er_err_plus, Er_is_outlier
    #     sigma_a, sigma_a_err_minus, sigma_a_err_plus, sigma_a_is_outlier
    renamed_station_keys = {
        'bsd': 'ssd',
        'bsd_err_minus': 'ssd_err_minus',
        'bsd_err_plus': 'ssd_err_plus'
    }
    list_sql_rename_station_keys = [
        f'ALTER TABLE Stations RENAME COLUMN {key} TO {value};'
        for key, value in renamed_station_keys.items()]
    sql_create_new_stations_table = (
        'CREATE TABLE IF NOT EXISTS StationsNew ('
        + '\n'.join(
            [f'{key} {value},' for key, value in STATIONS_TABLE.items()]
        )
        + 'PRIMARY KEY (' + ', '.join(STATIONS_PRIMARY_KEYS) + ')'
        + ');'
    )
    new_station_keys = [
        'Mo_is_outlier', 'Mw_is_outlier', 'fc_is_outlier', 't_star_is_outlier',
        'Qo_is_outlier', 'ssd_is_outlier', 'ra_is_outlier',
        'Er_err_minus', 'Er_err_plus', 'Er_is_outlier',
        'sigma_a', 'sigma_a_err_minus', 'sigma_a_err_plus',
        'sigma_a_is_outlier'
    ]
    station_keys = ', '.join([
        key for key in STATIONS_TABLE if key not in new_station_keys])
    sql_insert_new_station_keys = (
        f'INSERT INTO StationsNew ({station_keys}) '
        f'SELECT {station_keys} FROM Stations;'
    )
    sql_drop_old_stations_table = 'DROP TABLE Stations;'
    sql_rename_new_stations_table =\
        'ALTER TABLE StationsNew RENAME TO Stations;'

    # Events table:
    # New in version 2:
    #   - primary keys: evid, runid
    #   - renamed keys:
    #     bsd_mean -> ssd_mean,
    #     bsd_mean_err_minus -> ssd_mean_err_minus,
    #     bsd_mean_err_plus -> ssd_mean_err_plus,
    #     bsd_wmean -> ssd_wmean,
    #     bsd_wmean_err_minus -> ssd_wmean_err_minus,
    #     bsd_wmean_err_plus -> ssd_wmean_err_plus,
    #     bsd_pctl -> ssd_pctl,
    #     bsd_pctl_err_minus -> ssd_pctl_err_minus,
    #     bsd_pctl_err_plus -> ssd_pctl_err_plus
    #   - new keys:
    #     vp, vs, rho, wave_type,
    #     Er_wmean, Er_wmean_err_minus, Er_wmean_err_plus,
    #     Ml_wmean, Ml_wmean_err_minus, Ml_wmean_err_plus,
    #     sigma_a_mean, sigma_a_mean_err_minus, sigma_a_mean_err_plus,
    #     sigma_a_wmean, sigma_a_wmean_err_minus, sigma_a_wmean_err_plus,
    #     sigma_a_pctl, sigma_a_pctl_err_minus, sigma_a_pctl_err_plus
    renamed_event_keys = {
        'bsd_mean': 'ssd_mean',
        'bsd_mean_err_minus': 'ssd_mean_err_minus',
        'bsd_mean_err_plus': 'ssd_mean_err_plus',
        'bsd_wmean': 'ssd_wmean',
        'bsd_wmean_err_minus': 'ssd_wmean_err_minus',
        'bsd_wmean_err_plus': 'ssd_wmean_err_plus',
        'bsd_pctl': 'ssd_pctl',
        'bsd_pctl_err_minus': 'ssd_pctl_err_minus',
        'bsd_pctl_err_plus': 'ssd_pctl_err_plus'
    }
    list_sql_rename_event_keys = [
        f'ALTER TABLE Events RENAME COLUMN {key} TO {value};'
        for key, value in renamed_event_keys.items()]
    sql_create_new_events_table = (
        'CREATE TABLE IF NOT EXISTS EventsNew ('
        + '\n'.join(
            [f'{key} {value},' for key, value in EVENTS_TABLE.items()]
        )
        + 'PRIMARY KEY (' + ', '.join(EVENTS_PRIMARY_KEYS) + ')'
        + ');'
    )
    new_event_keys = [
        'vp', 'vs', 'rho', 'wave_type',
        'Er_wmean', 'Er_wmean_err_minus', 'Er_wmean_err_plus',
        'Ml_wmean', 'Ml_wmean_err_minus', 'Ml_wmean_err_plus',
        'sigma_a_mean', 'sigma_a_mean_err_minus', 'sigma_a_mean_err_plus',
        'sigma_a_wmean', 'sigma_a_wmean_err_minus', 'sigma_a_wmean_err_plus',
        'sigma_a_pctl', 'sigma_a_pctl_err_minus', 'sigma_a_pctl_err_plus'
    ]
    event_keys = ', '.join([
        key for key in EVENTS_TABLE if key not in new_event_keys])
    sql_insert_new_event_keys = (
        f'INSERT INTO EventsNew ({event_keys}) '
        f'SELECT {event_keys} FROM Events;'
    )
    sql_drop_old_events_table = 'DROP TABLE Events;'
    sql_rename_new_events_table = 'ALTER TABLE EventsNew RENAME TO Events;'

    # execute SQL statements
    try:
        # stations table
        for statement in list_sql_rename_station_keys:
            cursor.execute(statement)
        cursor.execute(sql_create_new_stations_table)
        cursor.execute(sql_insert_new_station_keys)
        cursor.execute(sql_drop_old_stations_table)
        cursor.execute(sql_rename_new_stations_table)
        # events table
        for statement in list_sql_rename_event_keys:
            cursor.execute(statement)
        cursor.execute(sql_create_new_events_table)
        cursor.execute(sql_insert_new_event_keys)
        cursor.execute(sql_drop_old_events_table)
        cursor.execute(sql_rename_new_events_table)
        cursor.execute('PRAGMA user_version = 2;')
    except Exception as db_err:
        sys.stderr.write(f'{db_err}\n')
        sys.exit(1)


def _overwrite_ok(db_file):
    """
    Check if db_file exists and ask for confirmation to overwrite it.

    :param db_file: SQLite database file
    :type db_file: str
    :return: True if overwrite is ok, False otherwise
    """
    if not os.path.exists(db_file):
        print(f'ERROR: {db_file} does not exist.')
        sys.exit(1)
    answer = input(
        f'Overwrite {db_file}?\n'
        f'(current file will be saved to {db_file}.bak) [y/N] ')
    return answer.lower() == 'y'


def update_db_file(db_file):
    """
    Update an existing SourceSpec database from a previous version.

    :param db_file: SQLite database file
    :type db_file: str
    """
    if not _overwrite_ok(db_file):
        return
    print(f'Updating {db_file}...')
    conn, cursor = _open_sqlite_db(db_file)
    db_version = _get_db_version(cursor, db_file)
    if db_version == 1:
        # create a backup copy
        shutil.copy2(db_file, f'{db_file}.bak')
        _version_1_to_2(cursor)
        print(
            f'{db_file} updated from version {db_version} '
            f'to version {DB_VERSION}.')
    elif db_version == DB_VERSION:
        print(f'{db_file} is already up-to-date.')
        sys.exit(0)
    else:
        print(f'ERROR: {db_file} has an unsupported version {db_version}.')
        sys.exit(1)
    conn.commit()
    conn.close()
