
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
        exit(1)
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
        exit(1)


def _version_1_to_2(cursor):
    """
    Update a version 1 database to version 2.

    :param cursor: SQLite cursor
    :type cursor: sqlite3.Cursor
    """
    # New in version 2:
    # Stations:
    #   - primary keys: stid, evid, runid
    #   - new keys: Mo_is_outlier, Mw_is_outlier, fc_is_outlier,
    #     t_star_is_outlier, Qo_is_outlier, bsd_is_outlier, ra_is_outlier,
    #     Er_is_outlier
    # Events:
    #   - primary keys: evid, runid
    sql_create_stations_table = (
        'CREATE TABLE IF NOT EXISTS StationsNew ('
        + '\n'.join(
            [f'{key} {value},' for key, value in STATIONS_TABLE.items()]
        )
        + 'PRIMARY KEY (' + ', '.join(STATIONS_PRIMARY_KEYS) + ')'
        + ');'
    )
    new_keys = [
        'Mo_is_outlier', 'Mw_is_outlier', 'fc_is_outlier', 't_star_is_outlier',
        'Qo_is_outlier', 'bsd_is_outlier', 'ra_is_outlier', 'Er_is_outlier'
    ]
    stations_keys = ', '.join([
        key for key in STATIONS_TABLE.keys() if key not in new_keys])
    sql_insert_stations = (
        f'INSERT INTO StationsNew ({stations_keys}) '
        f'SELECT {stations_keys} FROM Stations;'
    )
    sql_drop_stations = 'DROP TABLE Stations;'
    sql_rename_stations = 'ALTER TABLE StationsNew RENAME TO Stations;'
    sql_create_events_table = (
        'CREATE TABLE IF NOT EXISTS EventsNew ('
        + '\n'.join(
            [f'{key} {value},' for key, value in EVENTS_TABLE.items()]
        )
        + 'PRIMARY KEY (' + ', '.join(EVENTS_PRIMARY_KEYS) + ')'
        + ');'
    )
    events_keys = ', '.join(EVENTS_TABLE.keys())
    sql_insert_events = (
        f'INSERT INTO EventsNew ({events_keys}) '
        f'SELECT {events_keys} FROM Events;'
    )
    sql_drop_events = 'DROP TABLE Events;'
    sql_rename_events = 'ALTER TABLE EventsNew RENAME TO Events;'
    try:
        cursor.execute(sql_create_stations_table)
        cursor.execute(sql_insert_stations)
        cursor.execute(sql_drop_stations)
        cursor.execute(sql_rename_stations)
        cursor.execute(sql_create_events_table)
        cursor.execute(sql_insert_events)
        cursor.execute(sql_drop_events)
        cursor.execute(sql_rename_events)
        cursor.execute('PRAGMA user_version = 2;')
    except Exception as db_err:
        sys.stderr.write(f'{db_err}\n')
        exit(1)


def _overwrite_ok(db_file):
    """
    Check if db_file exists and ask for confirmation to overwrite it.

    :param db_file: SQLite database file
    :type db_file: str
    :return: True if overwrite is ok, False otherwise
    """
    if not os.path.exists(db_file):
        print(f'ERROR: {db_file} does not exist.')
        exit(1)
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
        exit(0)
    else:
        print(f'ERROR: {db_file} has an unsupported version {db_version}.')
        exit(1)
    conn.commit()
    conn.close()
