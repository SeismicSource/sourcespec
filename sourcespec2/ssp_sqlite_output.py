# -*- coding: utf-8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
SQLite output for source_spec.

:copyright:
    2013-2026 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
import os.path
import logging
import sqlite3
from sourcespec.ssp_setup import ssp_exit
from sourcespec.ssp_db_definitions import (
    DB_VERSION,
    STATIONS_TABLE, STATIONS_PRIMARY_KEYS, EVENTS_TABLE, EVENTS_PRIMARY_KEYS)
from sourcespec._version import get_versions
logger = logging.getLogger(__name__.rsplit('.', maxsplit=1)[-1])


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
    if db_version == DB_VERSION:
        return
    if db_version > DB_VERSION:
        logger.error(
            f'"{db_file}" has a newer database version: '
            f'"{db_version}" Current supported version is "{DB_VERSION}".'
        )
        ssp_exit(1)
    logger.error(
        f'"{db_file}" has an old database version: '
        f'"{db_version}" Current supported version is "{DB_VERSION}".'
    )
    logger.info(
        'Use the following command to update your database '
        '(the current database will be backed up):\n\n'
        f'  source_spec --updatedb {db_file}\n'
    )
    ssp_exit(1)


def _set_db_version(cursor):
    """
    Set database version.

    :param cursor: SQLite cursor
    :type cursor: sqlite3.Cursor
    """
    cursor.execute(f'PRAGMA user_version = {DB_VERSION:d}')


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
        'Use the following command to update your database '
        '(the current database will be backed up):\n\n'
        f'  source_spec --updatedb {db_file}\n'
    )
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


def _insert_station_row(cursor, db_file, statId, par, evid, runid):
    """
    Insert a row in the Stations table.

    :param cursor: SQLite cursor
    :type cursor: sqlite3.Cursor
    :param db_file: SQLite database file
    :type db_file: str
    :param statId: Station ID
    :type statId: str
    :param par: Station parameters
    :type par: ssp_data_types.StationParameters
    :param evid: Event ID
    :type evid: str
    :param runid: Run ID
    :type runid: str
    """
    station_row = dict(STATIONS_TABLE)
    station_row |= {
        'stid': statId,
        'evid': evid,
        'runid': runid,
        'Mo': getattr(par.Mo, 'value', None),
        'Mo_err_minus': (par.Mo.compact_uncertainty()[0] if par.Mo else None),
        'Mo_err_plus': (par.Mo.compact_uncertainty()[1] if par.Mo else None),
        'Mo_is_outlier': int(getattr(par.Mo, 'outlier', True)),
        'Mw': getattr(par.Mw, 'value', None),
        'Mw_err_minus': (par.Mw.compact_uncertainty()[0] if par.Mw else None),
        'Mw_err_plus': (par.Mw.compact_uncertainty()[1] if par.Mw else None),
        'Mw_is_outlier': int(getattr(par.Mw, 'outlier', True)),
        'fc': getattr(par.fc, 'value', None),
        'fc_err_minus': (par.fc.compact_uncertainty()[0] if par.fc else None),
        'fc_err_plus': (par.fc.compact_uncertainty()[1] if par.fc else None),
        'fc_is_outlier': int(getattr(par.fc, 'outlier', True)),
        't_star': getattr(par.t_star, 'value', None),
        't_star_err_minus': (
            par.t_star.compact_uncertainty()[0] if par.t_star else None
        ),
        't_star_err_plus': (
            par.t_star.compact_uncertainty()[1] if par.t_star else None
        ),
        't_star_is_outlier': int(getattr(par.t_star, 'outlier', True)),
        'Qo': getattr(par.Qo, 'value', None),
        'Qo_err_minus': (par.Qo.compact_uncertainty()[0] if par.Qo else None),
        'Qo_err_plus': (par.Qo.compact_uncertainty()[1] if par.Qo else None),
        'Qo_is_outlier': int(getattr(par.Qo, 'outlier', True)),
        'ssd': getattr(par.ssd, 'value', None),
        'ssd_err_minus': (
            par.ssd.compact_uncertainty()[0] if par.ssd else None
        ),
        'ssd_err_plus': (
            par.ssd.compact_uncertainty()[1] if par.ssd else None
        ),
        'ssd_is_outlier': int(getattr(par.ssd, 'outlier', True)),
        'ra': getattr(par.radius, 'value', None),
        'ra_err_minus': (
            par.radius.compact_uncertainty()[0] if par.radius else None
        ),
        'ra_err_plus': (
            par.radius.compact_uncertainty()[1] if par.radius else None
        ),
        'ra_is_outlier': int(getattr(par.radius, 'outlier', True)),
        'Er': getattr(par.Er, 'value', None),
        'Er_err_minus': (par.Er.compact_uncertainty()[0] if par.Er else None),
        'Er_err_plus': (par.Er.compact_uncertainty()[1] if par.Er else None),
        'Er_is_outlier': int(getattr(par.Er, 'outlier', True)),
        'sigma_a': getattr(par.sigma_a, 'value', None),
        'sigma_a_err_minus': (
            par.sigma_a.compact_uncertainty()[0] if par.sigma_a else None
        ),
        'sigma_a_err_plus': (
            par.sigma_a.compact_uncertainty()[1] if par.sigma_a else None
        ),
        'sigma_a_is_outlier': int(getattr(par.sigma_a, 'outlier', True)),
        'lon': par.longitude,
        'lat': par.latitude,
        'instr_type': par.instrument_type,
        'hypo_dist': par.hypo_dist_in_km,
        'epi_dist': par.epi_dist_in_km,
        'azimuth': par.azimuth,
        'spectral_snratio_mean': par.spectral_snratio_mean,
        'spectral_snratio_max': par.spectral_snratio_max,
        'rmsn': getattr(par, 'rmsn', None),
        'quality_of_fit': getattr(par, 'quality_of_fit', None),
        'ignored': par.ignored,
        'ignored_reason': getattr(par, 'ignored_reason', None),
    }
    columns = list(station_row.keys())
    row = tuple(station_row[col] for col in columns)
    sql_insert_into_stations = (
        f'INSERT OR REPLACE INTO Stations '
        f'({",".join(columns)}) VALUES({",".join("?" for _ in columns)});'
    )
    try:
        cursor.execute(sql_insert_into_stations, row)
    except Exception as msg:
        _log_db_write_error(msg, db_file)
        ssp_exit(1)


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
    for statId in sorted(stationpar.keys()):
        par = stationpar[statId]
        _insert_station_row(cursor, db_file, statId, par, evid, runid)


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


def _write_events_table(cursor, db_file, sspec_output, config):
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
    """
    event = config.event
    evid = event.event_id
    runid = config.options.run_id
    wave_type = config.wave_type
    n_input_stations = sspec_output.quality_info.n_input_stations
    n_input_spectra = sspec_output.quality_info.n_input_spectra
    n_spectra_inverted = sspec_output.quality_info.n_spectra_inverted
    azimuthal_gap_primary = sspec_output.quality_info.azimuthal_gap_primary
    azimuthal_gap_secondary = sspec_output.quality_info.azimuthal_gap_secondary
    rmsn_mean = sspec_output.quality_info.rmsn_mean
    quality_of_fit_mean = sspec_output.quality_info.quality_of_fit_mean
    spectral_dispersion_rmsn = \
        sspec_output.quality_info.spectral_dispersion_rmsn
    spectral_dispersion_score = \
        sspec_output.quality_info.spectral_dispersion_score
    means = sspec_output.mean_values()
    mean_errors = sspec_output.mean_uncertainties()
    mean_nobs = sspec_output.mean_nobs()
    wmeans = sspec_output.weighted_mean_values()
    wmean_errors = sspec_output.weighted_mean_uncertainties()
    wmean_nobs = sspec_output.weighted_mean_nobs()
    percentiles = sspec_output.percentiles_values()
    percentile_errors = sspec_output.percentiles_uncertainties()
    percentile_nobs = sspec_output.percentiles_nobs()
    run_completed = f'{config.end_of_run} {config.end_of_run_tz}'
    ssp_version = get_versions()['version']
    ev_lon = event.hypocenter.longitude.value_in_deg
    ev_lat = event.hypocenter.latitude.value_in_deg
    ev_depth = event.hypocenter.depth.value_in_km
    ev_origin_time = event.hypocenter.origin_time
    ev_vp = event.hypocenter.vp
    ev_vs = event.hypocenter.vs
    ev_rho = event.hypocenter.rho
    kp = config.kp
    ks = config.ks
    event_row = dict(EVENTS_TABLE)
    # Update event_row dictionary with event and run information
    event_row |= {
        # Event info
        'evid': evid,
        'runid': runid,
        'orig_time': str(ev_origin_time),
        'lon': float(ev_lon),
        'lat': float(ev_lat),
        'depth': float(ev_depth),
        'vp': float(ev_vp),
        'vs': float(ev_vs),
        'rho': float(ev_rho),
        'kp': kp,
        'ks': ks,
        # Inversion info
        'wave_type': wave_type,
        # Quality info
        'n_input_stations': n_input_stations,
        'n_input_spectra': n_input_spectra,
        'n_spectra_inverted': n_spectra_inverted,
        'azimuthal_gap_primary': azimuthal_gap_primary,
        'azimuthal_gap_secondary': azimuthal_gap_secondary,
        'rmsn_mean': rmsn_mean,
        'quality_of_fit_mean': quality_of_fit_mean,
        'spectral_dispersion_rmsn': spectral_dispersion_rmsn,
        'spectral_dispersion_score': spectral_dispersion_score,
        # Statistical info
        'nsigma': config.n_sigma,
        'mid_pct': config.mid_percentage,
        'lower_pct': config.lower_percentage,
        'upper_pct': config.upper_percentage,
        # Seismic moment
        'Mo_mean': means['Mo'],
        'Mo_mean_err_minus': mean_errors['Mo'][0],
        'Mo_mean_err_plus': mean_errors['Mo'][1],
        'Mo_mean_nobs': mean_nobs['Mo'],
        'Mo_wmean': wmeans['Mo'],
        'Mo_wmean_err_minus': wmean_errors['Mo'][0],
        'Mo_wmean_err_plus': wmean_errors['Mo'][1],
        'Mo_wmean_nobs': wmean_nobs['Mo'],
        'Mo_pctl': percentiles['Mo'],
        'Mo_pctl_err_minus': percentile_errors['Mo'][0],
        'Mo_pctl_err_plus': percentile_errors['Mo'][1],
        'Mo_pctl_nobs': percentile_nobs['Mo'],
        # Moment magnitude
        'Mw_mean': means['Mw'],
        'Mw_mean_err_minus': mean_errors['Mw'][0],
        'Mw_mean_err_plus': mean_errors['Mw'][1],
        'Mw_mean_nobs': mean_nobs['Mw'],
        'Mw_wmean': wmeans['Mw'],
        'Mw_wmean_err_minus': wmean_errors['Mw'][0],
        'Mw_wmean_err_plus': wmean_errors['Mw'][1],
        'Mw_wmean_nobs': wmean_nobs['Mw'],
        'Mw_pctl': percentiles['Mw'],
        'Mw_pctl_err_minus': percentile_errors['Mw'][0],
        'Mw_pctl_err_plus': percentile_errors['Mw'][1],
        'Mw_pctl_nobs': percentile_nobs['Mw'],
        # Corner frequency
        'fc_mean': means['fc'],
        'fc_mean_err_minus': mean_errors['fc'][0],
        'fc_mean_err_plus': mean_errors['fc'][1],
        'fc_mean_nobs': mean_nobs['fc'],
        'fc_wmean': wmeans['fc'],
        'fc_wmean_err_minus': wmean_errors['fc'][0],
        'fc_wmean_err_plus': wmean_errors['fc'][1],
        'fc_wmean_nobs': wmean_nobs['fc'],
        'fc_pctl': percentiles['fc'],
        'fc_pctl_err_minus': percentile_errors['fc'][0],
        'fc_pctl_err_plus': percentile_errors['fc'][1],
        'fc_pctl_nobs': percentile_nobs['fc'],
        # t-star
        't_star_mean': means['t_star'],
        't_star_mean_err_minus': mean_errors['t_star'][0],
        't_star_mean_err_plus': mean_errors['t_star'][1],
        't_star_mean_nobs': mean_nobs['t_star'],
        't_star_wmean': wmeans['t_star'],
        't_star_wmean_err_minus': wmean_errors['t_star'][0],
        't_star_wmean_err_plus': wmean_errors['t_star'][1],
        't_star_wmean_nobs': wmean_nobs['t_star'],
        't_star_pctl': percentiles['t_star'],
        't_star_pctl_err_minus': percentile_errors['t_star'][0],
        't_star_pctl_err_plus': percentile_errors['t_star'][1],
        't_star_pctl_nobs': percentile_nobs['t_star'],
        # Qo
        'Qo_mean': means['Qo'],
        'Qo_mean_err_minus': mean_errors['Qo'][0],
        'Qo_mean_err_plus': mean_errors['Qo'][1],
        'Qo_mean_nobs': mean_nobs['Qo'],
        'Qo_wmean': wmeans['Qo'],
        'Qo_wmean_err_minus': wmean_errors['Qo'][0],
        'Qo_wmean_err_plus': wmean_errors['Qo'][1],
        'Qo_wmean_nobs': wmean_nobs['Qo'],
        'Qo_pctl': percentiles['Qo'],
        'Qo_pctl_err_minus': percentile_errors['Qo'][0],
        'Qo_pctl_err_plus': percentile_errors['Qo'][1],
        'Qo_pctl_nobs': percentile_nobs['Qo'],
        # Source radius
        'ra_mean': means['radius'],
        'ra_mean_err_minus': mean_errors['radius'][0],
        'ra_mean_err_plus': mean_errors['radius'][1],
        'ra_mean_nobs': mean_nobs['radius'],
        'ra_wmean': wmeans['radius'],
        'ra_wmean_err_minus': wmean_errors['radius'][0],
        'ra_wmean_err_plus': wmean_errors['radius'][1],
        'ra_wmean_nobs': wmean_nobs['radius'],
        'ra_pctl': percentiles['radius'],
        'ra_pctl_err_minus': percentile_errors['radius'][0],
        'ra_pctl_err_plus': percentile_errors['radius'][1],
        'ra_pctl_nobs': percentile_nobs['radius'],
        # Static stress drop
        'ssd_mean': means['ssd'],
        'ssd_mean_err_minus': mean_errors['ssd'][0],
        'ssd_mean_err_plus': mean_errors['ssd'][1],
        'ssd_mean_nobs': mean_nobs['ssd'],
        'ssd_wmean': wmeans['ssd'],
        'ssd_wmean_err_minus': wmean_errors['ssd'][0],
        'ssd_wmean_err_plus': wmean_errors['ssd'][1],
        'ssd_wmean_nobs': wmean_nobs['ssd'],
        'ssd_pctl': percentiles['ssd'],
        'ssd_pctl_err_minus': percentile_errors['ssd'][0],
        'ssd_pctl_err_plus': percentile_errors['ssd'][1],
        'ssd_pctl_nobs': percentile_nobs['ssd'],
        # Radiated energy
        'Er_mean': means['Er'],
        'Er_mean_err_minus': mean_errors['Er'][0],
        'Er_mean_err_plus': mean_errors['Er'][1],
        'Er_mean_nobs': mean_nobs['Er'],
        'Er_wmean': wmeans['Er'],
        'Er_wmean_err_minus': wmean_errors['Er'][0],
        'Er_wmean_err_plus': wmean_errors['Er'][1],
        'Er_wmean_nobs': wmean_nobs['Er'],
        'Er_pctl': percentiles['Er'],
        'Er_pctl_err_minus': percentile_errors['Er'][0],
        'Er_pctl_err_plus': percentile_errors['Er'][1],
        'Er_pctl_nobs': percentile_nobs['Er'],
        # Apparent stress
        'sigma_a_mean': means['sigma_a'],
        'sigma_a_mean_err_minus': mean_errors['sigma_a'][0],
        'sigma_a_mean_err_plus': mean_errors['sigma_a'][1],
        'sigma_a_mean_nobs': mean_nobs['sigma_a'],
        'sigma_a_wmean': wmeans['sigma_a'],
        'sigma_a_wmean_err_minus': wmean_errors['sigma_a'][0],
        'sigma_a_wmean_err_plus': wmean_errors['sigma_a'][1],
        'sigma_a_wmean_nobs': wmean_nobs['sigma_a'],
        'sigma_a_pctl': percentiles['sigma_a'],
        'sigma_a_pctl_err_minus': percentile_errors['sigma_a'][0],
        'sigma_a_pctl_err_plus': percentile_errors['sigma_a'][1],
        'sigma_a_pctl_nobs': percentile_nobs['sigma_a'],
        # Local magnitude
        'Ml_mean': means.get('Ml', None),
        'Ml_mean_err_minus': mean_errors.get('Ml', (None, None))[0],
        'Ml_mean_err_plus': mean_errors.get('Ml', (None, None))[1],
        'Ml_mean_nobs': mean_nobs.get('Ml', None),
        'Ml_wmean': wmeans.get('Ml', None),
        'Ml_wmean_err_minus': wmean_errors.get('Ml', (None, None))[0],
        'Ml_wmean_err_plus': wmean_errors.get('Ml', (None, None))[1],
        'Ml_wmean_nobs': wmean_nobs.get('Ml', None),
        'Ml_pctl': percentiles.get('Ml', None),
        'Ml_pctl_err_minus': percentile_errors.get('Ml', (None, None))[0],
        'Ml_pctl_err_plus': percentile_errors.get('Ml', (None, None))[1],
        'Ml_pctl_nobs': percentile_nobs.get('Ml', None),
        # Run info
        'run_completed': run_completed,
        'sourcespec_version': ssp_version,
        'author_name': config.author_name,
        'author_email': config.author_email,
        'agency_full_name': config.agency_full_name,
        'agency_short_name': config.agency_short_name,
        'agency_url': config.agency_url
    }
    # Explicitly specify column names to avoid order issues
    columns = list(event_row.keys())
    row = tuple(event_row[col] for col in columns)
    sql_insert_into_events = (
        f'INSERT OR REPLACE INTO Events '
        f'({",".join(columns)}) VALUES({",".join("?" for _ in columns)});'
    )
    try:
        cursor.execute(sql_insert_into_events, row)
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
    _write_stations_table(cursor, db_file, sspec_output, config)
    # Commit changes
    conn.commit()
    # Create Events table
    _create_events_table(cursor, db_file)
    # Write event source parameters to database
    _write_events_table(cursor, db_file, sspec_output, config)
    # Commit changes and close database
    conn.commit()
    conn.close()
    logger.info(f'Output written to SQLite database: {db_file}')
