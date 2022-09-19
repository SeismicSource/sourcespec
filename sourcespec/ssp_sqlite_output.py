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

    # Open SQLite database
    conn = sqlite3.connect(database_file, timeout=60)
    c = conn.cursor()

    # Create Station table
    sql_create_stations_table = """CREATE TABLE IF NOT EXISTS Stations (
        stid, evid, runid,
        Mo, Mo_err_minus, Mo_err_plus,
        Mw, Mw_err_minus, Mw_err_plus,
        fc, fc_err_minus, fc_err_plus,
        t_star, t_star_err_minus, t_star_err_plus,
        Qo, Qo_err_minus, Qo_err_plus,
        bsd, bsd_err_minus, bsd_err_plus,
        ra, ra_err_minus, ra_err_plus,
        dist, azimuth, Er
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
        evid, runid, orig_time, lon, lat, depth, nobs,
        Mo, Mo_err_minus, Mo_err_plus,
        Mo_wavg, Mo_wavg_err_minus, Mo_wavg_err_plus,
        Mw, Mw_err, Mw_wavg, Mw_wavg_err,
        fc, fc_err_minus, fc_err_plus,
        fc_wavg, fc_wavg_err_minus, fc_wavg_err_plus,
        t_star, t_star_err,
        t_star_wavg, t_star_wavg_err,
        Qo, Qo_err,
        Qo_wavg, Qo_wavg_err,
        ra, ra_err_minus, ra_err_plus,
        ra_wavg, ra_wavg_err_minus, ra_wavg_err_plus,
        bsd, bsd_err_minus, bsd_err_plus,
        bsd_wavg, bsd_wavg_err_minus, bsd_wavg_err_plus,
        Er, Er_err_minus, Er_err_plus,
        Ml, Ml_err,
        run_completed, sourcespec_version,
        author_name, author_email,
        agency_full_name, agency_short_name, agency_url
    );"""
    c.execute(sql_create_events_table)
    means = sspec_output.mean_values()
    errors = sspec_output.mean_uncertainties()
    means_weight = sspec_output.weighted_mean_values()
    errors_weight = sspec_output.weighted_mean_uncertainties()
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
        evid, runid,
        str(hypo.origin_time),
        float(hypo.longitude), float(hypo.latitude), float(hypo.depth),
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
        means['radius'], *errors['radius'],
        means_weight['radius'], *errors_weight['radius'],
        means['bsd'], *errors['bsd'],
        means_weight['bsd'], *errors_weight['bsd'],
        means['Er'], *errors['Er'],
        means.get('Ml', None), errors.get('Ml', None),
        run_completed, ssp_version,
        config.author_name, config.author_email,
        config.agency_full_name, config.agency_short_name,
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
    logger.info('Output written to database: ' + database_file)
