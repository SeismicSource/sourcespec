# -*- coding: utf-8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Parse event metadata and picks in hypo71 format.

:copyright:
    2012-2024 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
import os
import logging
from datetime import datetime
from obspy import UTCDateTime
from ...setup import config, ssp_exit
from ...ssp_event import SSPEvent
from ...ssp_pick import SSPPick
logger = logging.getLogger(__name__.rsplit('.', maxsplit=1)[-1])


def _is_hypo71_picks(pick_file):
    """
    Check if a file is a hypo71 phase file.

    :param pick_file: path to hypo71 phase file
    :type pick_file: str

    :raises TypeError: if the file is not a hypo71 phase file
    """
    with open(pick_file, encoding='ascii') as fp:
        for line in fp:
            # remove newline
            line = line.replace('\n', '')
            # skip separator and empty lines
            stripped_line = line.strip()
            if stripped_line in ['10', '']:
                continue
            # Check if it is a pick line
            # 6th character should be alpha (phase name: P or S)
            # other character should be digits (date/time)
            if not (line[5].isalpha() and
                    line[9].isdigit() and
                    line[20].isdigit()):
                raise TypeError(f'{pick_file}: Not a hypo71 phase file')


def _correct_station_name(station):
    """
    Correct station name, based on a traceid map.

    :param station: station name
    :type station: str

    :return: corrected station name
    :rtype: str
    """
    if config.TRACEID_MAP is None:
        return station
    # get all the keys containing station name in it
    keys = [key for key in config.TRACEID_MAP if station == key.split('.')[1]]
    # then take just the first one
    try:
        key = keys[0]
    except IndexError:
        return station
    traceid = config.TRACEID_MAP[key]
    return traceid.split('.')[1]


def parse_hypo71_picks():
    """
    Parse hypo71 picks file

    :return: list of SSPPick objects
    :rtype: list
    """
    picks = []
    pick_file = config.options.pick_file
    if pick_file is None:
        return picks
    try:
        _is_hypo71_picks(pick_file)
    except Exception as err:
        logger.error(err)
        ssp_exit(1)
    with open(pick_file, encoding='ascii') as fp:
        for line in fp:
            # remove newline
            line = line.replace('\n', '')
            # skip separator and empty lines
            stripped_line = line.strip()
            if stripped_line in ['10', '']:
                continue
            # Check if it is a pick line
            # 6th character should be alpha (phase name: P or S)
            # other character should be digits (date/time)
            if not (line[5].isalpha() and
                    line[9].isdigit() and
                    line[20].isdigit()):
                continue
            pick = SSPPick()
            pick.station = line[:4].strip()
            pick.station = _correct_station_name(pick.station)
            pick.flag = line[4:5]
            pick.phase = line[5:6]
            pick.polarity = line[6:7]
            try:
                pick.quality = int(line[7:8])
            except ValueError:
                # If we cannot read pick quality,
                # we give the pick the lowest quality
                pick.quality = 4
            timestr = line[9:24]
            dt = datetime.strptime(timestr, '%y%m%d%H%M%S.%f')
            pick.time = UTCDateTime(dt)
            picks.append(pick)
            try:
                stime = line[31:36]
            except Exception:
                continue
            if stime.strip() == '':
                continue
            pick2 = SSPPick()
            pick2.station = pick.station
            pick2.flag = line[36:37]
            pick2.phase = line[37:38]
            pick2.polarity = line[38:39]
            try:
                pick2.quality = int(line[39:40])
            except ValueError:
                # If we cannot read pick quality,
                # we give the pick the lowest quality
                pick2.quality = 4
            # pick2.time has the same date, hour and minutes
            # than pick.time
            # We therefore make a copy of pick.time,
            # and set seconds and microseconds to 0
            pick2.time = pick.time.replace(second=0, microsecond=0)
            # finally we add stime
            pick2.time += float(stime)
            picks.append(pick2)
    return picks


def parse_hypo71_hypocenter(hypo_file, _):
    """
    Parse a hypo71 hypocenter file.

    :param hypo_file: path to hypo71 hypocenter file
    :type hypo_file: str
    :param _: unused (for consistency with other parsers)
    :type _: None

    :return: a tuple of (SSPEvent, picks)
    :rtype: tuple

    .. note::
        The returned picks list is empty, for consistency with other parsers.
    """
    with open(hypo_file, encoding='ascii') as fp:
        line = fp.readline()
        # Skip the first line if it contain characters in the first 10 digits:
        if any(c.isalpha() for c in line[:10]):
            line = fp.readline()
    timestr = line[:17]
    # There are two possible formats for the timestring. We try both of them
    try:
        dt = datetime.strptime(timestr, '%y%m%d %H %M%S.%f')
    except Exception:
        try:
            dt = datetime.strptime(timestr, '%y%m%d %H%M %S.%f')
        except Exception as e:
            raise ValueError('Cannot read origin time on first line.') from e
    ssp_event = SSPEvent()
    hypo = ssp_event.hypocenter
    hypo.origin_time = UTCDateTime(dt)
    lat = float(line[17:20])
    lat_deg = float(line[21:26])
    hypo.latitude.value_in_deg = lat + lat_deg / 60
    lon = float(line[26:30])
    lon_deg = float(line[31:36])
    hypo.longitude.value_in_deg = lon + lon_deg / 60
    hypo.depth.value = float(line[36:42])
    hypo.depth.units = 'km'
    evid = os.path.basename(hypo_file)
    evid = evid.replace('.phs', '').replace('.h', '').replace('.hyp', '')
    ssp_event.event_id = evid
    # empty picks list, for consistency with other parsers
    picks = []
    return ssp_event, picks
