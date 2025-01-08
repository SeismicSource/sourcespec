# -*- coding: utf-8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Read event and phase picks from an hypo2000 file.

:copyright:
    2012-2025 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
import os
import logging
from obspy import UTCDateTime
from ...ssp_event import SSPEvent
from ...ssp_pick import SSPPick
logger = logging.getLogger(__name__.rsplit('.', maxsplit=1)[-1])


def _parse_hypo2000_hypo_line(line):
    """
    Parse a line from a hypo2000 hypocenter file.

    :param line: line from hypo2000 hypocenter file
    :type line: str

    :return: SSPEvent object
    :rtype: SSPEvent
    """
    word = line.split()
    ssp_event = SSPEvent()
    hypo = ssp_event.hypocenter
    timestr = ' '.join(word[:3])
    hypo.origin_time = UTCDateTime(timestr)
    n = 3
    if word[n].isnumeric():
        # Check if word is integer
        # In this case the format should be: degrees and minutes
        latitude = float(word[n]) + float(word[n + 1]) / 60.
        n += 2
    elif 'N' in word[n] or 'S' in word[n]:
        # Check if there is N or S in the string
        # In this case the format should be: degrees and minutes
        _word = word[n].replace('N', ' ').replace('S', ' ').split()
        latitude = float(_word[0]) + float(_word[1]) / 60.
        n += 1
    else:
        # Otherwise latitude should be in float format
        try:
            latitude = float(word[n])
        except Exception as e:
            raise ValueError(f'cannot read latitude: {word[n]}') from e
        n += 1
    hypo.latitude.value_in_deg = latitude
    if word[n].isnumeric():
        # Check if word is integer
        # In this case the format should be: degrees and minutes
        longitude = float(word[n]) + float(word[n + 1]) / 60.
        n += 2
    elif 'E' in word[n] or 'W' in word[n]:
        # Check if there is E or W in the string
        # In this case the format should be: degrees and minutes
        _word = word[n].replace('E', ' ').replace('W', ' ').split()
        longitude = float(_word[0]) + float(_word[1]) / 60.
        n += 1
    else:
        # Otherwise longitude should be in float format
        try:
            longitude = float(word[n])
        except Exception as e:
            raise ValueError(f'cannot read longitude: {word[n]}') from e
        n += 1
    hypo.longitude.value_in_deg = longitude
    hypo.depth.value = float(word[n])
    # depth is in km, according to the hypo2000 manual
    hypo.depth.units = 'km'
    return ssp_event


def _parse_hypo2000_station_line(line, oldpick, origin_time):
    """
    Parse a line from a hypo2000 station file.

    :param line: line from hypo2000 station file
    :type line: str
    :param oldpick: previous pick
    :type oldpick: SSPPick
    :param origin_time: origin time
    :type origin_time: UTCDateTime

    :return: SSPPick object
    :rtype: SSPPick
    """
    if oldpick is not None:
        oldstation = oldpick.station
        oldnetwork = oldpick.network
        oldchannel = oldpick.channel
    else:
        oldstation = ''
        oldnetwork = ''
        oldchannel = ''
    pick = SSPPick()
    station = line[1:5].strip()
    pick.station = station or oldstation
    oldstation = pick.station
    network = line[6:8].strip()
    pick.network = network or oldnetwork
    oldnetwork = pick.network
    channel = line[9:12].strip()
    pick.channel = channel or oldchannel
    oldchannel = pick.channel
    # pick.flag = line[4:5]
    pick.phase = line[31:34].strip()
    if not pick.phase:
        raise ValueError('Cannot read pick phase')
    seconds = float(line[37:43])
    time = origin_time.replace(second=0, microsecond=0)
    pick.time = time + seconds
    return pick


def parse_hypo2000_file(hypo_file, _):
    """
    Parse a hypo2000 hypocenter file.

    :param hypo_file: path to hypo2000 hypocenter file
    :type hypo_file: str
    :param _: unused (for consistency with other parsers)
    :type _: None

    :return: a tuple of (SSPEvent, picks)
    :rtype: tuple
    """
    ssp_event = None
    picks = []
    hypo_line = False
    station_line = False
    oldpick = None
    with open(hypo_file, encoding='ascii') as fp:
        for n, line in enumerate(fp, start=1):
            word = line.split()
            if not word:
                continue
            # skip short lines, which are probably comments
            if len(line) < 50:
                continue
            if hypo_line:
                ssp_event = _parse_hypo2000_hypo_line(line)
                evid = os.path.basename(hypo_file)
                evid = evid.replace('.txt', '')
                ssp_event.event_id = evid
            if station_line and not ssp_event:
                raise TypeError('Could not find hypocenter data.')
            if station_line:
                try:
                    pick = _parse_hypo2000_station_line(
                        line, oldpick, ssp_event.hypocenter.origin_time)
                    oldpick = pick
                    picks.append(pick)
                except Exception as err:
                    logger.warning(
                        f'Error parsing line {n} in {hypo_file}: {err}')
                    continue
            if word[0] == 'YEAR':
                hypo_line = True
                continue
            hypo_line = False
            if word[0] == 'STA':
                station_line = True
    if not ssp_event:
        raise TypeError('Could not find hypocenter data.')
    return ssp_event, picks
