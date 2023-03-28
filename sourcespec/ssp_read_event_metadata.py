
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Read event metadata in QuakeML, HYPO71 or HYPOINVERSE format.

:copyright:
    2012-2023 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
import os
import contextlib
import logging
from datetime import datetime
from obspy import UTCDateTime
from obspy import read_events
from sourcespec.ssp_setup import ssp_exit, traceid_map
logger = logging.getLogger(__name__.split('.')[-1])


class Hypo():
    """A hypocenter object."""
    latitude = None
    longitude = None
    depth = None
    origin_time = None
    evid = None

    def __str__(self):
        return (
            f'latitude: {self.latitude}, longitude: {self.longitude}, '
            f'depth: {self.depth}, origin_time: {self.origin_time}, '
            f'evid: {self.evid}'
        )


class Pick():
    """A pick object."""
    station = None
    flag = None
    phase = None
    polarity = None
    quality = None
    time = None

    def __str__(self):
        return (
            f'station: {self.station}, flag: {self.flag}, '
            f'phase: {self.phase}, polarity: {self.polarity}, '
            f'quality: {self.quality}, time: {self.time}'
        )


def parse_qml(qml_file, evid=None):
    """Parse event metadata and picks from a QuakeML file."""
    hypo = None
    picks = []
    if qml_file is None:
        return hypo, picks
    try:
        ev = _get_event_from_qml(qml_file, evid)
        hypo = _get_hypo_from_event(ev)
        picks = _get_picks_from_event(ev)
    except Exception as err:
        logger.error(err)
        ssp_exit(1)
    # See if there is a focal mechanism with nodal planes
    with contextlib.suppress(Exception):
        _get_focal_mechanism_from_event(ev, hypo)
        logger.info('Found focal mechanism in QuakeML file')
    return hypo, picks


def _get_event_from_qml(qml_file, evid=None):
    cat = read_events(qml_file)
    # If evid is None, return the first event in the catalog
    # otherwise return the event with the given evid
    return (
        [e for e in cat if evid in str(e.resource_id)][0]
        if evid is not None
        else cat[0]
    )


def _get_hypo_from_event(ev):
    # See if there is a preferred origin...
    origin = ev.preferred_origin()
    # ...or just use the first one
    if origin is None:
        origin = ev.origins[0]
    hypo = Hypo()
    hypo.latitude = origin.latitude
    hypo.longitude = origin.longitude
    hypo.depth = origin.depth/1000.
    hypo.origin_time = origin.time
    hypo.evid = ev.resource_id.id.split('/')[-1].split('=')[-1]
    return hypo


def _get_focal_mechanism_from_event(ev, hypo):
    fm = ev.focal_mechanisms[0]
    nodal_plane = fm.nodal_planes.nodal_plane_1
    hypo.strike = nodal_plane.strike
    hypo.dip = nodal_plane.dip
    hypo.rake = nodal_plane.rake


def _get_picks_from_event(ev):
    picks = []
    for pck in ev.picks:
        pick = Pick()
        pick.station = pck.waveform_id.station_code
        pick.network = pck.waveform_id.network_code
        pick.channel = pck.waveform_id.channel_code
        if pck.waveform_id.location_code is not None:
            pick.location = pck.waveform_id.location_code
        else:
            pick.location = ''
        if pck.onset == 'emergent':
            pick.flag = 'E'
        elif pck.onset == 'impulsive':
            pick.flag = 'I'
        try:
            pick.phase = pck.phase_hint[:1]
        except Exception:
            # ignore picks with no phase hint
            continue
        if pck.polarity == 'negative':
            pick.polarity = 'D'
        elif pck.polarity == 'positive':
            pick.polarity = 'U'
        pick.time = pck.time
        picks.append(pick)
    return picks


def _parse_hypo71_hypocenter(hypo_file):
    with open(hypo_file) as fp:
        line = fp.readline()
        # Skip the first line if it contain characters in the first 10 digits:
        if any(c.isalpha() for c in line[:10]):
            line = fp.readline()
    hypo = Hypo()
    timestr = line[:17]
    # There are two possible formats for the timestring. We try both of them
    try:
        dt = datetime.strptime(timestr, '%y%m%d %H %M%S.%f')
    except Exception:
        dt = datetime.strptime(timestr, '%y%m%d %H%M %S.%f')
    hypo.origin_time = UTCDateTime(dt)
    lat = float(line[17:20])
    lat_deg = float(line[21:26])
    hypo.latitude = lat + lat_deg/60
    lon = float(line[26:30])
    lon_deg = float(line[31:36])
    hypo.longitude = lon + lon_deg/60
    hypo.depth = float(line[36:42])
    evid = os.path.basename(hypo_file)
    evid = evid.replace('.phs', '').replace('.h', '').replace('.hyp', '')
    hypo.evid = evid
    # empty picks list, for consistency with other parsers
    picks = []
    return hypo, picks


def _parse_hypo2000_hypo_line(line):
    word = line.split()
    hypo = Hypo()
    timestr = ' '.join(word[:3])
    hypo.origin_time = UTCDateTime(timestr)
    n = 3
    if word[n].isnumeric():
        # Check if word is integer
        # In this case the format should be: degrees and minutes
        latitude = float(word[n]) + float(word[n+1])/60.
        n += 2
    elif 'N' in word[n] or 'S' in word[n]:
        # Check if there is N or S in the string
        # In this case the format should be: degrees and minutes
        _word = word[n].replace('N', ' ').replace('S', ' ').split()
        latitude = float(_word[0]) + float(_word[1])/60.
        n += 1
    else:
        # Otherwise latitude should be in float format
        try:
            latitude = float(word[n])
        except Exception as e:
            raise ValueError(f'cannot read latitude: {word[n]}') from e
        n += 1
    hypo.latitude = latitude
    if word[n].isnumeric():
        # Check if word is integer
        # In this case the format should be: degrees and minutes
        longitude = float(word[n]) + float(word[n+1])/60.
        n += 2
    elif 'E' in word[n] or 'W' in word[n]:
        # Check if there is E or W in the string
        # In this case the format should be: degrees and minutes
        _word = word[n].replace('E', ' ').replace('W', ' ').split()
        longitude = float(_word[0]) + float(_word[1])/60.
        n += 1
    else:
        # Otherwise longitude should be in float format
        try:
            longitude = float(word[n])
        except Exception as e:
            raise ValueError(f'cannot read longitude: {word[n]}') from e
        n += 1
    hypo.longitude = longitude
    # depth is in km, according to the hypo2000 manual
    hypo.depth = float(word[n])
    return hypo


def _parse_hypo2000_station_line(line, oldpick, origin_time):
    if oldpick is not None:
        oldstation = oldpick.station
        oldnetwork = oldpick.network
        oldchannel = oldpick.channel
    else:
        oldstation = ''
        oldnetwork = ''
        oldchannel = ''
    pick = Pick()
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


def _parse_hypo2000_file(hypo_file):
    hypo = Hypo()
    picks = []
    hypo_line = False
    station_line = False
    oldpick = None
    for n, line in enumerate(open(hypo_file), start=1):
        word = line.split()
        if not word:
            continue
        # skip short lines, which are probably comments
        if len(line) < 50:
            continue
        if hypo_line:
            hypo = _parse_hypo2000_hypo_line(line)
            evid = os.path.basename(hypo_file)
            evid = evid.replace('.txt', '')
            hypo.evid = evid
        if station_line:
            try:
                pick = _parse_hypo2000_station_line(
                    line, oldpick, hypo.origin_time)
                oldpick = pick
                picks.append(pick)
            except Exception as err:
                logger.warning(f'Error parsing line {n} in {hypo_file}: {err}')
                continue
        if word[0] == 'YEAR':
            hypo_line = True
            continue
        hypo_line = False
        if word[0] == 'STA':
            station_line = True
    if not hypo:
        raise TypeError('Could not find hypocenter data.')
    return hypo, picks


def parse_hypo_file(hypo_file):
    """Parse a hypo71 or hypo2000 hypocenter file."""
    err_msgs = []
    parsers = {
        'hypo71': _parse_hypo71_hypocenter,
        'hypo2000': _parse_hypo2000_file,
    }
    for format, parser in parsers.items():
        try:
            hypo, picks = parser(hypo_file)
            return hypo, picks
        except Exception as err:
            msg = f'{hypo_file}: Not a {format} hypocenter file'
            err_msgs.append(msg)
            msg = f'Parsing error: {err}'
            err_msgs.append(msg)
    # If we arrive here, the file was not recognized as valid
    for msg in err_msgs:
        logger.error(msg)
    ssp_exit(1)


def _is_hypo71_picks(pick_file):
    for line in open(pick_file):
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
            msg = f'{pick_file}: Not a hypo71 phase file'
            raise TypeError(msg)


def _correct_station_name(station):
    """
    Correct station name, based on a traceid map.

    :param station: station name

    :return: corrected station name
    """
    if traceid_map is None:
        return station
    # get all the keys containing station name in it
    keys = [key for key in traceid_map if station == key.split('.')[1]]
    # then take just the first one
    try:
        key = keys[0]
    except IndexError:
        return station
    traceid = traceid_map[key]
    return traceid.split('.')[1]


def parse_hypo71_picks(config):
    """
    Parse hypo71 picks file

    :param config: Config object

    :return: list of Pick objects
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
    for line in open(pick_file):
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
        pick = Pick()
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
        pick2 = Pick()
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
