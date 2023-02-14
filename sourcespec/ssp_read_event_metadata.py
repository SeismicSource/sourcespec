
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
import logging
from datetime import datetime
from obspy import UTCDateTime
from obspy import read_events
from sourcespec.ssp_setup import ssp_exit
logger = logging.getLogger(__name__.split('.')[-1])


class Hypo():
    """A hypocenter object."""
    latitude = None
    longitude = None
    depth = None
    origin_time = None
    evid = None


class Pick():
    """A pick object."""
    station = None
    flag = None
    phase = None
    polarity = None
    quality = None
    time = None


def parse_qml(qml_file, evid=None):
    if qml_file is None:
        return None, None
    hypo = Hypo()
    try:
        cat = read_events(qml_file)
    except Exception as err:
        logger.error(err)
        ssp_exit(1)
    if evid is not None:
        ev = [e for e in cat if evid in str(e.resource_id)][0]
    else:
        # just take the first event
        ev = cat[0]
    # See if there is a preferred origin...
    origin = ev.preferred_origin()
    # ...or just use the first one
    if origin is None:
        origin = ev.origins[0]
    hypo.origin_time = origin.time
    hypo.latitude = origin.latitude
    hypo.longitude = origin.longitude
    hypo.depth = origin.depth/1000.
    hypo.evid = ev.resource_id.id.split('/')[-1].split('=')[-1]

    # See if there is a focal mechanism with nodal planes
    try:
        fm = ev.focal_mechanisms[0]
        nodal_plane = fm.nodal_planes.nodal_plane_1
        hypo.strike = nodal_plane.strike
        hypo.dip = nodal_plane.dip
        hypo.rake = nodal_plane.rake
        logger.info('Found focal mechanism in QuakeML file')
    except Exception:
        pass

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
            pick.phase = pck.phase_hint[0:1]
        except Exception:
            # ignore picks with no phase hint
            continue
        if pck.polarity == 'positive':
            pick.polarity = 'U'
        elif pck.polarity == 'negative':
            pick.polarity = 'D'
        pick.time = pck.time
        picks.append(pick)

    return hypo, picks


def _parse_hypo71_hypocenter(hypo_file):
    with open(hypo_file) as fp:
        line = fp.readline()
        # Skip the first line if it contains
        # characters in the first 10 digits:
        if any(c.isalpha() for c in line[0:10]):
            line = fp.readline()
    hypo = Hypo()
    timestr = line[0:17]
    # There are two possible formats for the timestring.
    # We try both of them
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
    return hypo


def _parse_hypo2000_hypo_line(line):
    word = line.split()
    hypo = Hypo()
    timestr = ' '.join(word[0:3])
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
        except Exception:
            msg = 'cannot read latitude: {}'.format(word[n])
            raise Exception(msg)
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
        except Exception:
            msg = 'cannot read longitude: {}'.format(word[n])
            raise Exception(msg)
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
        oldstation = ""
        oldnetwork = ""
        oldchannel = ""
    pick = Pick()
    station = line[1:5].strip()
    pick.station = station if station else oldstation
    oldstation = pick.station
    network = line[6:8].strip()
    pick.network = network if network else oldnetwork
    oldnetwork = pick.network
    channel = line[9:12].strip()
    pick.channel = channel if channel else oldchannel
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
    picks = list()
    hypo_line = False
    station_line = False
    oldpick = None
    for n, line in enumerate(open(hypo_file)):
        word = line.split()
        if not word:
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
    picks = None
    err_msgs = []
    try:
        hypo = _parse_hypo71_hypocenter(hypo_file)
        return hypo, picks
    except Exception as err:
        msg = '{}: Not a hypo71 hypocenter file'.format(hypo_file)
        err_msgs.append(msg)
        msg = 'Parsing error: ' + str(err)
        err_msgs.append(msg)
    try:
        hypo, picks = _parse_hypo2000_file(hypo_file)
        return hypo, picks
    except Exception as err:
        msg = '{}: Not a hypo2000 hypocenter file'.format(hypo_file)
        err_msgs.append(msg)
        msg = 'Parsing error: ' + str(err)
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
        if stripped_line == '10' or stripped_line == '':
            continue
        # Check if it is a pick line
        # 6th character should be alpha (phase name: P or S)
        # other character should be digits (date/time)
        if not (line[5].isalpha() and
                line[9].isdigit() and
                line[20].isdigit()):
            msg = '{}: Not a hypo71 phase file'.format(pick_file)
            raise Exception(msg)


def _correct_station_name(station, traceid_file):
    if traceid_file is None:
        return station
    correct_traceids = _get_correct_traceids(traceid_file)
    # get all the keys containing station name in it
    keys = [key for key in correct_traceids
            if station == key.split('.')[1]]
    # then take just the first one
    try:
        key = keys[0]
    except IndexError:
        return station
    traceid = correct_traceids[key]
    correct_station = traceid.split('.')[1]
    return correct_station


def parse_hypo71_picks(config):
    pick_file = config.options.pick_file
    if pick_file is None:
        return None

    try:
        _is_hypo71_picks(pick_file)
    except Exception as err:
        logger.error(err)
        ssp_exit(1)
    picks = []
    for line in open(pick_file):
        # remove newline
        line = line.replace('\n', '')
        # skip separator and empty lines
        stripped_line = line.strip()
        if stripped_line == '10' or stripped_line == '':
            continue
        # Check if it is a pick line
        # 6th character should be alpha (phase name: P or S)
        # other character should be digits (date/time)
        if not (line[5].isalpha() and
                line[9].isdigit() and
                line[20].isdigit()):
            continue

        pick = Pick()
        pick.station = line[0:4].strip()
        pick.station = \
            _correct_station_name(pick.station, config.traceid_mapping_file)
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
