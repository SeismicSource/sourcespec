# -*- coding: utf-8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Read event metadata in QuakeML, SourceSpec Event File, HYPO71 or
HYPOINVERSE format.

:copyright:
    2012-2026 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
import os
import warnings
import contextlib
import logging
import re
import xml.etree.ElementTree as ET
from datetime import datetime
import yaml
from obspy import UTCDateTime
from obspy import read_events
from sourcespec.ssp_setup import ssp_exit, TRACEID_MAP
from sourcespec.ssp_event import SSPEvent
from sourcespec.ssp_pick import SSPPick
logger = logging.getLogger(__name__.rsplit('.', maxsplit=1)[-1])


def parse_qml(config):
    """Parse event metadata and picks from a QuakeML file."""
    ssp_event = None
    picks = []
    qml_file = config.options.qml_file
    event_id = config.options.evid
    # No need to parse event name from QuakeML if event name is given
    # in the command line
    if config.options.evname is None:
        qml_event_description = config.qml_event_description
        qml_event_description_regex = config.qml_event_description_regex
    else:
        qml_event_description = False
        qml_event_description_regex = None
    if qml_file is None:
        return ssp_event, picks
    try:
        qml_event = _get_event_from_qml(qml_file, event_id)
        ssp_event, qml_origin = _parse_qml_event(
            qml_event,
            parse_event_name_from_description=qml_event_description,
            event_description_regex=qml_event_description_regex)
        picks = _parse_picks_from_qml_event(qml_event, qml_origin)
    except Exception as err:
        logger.error(err)
        ssp_exit(1)
    log_messages = []
    with contextlib.suppress(Exception):
        _parse_moment_tensor_from_qml_event(qml_event, ssp_event)
        # Compute focal mechanism, scalar moment and magnitude.
        # They will be overwritten later on, if they are found in the
        # QuakeML file.
        ssp_event.focal_mechanism.from_moment_tensor(ssp_event.moment_tensor)
        ssp_event.scalar_moment.from_moment_tensor(ssp_event.moment_tensor)
        ssp_event.magnitude.from_scalar_moment(ssp_event.scalar_moment)
    with contextlib.suppress(Exception):
        _parse_scalar_moment_from_qml_event(qml_event, ssp_event)
        # Compute magnitude from scalar moment. It will be overwritten later
        # on, if it is found in the QuakeML file.
        ssp_event.magnitude.from_scalar_moment(ssp_event.scalar_moment)
    with contextlib.suppress(Exception):
        _parse_magnitude_from_qml_event(qml_event, ssp_event)
    with contextlib.suppress(Exception):
        _parse_focal_mechanism_from_qml_event(qml_event, ssp_event)
    for msg in log_messages:
        logger.info(msg)
    return ssp_event, picks


def _get_event_from_qml(qml_file, event_id=None):
    cat = read_events(qml_file)
    if event_id is not None:
        _qml_events = [ev for ev in cat if event_id in str(ev.resource_id)]
        try:
            qml_event = _qml_events[0]
        except IndexError as e:
            raise ValueError(
                f'Event {event_id} not found in {qml_file}') from e
    else:
        qml_event = cat[0]
        if len(cat) > 1:
            logger.warning(
                f'Found {len(cat)} events in {qml_file}. '
                'Using the first one.')
    return qml_event


def _get_evid_from_resource_id(resource_id):
    """
    Get evid from resource_id.

    :param resource_id: resource_id string
    :returns: evid string
    """
    evid = resource_id
    if '/' in evid:
        evid = resource_id.split('/')[-1]
    if '?' in evid:
        evid = resource_id.split('?')[-1]
    if '&' in evid:
        evid = evid.split('&')[0]
    if '=' in evid:
        evid = evid.split('=')[-1]
    return evid


def _parse_qml_event(
        qml_event,
        parse_event_name_from_description=False, event_description_regex=None):
    ssp_event = SSPEvent()
    ssp_event.event_id = _get_evid_from_resource_id(
        str(qml_event.resource_id.id))
    if parse_event_name_from_description:
        try:
            ssp_event.name = str(qml_event.event_descriptions[0].text)
        except IndexError:
            logger.warning(
                'QuakeML file does not contain an event description.')
    if ssp_event.name and event_description_regex:
        pattern = re.compile(event_description_regex)
        match = pattern.search(ssp_event.name)
        if match:
            name = match.group()
            # capitalize first letter
            name = name[0].upper() + name[1:]
            ssp_event.name = name
    # See if there is a preferred origin...
    origin = qml_event.preferred_origin()
    # ...or just use the first one
    if origin is None:
        origin = qml_event.origins[0]
    ssp_event.hypocenter.longitude.value_in_deg = origin.longitude
    ssp_event.hypocenter.latitude.value_in_deg = origin.latitude
    ssp_event.hypocenter.depth.value = origin.depth
    ssp_event.hypocenter.depth.units = 'm'
    ssp_event.hypocenter.origin_time = origin.time
    return ssp_event, origin


def _parse_magnitude_from_qml_event(qml_event, ssp_event):
    mag = qml_event.preferred_magnitude() or qml_event.magnitudes[0]
    ssp_event.magnitude.value = mag.mag
    ssp_event.magnitude.mag_type = mag.magnitude_type


def _parse_scalar_moment_from_qml_event(qml_event, ssp_event):
    fm = qml_event.preferred_focal_mechanism() or qml_event.focal_mechanisms[0]
    ssp_event.scalar_moment.value = fm.moment_tensor.scalar_moment
    ssp_event.scalar_moment.units = 'N-m'


def _parse_moment_tensor_from_qml_event(qml_event, ssp_event):
    fm = qml_event.preferred_focal_mechanism() or qml_event.focal_mechanisms[0]
    mt = fm.moment_tensor.tensor
    ssp_event.moment_tensor.m_rr = mt.m_rr
    ssp_event.moment_tensor.m_tt = mt.m_tt
    ssp_event.moment_tensor.m_pp = mt.m_pp
    ssp_event.moment_tensor.m_rt = mt.m_rt
    ssp_event.moment_tensor.m_rp = mt.m_rp
    ssp_event.moment_tensor.m_tp = mt.m_tp
    ssp_event.moment_tensor.units = 'N-m'


def _parse_focal_mechanism_from_qml_event(qml_event, ssp_event):
    fm = qml_event.focal_mechanisms[0]
    nodal_plane = fm.nodal_planes.nodal_plane_1
    ssp_event.focal_mechanism.strike = nodal_plane.strike
    ssp_event.focal_mechanism.dip = nodal_plane.dip
    ssp_event.focal_mechanism.rake = nodal_plane.rake


def _parse_picks_from_qml_event(ev, origin):
    picks = []
    for pck in ev.picks:
        pick = SSPPick()
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
        pick.phase = None
        if pck.phase_hint:
            pick.phase = pck.phase_hint[:1]
        else:
            # try to get the phase from the arrival object that uses this pick
            pick_id = pck.resource_id.id
            arrivals = [
                arr for arr in origin.arrivals if arr.pick_id.id == pick_id]
            if arrivals:
                pick.phase = arrivals[0].phase[:1]
        if pick.phase is None:
            # ignore picks with no phase hint
            continue
        if pck.polarity == 'negative':
            pick.polarity = 'D'
        elif pck.polarity == 'positive':
            pick.polarity = 'U'
        pick.time = pck.time
        picks.append(pick)
    return picks


def _parse_hypo71_hypocenter(hypo_file, _):
    """
    Parse a hypo71 hypocenter file.

    :param hypo_file: path to hypo71 hypocenter file
    :param _: unused (for consistency with other parsers)

    :return: a tuple of (SSPEvent, picks)
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


def _parse_hypo2000_hypo_line(line):
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


def _parse_hypo2000_file(hypo_file, _):
    """
    Parse a hypo2000 hypocenter file.

    :param hypo_file: path to hypo2000 hypocenter file
    :param _: unused (for consistency with other parsers)

    :return: a tuple of (SSPEvent, picks)
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


def _parse_source_spec_event_file(event_file, event_id=None):
    """
    Parse a SourceSpec Event File, which is a YAML file.

    :param event_file: path to SourceSpec event file
    :param evid: event id

    :return: SSPEvent object
    """
    try:
        with open(event_file, encoding='utf-8') as fp:
            events = yaml.safe_load(fp)
    except Exception as e:
        raise TypeError('Not a valid YAML file.') from e
    # XML is valid YAML, but obviously not a SourceSpec Event File
    try:
        root = ET.fromstring(events)
    except Exception:
        root = None
    if root:
        raise TypeError(
            'The file is an XML file, not a YAML file. '
            'Try reading it with the "-q" option (QuakeML format).')
    # raise TypeError if events is not a list
    if not isinstance(events, list):
        raise TypeError(
            'This is a valid YAML file, but it does not contain the key: '
            '"- event_id:", preceded by a dash (-).')
    # make sure that all the event_id are a string
    for ev in events:
        ev['event_id'] = str(ev['event_id'])
    if event_id is not None:
        _events = [ev for ev in events if ev.get('event_id') == event_id]
        try:
            event = _events[0]
        except IndexError as e:
            raise ValueError(
                f'Event {event_id} not found in {event_file}') from e
    else:
        event = events[0]
        if len(events) > 1:
            logger.warning(
                f'Found {len(events)} events in {event_file}. '
                'Using the first one.')
    try:
        with warnings.catch_warnings(record=True) as w:
            ssp_event = SSPEvent(event)
            if len(w) > 0:
                logger.warning(f'Warnings while parsing {event_file}:')
            for warning in w:
                logger.warning(warning.message)
    except Exception as e:
        raise TypeError(
            'This is a valid YAML file, but the following error occurred: '
            f'{e}.'
        ) from e
    # empty picks list, for consistency with other parsers
    picks = []
    return ssp_event, picks


# pylint: disable=inconsistent-return-statements
def parse_hypo_file(hypo_file, event_id=None):
    """
    Parse a SourceSpec Event File, hypo71 or hypo2000 hypocenter file.

    :param hypo_file:
        Path to the hypocenter file.
    :returns:
        A tuple of (SSPEvent, picks, format).
    """
    if not os.path.exists(hypo_file):
        logger.error(f'{hypo_file}: No such file or directory')
        ssp_exit(1)
    err_msgs = []
    parsers = {
        'ssp_event_file': _parse_source_spec_event_file,
        'hypo71': _parse_hypo71_hypocenter,
        'hypo2000': _parse_hypo2000_file,
    }
    format_strings = {
        'ssp_event_file': 'SourceSpec Event File',
        'hypo71': 'hypo71 hypocenter file',
        'hypo2000': 'hypo2000 hypocenter file',
    }
    for file_format, parser in parsers.items():
        try:
            ssp_event, picks = parser(hypo_file, event_id)
            return ssp_event, picks, file_format
        except Exception as err:
            format_str = format_strings[file_format]
            msg = f'{hypo_file}: Not a {format_str}'
            err_msgs.append(msg)
            msg = f'Parsing error: {err}'
            err_msgs.append(msg)
    # If we arrive here, the file was not recognized as valid
    for msg in err_msgs:
        logger.error(msg)
    ssp_exit(1)


def override_event_depth(ssp_event, depth_in_km):
    """
    Override event depth if specified in command line arguments.

    :param ssp_event: SSPEvent object
    :param depth_in_km: depth value in kilometers
    """
    if depth_in_km is not None:
        ssp_event.hypocenter.depth.value = depth_in_km
        ssp_event.hypocenter.depth.units = 'km'
        logger.info(f'Overriding event depth to {depth_in_km} km')


def _is_hypo71_picks(pick_file):
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

    :return: corrected station name
    """
    if TRACEID_MAP is None:
        return station
    # get all the keys containing station name in it
    keys = [key for key in TRACEID_MAP if station == key.split('.')[1]]
    # then take just the first one
    try:
        key = keys[0]
    except IndexError:
        return station
    traceid = TRACEID_MAP[key]
    return traceid.split('.')[1]


def parse_hypo71_picks(config):
    """
    Parse hypo71 picks file

    :param config: Config object

    :return: list of SSPPick objects
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
