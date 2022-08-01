# -*- coding: utf-8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Read traces in multiple formats of data and metadata.

:copyright:
    2012 Claudio Satriano <satriano@ipgp.fr>

    2013-2014 Claudio Satriano <satriano@ipgp.fr>,
              Emanuela Matrullo <matrullo@geologie.ens.fr>

    2015-2022 Claudio Satriano <satriano@ipgp.fr>,
              Sophie Lambotte <sophie.lambotte@unistra.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
import sys
import os
import io
import re
import logging
import warnings
import shutil
import tarfile
import tempfile
import json
from datetime import datetime
from obspy import read
from obspy.core import Stream, UTCDateTime
from obspy.core.util import AttribDict
from obspy import read_inventory
from obspy.io.sac import attach_paz
from obspy.core.inventory import Inventory, Network, Station, Channel, Response
from obspy import read_events
from sourcespec.ssp_setup import ssp_exit
from sourcespec.ssp_util import get_vel
logger = logging.getLogger(__name__.split('.')[-1])


class Pick(AttribDict):
    """A pick object."""

    def __init__(self):
        self.station = None
        self.flag = None
        self.phase = None
        self.polarity = None
        self.quality = None
        self.time = None


# TRACE MANIPULATION ----------------------------------------------------------
correct_traceids = None


def _get_correct_traceids(traceid_file):
    global correct_traceids
    if correct_traceids is None:
        with open(traceid_file, 'r') as fp:
            correct_traceids = json.loads(fp.read())
    return correct_traceids


def _correct_traceid(trace, traceid_file):
    if traceid_file is None:
        return
    try:
        correct_traceids = _get_correct_traceids(traceid_file)
    except Exception:
        msg = ('traceid mapping file "{}" not found '
               'or not in json format'.format(traceid_file))
        logger.error(msg)
        ssp_exit(1)
    try:
        traceid = correct_traceids[trace.get_id()]
        net, sta, loc, chan = traceid.split('.')
        trace.stats.network = net
        trace.stats.station = sta
        trace.stats.location = loc
        trace.stats.channel = chan
    except KeyError:
        pass


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


def _compute_sensitivity(trace, config):
    # Securize the string before calling eval()
    # see https://stackoverflow.com/a/25437733/2021880
    inp = re.sub(r'\.(?![0-9])', '', config.sensitivity)
    # Check if string contains letters, meaning that
    # it must contain SAC header fields and trace must be in SAC format
    namespace = None
    if re.search(r'[a-zA-Z]', inp):
        try:
            namespace = trace.stats.sac
        except Exception:
            raise Exception(
                '{}: trace must be in SAC format'.format(trace.id))
    try:
        sensitivity = eval(inp, {}, namespace)
    except NameError as msg:
        hdr_field = str(msg).split()[1]
        logger.error('SAC header field {} does not exist'.format(hdr_field))
        ssp_exit(1)
    return sensitivity


def _validate_coords(coords):
    """
    Return ``None`` if lon==lat==0 and depth==elevation==123456.

    Those values are given when reading an ``Inventory`` object from a
    RESP or PAZ file.
    """
    lon = coords.longitude
    lat = coords.latitude
    depth = coords.local_depth
    elevation = coords.elevation
    if lon == lat == 0 and depth == elevation == 123456:
        return None
    return coords


def _add_paz_and_coords(trace, inventory, config):
    traceid = trace.get_id()
    # If we already know that traceid is skipped, raise a silent exception
    if traceid in _add_paz_and_coords.skipped:
        raise Exception()
    trace.stats.paz = None
    trace.stats.coords = None
    time = trace.stats.starttime
    if isinstance(inventory, Inventory):
        try:
            with warnings.catch_warnings(record=True) as warns:
                # get_sacpz() can issue warnings on more than one PAZ found,
                # so let's catch those warnings and log them properly
                # warnings.filterwarnings('ignore', message='Found more than')
                # warnings.filterwarnings('ignore', message='More than')
                sacpz = inventory.get_response(traceid, time).get_sacpz()
                for w in warns:
                    msg = str(w.message)
                    logger.warning(
                        '{}: {} Time: {}'.format(traceid, msg, time))
            attach_paz(trace, io.StringIO(sacpz))
            paz = trace.stats.paz
            coords = AttribDict(inventory.get_coordinates(traceid, time))
            coords = _validate_coords(coords)
        except Exception as msg:
            logger.warning('{}: {} Time: {}'.format(traceid, msg, time))
            pass
    try:
        trace.stats.paz = paz
        # elevation is in meters
        coords.elevation /= 1000.
        trace.stats.coords = coords
    except Exception:
        pass
    # If a "sensitivity" config option is provided, override the paz computed
    # from the "Inventory" object (if any)
    if config.sensitivity is not None:
        # instrument constants
        paz = AttribDict()
        paz.sensitivity = _compute_sensitivity(trace, config)
        paz.poles = []
        paz.zeros = []
        paz.gain = 1
        trace.stats.paz = paz
    # If we still don't have trace coordinates,
    # we try to get them from SAC header
    if trace.stats.coords is None:
        try:
            stla = trace.stats.sac.stla
            stlo = trace.stats.sac.stlo
            try:
                stel = trace.stats.sac.stel
                # elevation is in meters in SAC header:
                stel /= 1000.
            except AttributeError:
                stel = 0.
            coords = AttribDict()
            coords.elevation = stel
            coords.latitude = stla
            coords.longitude = stlo
            trace.stats.coords = coords
        except AttributeError:
            pass
    # Still no coords? Raise an exception
    if trace.stats.coords is None:
        _add_paz_and_coords.skipped.append(traceid)
        raise Exception(
            '{}: could not find coords for trace: skipping trace'.format(
                traceid))
    if trace.stats.coords.latitude == trace.stats.coords.longitude == 0:
        logger.warning(
            '{}: trace has latitude and longitude equal to zero!'.format(
                traceid))


# list to keep track of skipped traces
_add_paz_and_coords.skipped = list()


def _add_instrtype(trace, config):
    instrtype = None
    trace.stats.instrtype = None

    # First, try to get the instrtype from channel name
    chan = trace.stats.channel
    if len(chan) > 2:
        band_code = chan[0]
        instr_code = chan[1]
    else:
        band_code = None
        instr_code = None
    # SEED standard instrument codes:
    # https://ds.iris.edu/ds/nodes/dmc/data/formats/seed-channel-naming/
    instr_codes_vel = ['H', 'L']
    instr_codes_acc = ['N', ]
    # User-defined instrument codes:
    instr_code_acc_user = config.instrument_code_acceleration
    instr_code_vel_user = config.instrument_code_velocity
    # Remove user-defined instrument codes if they conflict
    # with another instrument
    try:
        instr_codes_vel.remove(instr_code_acc_user)
    except ValueError:
        pass
    try:
        instr_codes_acc.remove(instr_code_vel_user)
    except ValueError:
        pass
    # Add user-defined instrument codes
    if instr_code_vel_user is not None:
        instr_codes_vel.append(instr_code_vel_user)
    if instr_code_acc_user is not None:
        instr_codes_acc.append(instr_code_acc_user)
    if instr_code in instr_codes_vel:
        # SEED standard band codes from higher to lower sampling rate
        # https://ds.iris.edu/ds/nodes/dmc/data/formats/seed-channel-naming/
        if band_code in ['G', 'D', 'E', 'S']:
            instrtype = 'shortp'
        if band_code in ['F', 'C', 'H', 'B']:
            instrtype = 'broadb'
    if instr_code in instr_codes_acc:
        instrtype = 'acc'

    # If, not possible, let's see if there is an instrument
    # name in "kinst" (ISNet format).
    # In this case, we define band and instrument codes
    # a posteriori.
    if instrtype is None:
        try:
            instr = trace.stats.sac.kinst
            if 'CMG-5T' in instr:
                instrtype = 'acc'
                band_code = 'H'
                instr_code = 'N'
            elif 'TRILLIUM' or 'CMG-40T' in instr:
                instrtype = 'broadb'
                band_code = 'H'
                instr_code = 'H'
            elif 'S13J' in instr:
                instrtype = 'shortp'
                band_code = 'S'
                instr_code = 'H'
            elif 'KS2000ED' in instr:
                instrtype = 'shortp'
                band_code = 'S'
                instr_code = 'H'
            else:
                return
        except AttributeError:
            return
        orientation = trace.stats.channel[-1]
        trace.stats.channel = ''.join((band_code, instr_code, orientation))
    trace.stats.instrtype = instrtype


def _add_hypocenter(trace, hypo):
    if hypo is None:
        # Try to get hypocenter information from the SAC header
        try:
            evla = trace.stats.sac.evla
            evlo = trace.stats.sac.evlo
            evdp = trace.stats.sac.evdp
            begin = trace.stats.sac.b
        except AttributeError:
            return

        try:
            tori = trace.stats.sac.o
            origin_time = trace.stats.starttime + tori - begin
        except AttributeError:
            origin_time = None

        if origin_time is not None:
            # make a copy of origin_time and round it to the nearest second
            _second = origin_time.second
            if origin_time.microsecond >= 500000:
                _second += 1
            _microsecond = 0
            _evid_time = origin_time.replace(
                second=_second, microsecond=_microsecond)
        else:
            # make a copy of starttime and round it to the nearest minute
            _starttime = trace.stats.starttime
            _minute = _starttime.minute
            if _starttime.second >= 30:
                _minute += 1
            _second = 0
            _microsecond = 0
            _evid_time = _starttime.replace(
                minute=_minute, second=_second, microsecond=_microsecond)

        hypo = AttribDict()
        hypo.origin_time = origin_time
        try:
            kevnm = trace.stats.sac.kevnm
            # if string is empty, raise Exception
            if not kevnm:
                raise Exception
            # if string has spaces, then kevnm is not a code,
            # so raise Exception
            if ' ' in kevnm:
                raise Exception
            hypo.evid = kevnm
        except Exception:
            hypo.evid = _evid_time.strftime('%Y%m%d_%H%M%S')
        hypo.latitude = evla
        hypo.longitude = evlo
        hypo.depth = evdp
    trace.stats.hypo = hypo


def _get_picks_from_SAC(trace):
    trace_picks = []
    station = trace.stats.station
    fields = ('a', 't0', 't1', 't2', 't3', 't4',
              't5', 't6', 't7', 't8', 't9')
    times = []
    labels = []
    for key in fields:
        try:
            times.append(trace.stats.sac[key])
        except KeyError:
            times.append(None)
        # now look at labels (ka, kt0, ...)
        key = 'k' + key
        try:
            labels.append(trace.stats.sac[key].strip())
        except KeyError:
            labels.append(None)
    for time, label, field in zip(times, labels, fields):
        if time is None:
            continue
        pick = Pick()
        pick.station = station
        begin = trace.stats.sac.b
        pick.time = trace.stats.starttime + time - begin
        if label is not None and len(label) == 4:
            pick.flag = label[0]
            pick.phase = label[1]
            pick.polarity = label[2]
            pick.quality = label[3]
        else:
            if field == 'a':
                pick.phase = 'P'
            elif field == 't0':
                pick.phase = 'S'
            else:
                pick.phase = 'X'
        trace_picks.append(pick)
    return trace_picks


def _add_picks(trace, picks):
    trace_picks = []
    station = trace.stats.station
    if picks is None:
        # try to get picks from SAC header
        if trace.stats._format == 'SAC':
            trace_picks = _get_picks_from_SAC(trace)
    else:
        for pick in picks:
            if pick.station == station:
                trace_picks.append(pick)
    trace.stats.picks = trace_picks
    # Create empty dicts for arrivals, travel_times and takeoff angles.
    # They will be used later.
    trace.stats.arrivals = dict()
    trace.stats.travel_times = dict()
    trace.stats.takeoff_angles = dict()


def _complete_picks(st):
    """Add component-specific picks to all components."""
    for station in set(tr.stats.station for tr in st):
        st_sel = st.select(station=station)
        # 'code' is band+instrument code
        for code in set(tr.stats.channel[:-1] for tr in st_sel):
            st_sel2 = st_sel.select(channel=code + '?')
            # Select default P and S picks as the first in list
            all_picks = [pick for tr in st_sel2 for pick in tr.stats.picks]
            default_P_pick = [pick for pick in all_picks
                              if pick.phase == 'P'][0:1]
            default_S_pick = [pick for pick in all_picks
                              if pick.phase == 'S'][0:1]
            for tr in st_sel2:
                # Attribute default picks to components without picks
                if len([pick for pick in tr.stats.picks
                        if pick.phase == 'P']) == 0:
                    tr.stats.picks += default_P_pick
                if len([pick for pick in tr.stats.picks
                        if pick.phase == 'S']) == 0:
                    tr.stats.picks += default_S_pick
# -----------------------------------------------------------------------------


# FILE PARSING ----------------------------------------------------------------
def _read_metadata(path):
    """
    Read metadata into an ObsPy ``Inventory`` object.
    """
    if path is None:
        return None
    logger.info('Reading station metadata...')
    inventory = Inventory()
    if os.path.isdir(path):
        filelist = [os.path.join(path, file) for file in os.listdir(path)]
    else:
        filelist = [path, ]
    for file in sorted(filelist):
        if os.path.isdir(file):
            # we do not enter into subdirs of "path"
            continue
        logger.info('Reading station metadata from file: {}'.format(file))
        try:
            inventory += read_inventory(file)
        except Exception:
            msg1 = 'Unable to parse file "{}" as Inventory'.format(file)
            try:
                inventory += _read_paz_file(file)
            except Exception as msg2:
                logger.warning(msg1)
                logger.warning(msg2)
                continue
    if not inventory:
        inventory = None
    logger.info('Reading station metadata: done')
    return inventory


def _parse_paz_file(file):
    lines = iter(open(file, 'r'))
    zeros = []
    poles = []
    constant = None
    linenumber = 0
    try:
        for line in lines:
            linenumber += 1
            word = line.split()
            if word[0] == 'ZEROS':
                nzeros = int(word[1])
                for _ in range(nzeros):
                    linenumber += 1
                    _zero = complex(*map(float, next(lines).split()))
                    zeros.append(_zero)
            if word[0] == 'POLES':
                npoles = int(word[1])
                for _ in range(npoles):
                    linenumber += 1
                    _pole = complex(*map(float, next(lines).split()))
                    poles.append(_pole)
            if word[0] == 'CONSTANT':
                constant = float(word[1])
    except Exception:
        msg = 'Unable to parse file "{}" as PAZ. '.format(file)
        msg += 'Parse error at line {}'.format(linenumber)
        raise TypeError(msg)
    if constant is None:
        msg = 'Unable to parse file "{}" as PAZ. '.format(file)
        msg += 'Cannot find a "CONSTANT" value'
        raise TypeError(msg)
    return zeros, poles, constant


def _read_paz_file(file):
    """
    Read a paz file into an ``Inventory``object.

    Limitations:
    (1) paz file must have ".pz" or ".paz" suffix (or no suffix)
    (2) paz file name (without prefix and suffix) *has* to have
        the trace_id (NET.STA.LOC.CHAN) of the corresponding trace
        in the last part of his name
        (e.g., 20110208_1600.NOW.IV.CRAC.00.EHZ.paz)
    """
    bname = os.path.basename(file)
    # strip .pz suffix, if there
    bname = re.sub('.pz$', '', bname)
    # strip .paz suffix, if there
    bname = re.sub('.paz$', '', bname)
    # we assume that the last four fields of bname
    # (separated by '.') are the trace_id
    trace_id = '.'.join(bname.split('.')[-4:])
    zeros, poles, constant = _parse_paz_file(file)
    resp = Response().from_paz(
        zeros, poles, stage_gain=1, input_units='M/S', output_units='COUNTS')
    resp.instrument_sensitivity.value = constant
    net, sta, loc, chan = trace_id.split('.')
    channel = Channel(
        code=chan, location_code=loc, response=resp,
        latitude=0, longitude=0, elevation=123456, depth=123456)
    station = Station(
        code=sta, channels=[channel, ],
        latitude=0, longitude=0, elevation=123456)
    network = Network(code=net, stations=[station, ])
    inv = Inventory(networks=[network, ])
    return inv


def _parse_qml(qml_file, evid=None):
    if qml_file is None:
        return None, None

    hypo = AttribDict()
    hypo.latitude = None
    hypo.longitude = None
    hypo.depth = None
    hypo.origin_time = None
    hypo.evid = None

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
    hypo = AttribDict()
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
    hypo = AttribDict()
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


def _parse_hypo2000_station_line(line, oldpick, orig_time):
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
        raise Exception()
    seconds = float(line[37:43])
    time = orig_time.replace(second=0, microsecond=0)
    pick.time = time + seconds
    return pick


def _parse_hypo2000_file(hypo_file):
    hypo = AttribDict()
    picks = list()
    hypo_line = False
    station_line = False
    oldpick = None
    for line in open(hypo_file):
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
                    line, oldpick, hypo.orig_time)
                oldpick = pick
                picks.append(pick)
            except Exception:
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


def _parse_hypo_file(hypo_file):
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


def _parse_hypo71_picks(config):
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


def _hypo_vel(hypo, config):
    hypo.vp = get_vel(hypo.longitude, hypo.latitude, hypo.depth, 'P', config)
    hypo.vs = get_vel(hypo.longitude, hypo.latitude, hypo.depth, 'S', config)


def _build_filelist(path, filelist, tmpdir):
    if os.path.isdir(path):
        listing = os.listdir(path)
        for filename in listing:
            fullpath = os.path.join(path, filename)
            _build_filelist(fullpath, filelist, tmpdir)
    else:
        try:
            open(path)
        except IOError as err:
            logger.error(err)
            return
        if tarfile.is_tarfile(path) and tmpdir is not None:
            tar = tarfile.open(path)
            tar.extractall(path=tmpdir)
            tar.close()
        else:
            filelist.append(path)
# -----------------------------------------------------------------------------


# Public interface:
def read_traces(config):
    """Read traces, store waveforms and metadata."""
    # read metadata into an ObsPy ``Inventory`` object
    inventory = _read_metadata(config.station_metadata)

    hypo = picks = None
    # parse hypocenter file
    if config.options.hypo_file is not None:
        hypo, picks = _parse_hypo_file(config.options.hypo_file)
    # parse pick file
    if config.options.pick_file is not None:
        picks = _parse_hypo71_picks(config)
    # parse QML file
    if config.options.qml_file is not None:
        hypo, picks = _parse_qml(config.options.qml_file, config.options.evid)

    # finally, read traces
    logger.info('Reading traces...')
    # phase 1: build a file list
    # ph 1.1: create a temporary dir and run '_build_filelist()'
    #         to move files to it and extract all tar archives
    tmpdir = tempfile.mkdtemp()
    filelist = []
    for trace_path in config.options.trace_path:
        _build_filelist(trace_path, filelist, tmpdir)
    # ph 1.2: rerun '_build_filelist()' in tmpdir to add to the
    #         filelist all the extraceted files
    listing = os.listdir(tmpdir)
    for filename in listing:
        fullpath = os.path.join(tmpdir, filename)
        _build_filelist(fullpath, filelist, None)
    # phase 2: build a stream object from the file list
    orientation_codes = config.vertical_channel_codes +\
        config.horizontal_channel_codes_1 +\
        config.horizontal_channel_codes_2
    st = Stream()
    for filename in sorted(filelist):
        try:
            tmpst = read(filename, fsize=False)
        except Exception:
            logger.warning('%s: Unable to read file as a trace: '
                           'skipping' % filename)
            continue
        for trace in tmpst.traces:
            orientation = trace.stats.channel[-1]
            if orientation not in orientation_codes:
                logger.warning('%s: Unknown channel orientation: "%s": '
                               'skipping trace' % (trace.id, orientation))
                continue
            if config.options.station is not None:
                if not trace.stats.station == config.options.station:
                    continue
            _correct_traceid(trace, config.traceid_mapping_file)
            try:
                _add_paz_and_coords(trace, inventory, config)
                _add_instrtype(trace, config)
                _add_hypocenter(trace, hypo)
                _add_picks(trace, picks)
            except Exception as err:
                # only warn if error message is not empty
                if str(err):
                    logger.warning(err)
                continue
            st.append(trace)
    shutil.rmtree(tmpdir)

    logger.info('Reading traces: done')
    if len(st.traces) == 0:
        logger.info('No trace loaded')
        ssp_exit()

    _complete_picks(st)

    # if hypo is still None, get it from first trace
    if hypo is None:
        try:
            hypo = st[0].stats.hypo
        except AttributeError:
            logger.error('No hypocenter information found.')
            sys.stderr.write(
                '\n'
                'Use "-q" or "-H" options to provide hypocenter information\n'
                'or add hypocenter information to the SAC file header\n'
                '(if you use the SAC format).\n'
            )
            ssp_exit()
    # add vs to hypo
    _hypo_vel(hypo, config)
    # add hypo to config file
    config.hypo = hypo

    st.sort()
    return st
