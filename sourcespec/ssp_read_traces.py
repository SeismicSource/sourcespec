# -*- coding: utf-8 -*-
"""
Read traces in multiple formats of data and metadata.

:copyright:
    2012 Claudio Satriano <satriano@ipgp.fr>

    2013-2014 Claudio Satriano <satriano@ipgp.fr>,
              Emanuela Matrullo <matrullo@geologie.ens.fr>

    2015-2021 Claudio Satriano <satriano@ipgp.fr>,
              Sophie Lambotte <sophie.lambotte@unistra.fr>
:license:
    CeCILL Free Software License Agreement, Version 2.1
    (http://www.cecill.info/index.en.html)
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import sys
import os
import io
import re
import logging
import warnings
import shutil
import tarfile
import tempfile
try:
    import cPickle as pickle
except ImportError:
    import pickle
import json
from datetime import datetime
from obspy.core import Stream, read, UTCDateTime
from obspy.core.util import AttribDict
from obspy.io.xseed import Parser
from obspy.io.xseed.utils import SEEDParserException
from obspy import read_inventory
from obspy.core.inventory import Inventory
from obspy.io.sac import attach_paz
from obspy.core import Trace
from obspy.geodetics import gps2dist_azimuth
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
    correct_traceids = _get_correct_traceids(traceid_file)
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


def _add_paz_and_coords(trace, dataless, paz_dict=None):
    trace.stats.paz = None
    trace.stats.coords = None
    traceid = trace.get_id()
    time = trace.stats.starttime
    # We first look into the dataless dictionary, if available
    if isinstance(dataless, dict):
        for sp in dataless.values():
            # Check first if our traceid is in the dataless file
            if traceid not in str(sp):
                continue
            try:
                paz = AttribDict(sp.get_paz(traceid, time))
                coords = AttribDict(sp.get_coordinates(traceid, time))
            except SEEDParserException as err:
                logger.error('%s time: %s' % (err, str(time)))
                pass
    elif isinstance(dataless, Inventory):
        try:
            with warnings.catch_warnings(record=True) as warns:
                # get_sacpz() can issue warnings on more than one PAZ found,
                # so let's catch those warnings and log them properly
                sacpz = dataless.get_response(traceid, time).get_sacpz()
                for w in warns:
                    message = str(w.message)
                    logger.warning('%s: %s' % (traceid, message))
            attach_paz(trace, io.StringIO(sacpz))
            paz = trace.stats.paz
            coords = AttribDict(dataless.get_coordinates(traceid, time))
        except Exception as err:
            logger.error('%s traceid: %s time: %s' % (err, traceid, str(time)))
            pass
    try:
        trace.stats.paz = paz
        # elevation is in meters in the dataless
        coords.elevation /= 1000.
        trace.stats.coords = coords
    except Exception:
        pass
    # If we couldn't find any PAZ in the dataless dictionary,
    # we try to attach paz from the paz dictionary passed
    # as argument
    if trace.stats.paz is None and paz_dict is not None:
        # Look for traceid or for a generic paz
        net, sta, loc, chan = trace.id.split('.')
        ids = [trace.id,
               '.'.join(('__', '__', '__', '__')),
               '.'.join((net, '__', '__', '__')),
               '.'.join((net, sta, '__', '__')),
               '.'.join((net, sta, loc, '__')),
               'default']
        for id in ids:
            try:
                paz = paz_dict[id]
                trace.stats.paz = paz
            except KeyError:
                pass
    # If we're still out of luck,
    # we try to build the sensitivity from the
    # user2 and user3 header fields (ISNet format)
    if trace.stats.paz is None and trace.stats.format == 'ISNet':
        try:
            # instrument constants
            u2 = trace.stats.sac.user2
            u3 = trace.stats.sac.user3
            paz = AttribDict()
            paz.sensitivity = u3/u2
            paz.poles = []
            paz.zeros = []
            paz.gain = 1
            trace.stats.paz = paz
        except AttributeError:
            pass
    # Still no paz? Antilles or IPOC format!
    if (trace.stats.paz is None and
        (trace.stats.format == 'Antilles' or
         trace.stats.format == 'IPOC')):
        paz = AttribDict()
        paz.sensitivity = 1
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
        raise Exception(
            '%s: could not find coords for trace: skipping trace' % traceid)


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
        if band_code in ['E', 'S']:
            instrtype = 'shortp'
        if band_code in ['B', 'C', 'H']:
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
    # we need to lazy-import here, so that OBSPY_VERSION is defined
    from sourcespec.ssp_setup import OBSPY_VERSION
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
            if OBSPY_VERSION > (1, 1, 1):
                # UTCDateTime objects will become immutable in future
                # versions of ObsPy
                _evid_time = origin_time.replace(
                    second=_second, microsecond=_microsecond)
            else:
                # For old versions, UTCDateTime objects are mutable
                _evid_time = UTCDateTime(origin_time)
                _evid_time.second = _second
                _evid_time.microsecond = _microsecond
        else:
            # make a copy of starttime and round it to the nearest minute
            _starttime = trace.stats.starttime
            _minute = _starttime.minute
            if _starttime.second >= 30:
                _minute += 1
            _second = 0
            _microsecond = 0
            if OBSPY_VERSION > (1, 1, 1):
                # UTCDateTime objects will become immutable in future
                # versions of ObsPy
                _evid_time = _starttime.replace(
                    minute=_minute, second=_second, microsecond=_microsecond)
            else:
                # For old versions, UTCDateTime objects are mutable
                _evid_time = UTCDateTime(_starttime)
                _evid_time.minute = _minute
                _evid_time.second = _second
                _evid_time.microsecond = _microsecond

        hypo = AttribDict()
        hypo.origin_time = origin_time
        hypo.evid = _evid_time.strftime('%Y%m%d_%H%M%S')
        hypo.latitude = evla
        hypo.longitude = evlo
        hypo.depth = evdp
    trace.stats.hypo = hypo
    _, _, baz = gps2dist_azimuth(
        hypo.latitude, hypo.longitude,
        trace.stats.coords.latitude, trace.stats.coords.longitude)
    trace.stats.back_azimuth = baz


def _add_picks(trace, picks):
    stat_picks = []
    station = trace.stats.station

    if picks is None:
        # try to get picks from SAC header
        if trace.stats._format != 'SAC':
            return

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
            stat_picks.append(pick)

    else:
        for pick in picks:
            if pick.station == station:
                stat_picks.append(pick)

    trace.stats.picks = stat_picks
    # Create an empty dict for arrivals.
    # It will be used later.
    trace.stats.arrivals = dict()


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
def _read_dataless(path):
    if path is None:
        return None

    # Try to read the file as StationXML
    if not os.path.isdir(path):
        try:
            inv = read_inventory(path)
            return inv
        except TypeError:
            pass
        except IOError as err:
            logger.error(err)
            ssp_exit()

    logger.info('Reading dataless...')
    dataless = dict()
    if os.path.isdir(path):
        listing = os.listdir(path)
        for filename in listing:
            fullpath = os.path.join(path, filename)
            try:
                dataless[filename] = Parser(fullpath)
            except Exception:
                continue
        # TODO: manage the case in which "path" is a file name
    logger.info('Reading dataless: done')
    return dataless


def _read_paz(path):
    """
    Read a directory with paz files or a single file.

    Limitations:
    (1) directory must contain *only* paz files
    (2) paz file can optionally have ".pz" or ".paz" suffixes
    (3) paz file name (without prefix and suffix) *has* to have
        the trace_id (NET.STA.LOC.CHAN) of the corresponding trace
        in the last part of his name
        (e.g., 20110208_1600.NOW.IV.CRAC.00.EHZ.paz)
    """
    if path is None:
        return None

    logger.info('Reading PAZ...')
    paz = dict()
    if os.path.isdir(path):
        listing = os.listdir(path)
        # check if files have a common prefix: we will strip it later
        prefix = os.path.commonprefix(listing)
        for filename in listing:
            fullpath = os.path.join(path, filename)
            try:
                # This is a horrible hack!
                # Since attach_paz needs a trace,
                # we create a trace and then, later,
                # we just retrieve the paz object
                # from the trace ;)
                tr = Trace()
                attach_paz(tr, fullpath)
                bname = os.path.basename(filename)
                # strip .pz suffix, if there
                bname = re.sub('.pz$', '', bname)
                # strip .paz suffix, if there
                bname = re.sub('.paz$', '', bname)
                # and strip any common prefix
                bname = re.sub('^' + prefix, '', bname)
                # we assume that the last four fields of bname
                # (separated by '.') are the trace_id
                trace_id = '.'.join(bname.split('.')[-4:])
                paz[trace_id] = tr.stats.paz.copy()
            except IOError:
                continue
    elif os.path.isfile(path):
        # If a filename is provided, store it as
        # 'default' paz.
        filename = path
        tr = Trace()
        attach_paz(tr, filename)
        paz['default'] = tr.stats.paz.copy()
    logger.info('Reading PAZ: done')
    return paz


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
    hypo.evid = ev.resource_id.id.split('/')[-1]

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


def _parse_hypocenter(hypo_file):
    if hypo_file is None:
        return None

    hypo = AttribDict()
    hypo.latitude = None
    hypo.longitude = None
    hypo.depth = None
    hypo.origin_time = None
    hypo.evid = None

    if isinstance(hypo_file, str):
        try:
            with open(hypo_file) as fp:
                # Corinth hypocenter file format:
                # TODO: check file format
                line = fp.readline()
                # Skip the first line if it contains
                # characters in the first 10 digits:
                if any(c.isalpha() for c in line[0:10]):
                    line = fp.readline()
        except IOError as err:
            logger.error(err)
            ssp_exit(1)

        timestr = line[0:17]
        # There are two possible formats for the timestring.
        # We try both of them
        try:
            dt = datetime.strptime(timestr, '%y%m%d %H %M%S.%f')
        except ValueError:
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

    else:  # FIXME: put a condition here!
        ev = hypo_file  # FIXME: improve this!
        hypo.latitude = ev.latitude
        hypo.longitude = ev.longitude
        hypo.depth = ev.depth
        hypo.origin_time = ev.utcdate
        hypo.evid = ev.event_id

    return hypo


def _is_hypo_format(fp):
    for line in fp.readlines():
        # remove newline
        line = line.replace('\n', '')
        # skip separator and empty lines
        stripped_line = line.strip()
        if stripped_line == '10' or stripped_line == '':
            continue
        # Check if it is a pick line
        # 6th character should be alpha (phase name: P or S)
        # other character should be digits (date/time)
        if (line[5].isalpha() and
                line[9].isdigit() and
                line[20].isdigit()):
            fp.seek(0)  # rewind file
            return True
        else:
            fp.seek(0)  # rewind file
            return False


# TODO: def _is_NLL_format(fp):


def _parse_picks(config):
    # we need to lazy-import here, so that OBSPY_VERSION is defined
    from sourcespec.ssp_setup import OBSPY_VERSION
    pick_file = config.options.pick_file
    if pick_file is None:
        return None

    try:
        fp = open(pick_file)
        if not _is_hypo_format(fp):
            raise Exception('%s: Not a phase file' % pick_file)
        lines = fp.readlines()
        fp.close()
    except Exception as err:
        logger.error(err)
        ssp_exit(1)
    picks = []
    for line in lines:
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
            _correct_station_name(pick.station, config.traceids)
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
        if OBSPY_VERSION > (1, 1, 1):
            # UTCDateTime objects will become immutable in future
            # versions of ObsPy
            pick2.time = pick.time.replace(second=0, microsecond=0)
        else:
            # For old versions, UTCDateTime objects are mutable
            pick2.time = UTCDateTime(pick.time)
            pick2.time.second = 0
            pick2.time.microsecond = 0
        # finally we add stime
        pick2.time += float(stime)
        picks.append(pick2)
    return picks


def _hypo_vel(hypo, config):
    vs = get_vel(hypo.longitude, hypo.latitude, hypo.depth,
                 'S', config.NLL_model_dir)
    if vs is None:
        vs = config.vs
    hypo.vs = vs


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


# PATH DISCOVERY --------------------------------------------------------------
# We try to guess the path of the hypo and pick file from the data dir
# This applies (for the moment) only to the Corinth format
def _set_hypo_file_path(config):
    if config.options.hypo_file is not None:
        return
    # try with the basename of the datadir
    if os.path.isdir(config.options.trace_path[0]):
        hypo_file = config.options.trace_path[0] + '.phs.h'
        try:
            open(hypo_file)
            config.options.hypo_file = hypo_file
        except Exception:
            pass
    return


def _set_pick_file_path(config):
    if config.options.pick_file is not None:
        return
    # try with the basename of the datadir
    if os.path.isdir(config.options.trace_path[0]):
        pick_file = config.options.trace_path[0] + '.phs'
        try:
            open(pick_file)
            config.options.pick_file = pick_file
        except Exception:
            pass
    return
# -----------------------------------------------------------------------------


# Public interface:
def read_traces(config):
    """Read traces, store waveforms and metadata."""
    # read dataless
    dataless = _read_dataless(config.dataless)
    # read PAZ (normally this is an alternative to dataless)
    paz = _read_paz(config.paz)

    # parse hypocenter file
    _set_hypo_file_path(config)
    hypo = _parse_hypocenter(config.options.hypo_file)
    # parse pick file
    _set_pick_file_path(config)
    picks = _parse_picks(config)

    # parse QML file
    if hypo is None:
        hypo, picks = _parse_qml(config.options.qml_file,
                                 config.options.evid)

    # finally, read traces
    # traces can be defined in a pickle catalog (Antilles format)...
    if config.pickle_catalog:
        sys.path.append(config.pickle_classpath)
        with open(config.pickle_catalog, 'rb') as fp:
            catalog = pickle.load(fp)
        event = [ev for ev in catalog if ev.event_id == config.options.evid][0]
        hypo = _parse_hypocenter(event)
        st = Stream()
        for trace in event.traces:
            if config.options.station is not None:
                if not trace.station == config.options.station:
                    continue
            try:
                tmpst = read(trace.trace_file, fsize=False)
            except Exception as err:
                logger.error(err)
                continue
            for trace in tmpst.traces:
                st.append(trace)
                trace.stats.format = config.trace_format
                _add_paz_and_coords(trace, dataless, paz)
                _add_hypocenter(trace, hypo)
                _add_picks(trace, picks)  # FIXME: actually add picks!
                # _add_instrtype(trace)
                trace.stats.instrtype = 'acc'  # FIXME
    # ...or in standard files (CRL and ISNet)
    else:
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
                trace.stats.format = config.trace_format
                _correct_traceid(trace, config.traceids)
                try:
                    _add_paz_and_coords(trace, dataless, paz)
                    _add_instrtype(trace, config)
                    _add_hypocenter(trace, hypo)
                    _add_picks(trace, picks)
                except Exception as err:
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
            logger.error('No hypocenter information found')
            ssp_exit()
    # add vs to hypo
    _hypo_vel(hypo, config)
    # add hypo to config file
    config.hypo = hypo

    st.sort()
    return st
