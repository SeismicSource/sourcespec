# -*- coding: utf-8 -*-
# ssp_read_traces.py
#
# All the functions whose name is between "__" are intended to be private
# (c) 2012 Claudio Satriano <satriano@ipgp.fr>
# (c) 2013-2014 Claudio Satriano <satriano@ipgp.fr>,
#               Emanuela Matrullo <matrullo@geologie.ens.fr>
'''
Read traces in multiple formats of data and metadata.
'''
from __future__ import division
import sys
import os
import re
import logging
import shutil
import tarfile
import tempfile
import cPickle as pickle
from datetime import datetime
from imp import load_source
from obspy.core import Stream, read, UTCDateTime
from obspy.core.util import AttribDict
from obspy.xseed import Parser
from obspy.xseed.utils import SEEDParserException
from ssp_setup import ssp_exit


class Pick(AttribDict):
    '''
    A pick object.
    '''
    def __init__(self):
        self.station = None
        self.flag = None
        self.phase = None
        self.polarity = None
        self.quality = None
        self.time = None


# TRACE MANIPULATION ----------------------------------------------------------
correct_traceids = None

def __correct_traceid__(trace, traceid_file):
    if traceid_file is None:
        return
    global correct_traceids
    if correct_traceids is None:
        #FIXME: load_source is not secure!
        trids = load_source('traceids', traceid_file)
        correct_traceids = trids.__correct_traceid_dict__
    try:
        traceid = correct_traceids[trace.getId()]
        net, sta, loc, chan = traceid.split('.')
        trace.stats.network = net
        trace.stats.station = sta
        trace.stats.location = loc
        trace.stats.channel = chan
    except KeyError:
        pass


def __correct_station_name__(station, traceid_file):
    if traceid_file is None:
        return
    global correct_traceids
    if correct_traceids is None:
        #FIXME: load_source is not secure!
        trids = load_source('traceids', traceid_file)
        correct_traceids = trids.__correct_traceid_dict__
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


def __add_paz_and_coords__(trace, dataless, paz_dict=None):
    trace.stats.paz = None
    trace.stats.coords = None
    traceid = trace.getId()
    time = trace.stats.starttime
    # We first look into the dataless dictionary, if available
    if dataless != None:
        for sp in dataless.values():
            # Check first if our traceid is in the dataless file
            if not traceid in str(sp):
                continue
            try:
                paz = sp.getPAZ(traceid, time)
                coords = AttribDict(sp.getCoordinates(traceid, time))
                # elevation is in meters in the dataless
                coords.elevation /= 1000.
            except SEEDParserException, message:
                logging.error("%s time: %s" % (message, str(time)))
                pass
    try:
        trace.stats.paz = paz
        trace.stats.coords = coords
    except:
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
            if u2 == -12345 or u3 == -12345:
                raise AttributeError
            paz = AttribDict()
            paz.sensitivity = u3/u2
            paz.poles = []
            paz.zeros = []
            paz.gain = 1
            trace.stats.paz = paz
        except AttributeError:
            pass
    # Still no paz? Antilles or IPOC format!
    if trace.stats.paz is None\
            and (trace.stats.format == 'Antilles'
                 or trace.stats.format == 'IPOC'):
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
            stel = trace.stats.sac.stel
            if stla == -12345 or stlo == -12345:
                raise AttributeError
            if stel == -12345:
                stel = 0.
            # elevation is in meters in SAC header:
            stel /= 1000.
            coords = AttribDict()
            coords.elevation = stel
            coords.latitude = stla
            coords.longitude = stlo
            trace.stats.coords = coords
        except AttributeError:
            pass


def __add_instrtype__(trace):
    instrtype = None
    trace.stats.instrtype = None

    # First, try to get the instrtype from channel name
    chan = trace.stats.channel
    try:
        band_code = chan[0]
    except IndexError:
        band_code = None
    try:
        instr_code = chan[1]
    except IndexError:
        instr_code = None
    if instr_code == 'H' or instr_code == 'L':
        if band_code == 'E' or band_code == 'S':
            instrtype = 'shortp'
        if band_code == 'H':
            instrtype = 'broadb'
    if instr_code == 'N' or instr_code =='L':
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


def __add_hypocenter__(trace, hypo):
    if hypo is None:
        # Try to get hypocenter information from the SAC header
        try:
            evla = trace.stats.sac.evla
            evlo = trace.stats.sac.evlo
            evdp = trace.stats.sac.evdp
            tori = trace.stats.sac.o
            begin = trace.stats.sac.b
            if (evla == -12345 or evlo == -12345
                    or evdp == -12345 or begin == -12345):
                raise AttributeError

            hypo = AttribDict()
            if tori == -12345:
                hypo.origin_time = None
                hypo.evid = trace.stats.starttime.strftime("%Y%m%d_%H%M%S")
            else:
                hypo.origin_time = trace.stats.starttime + tori - begin
                hypo.evid = hypo.origin_time.strftime("%Y%m%d_%H%M%S")
            hypo.latitude = evla
            hypo.longitude = evlo
            hypo.depth = evdp
        except AttributeError:
            return
    trace.stats.hypo = hypo


def __add_picks__(trace, picks):
    stat_picks = []
    station = trace.stats.station

    if picks is None:
        # try to get picks from SAC header
        if trace.stats._format != 'SAC':
            return

        times = []
        times.append(trace.stats.sac.a)
        times.append(trace.stats.sac.t0)
        times.append(trace.stats.sac.t1)
        times.append(trace.stats.sac.t2)
        times.append(trace.stats.sac.t3)
        times.append(trace.stats.sac.t4)
        times.append(trace.stats.sac.t5)
        times.append(trace.stats.sac.t6)
        times.append(trace.stats.sac.t7)
        times.append(trace.stats.sac.t8)
        times.append(trace.stats.sac.t9)
        labels = []
        labels.append(trace.stats.sac.ka.strip())
        labels.append(trace.stats.sac.kt0.strip())
        labels.append(trace.stats.sac.kt1.strip())
        labels.append(trace.stats.sac.kt2.strip())
        labels.append(trace.stats.sac.kt3.strip())
        labels.append(trace.stats.sac.kt4.strip())
        labels.append(trace.stats.sac.kt5.strip())
        labels.append(trace.stats.sac.kt6.strip())
        labels.append(trace.stats.sac.kt7.strip())
        labels.append(trace.stats.sac.kt8.strip())
        labels.append(trace.stats.sac.kt9.strip())
        fields = ['a', 't0', 't1', 't2', 't3', 't4',
                  't5', 't6', 't7', 't8', 't9']

        for time, label, field in zip(times, labels, fields):
            if time == -12345:
                continue

            pick = Pick()
            pick.station = station
            begin = trace.stats.sac.b
            pick.time = trace.stats.starttime + time - begin
            if len(label) == 4:
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


def __complete_picks__(st):
    '''
    For every station/instrument, adds component-specific picks
    to all components.
    '''
    for station in set(tr.stats.station for tr in st):
        st_sel = st.select(station=station)
        # 'code' is band+instrument code
        for code in set(tr.stats.channel[0:2] for tr in st_sel):
            st_sel2 = st_sel.select(channel=code + '?')
            # Select default P and S picks as the first in list
            all_picks = [pick for tr in st_sel2 for pick in tr.stats.picks]
            default_P_pick = [pick for pick in all_picks if pick.phase == 'P'][0:1]
            default_S_pick = [pick for pick in all_picks if pick.phase == 'S'][0:1]
            for tr in st_sel2:
                # Attribute default picks to components without picks
                if len([pick for pick in tr.stats.picks if pick.phase == 'P']) == 0:
                    tr.stats.picks += default_P_pick
                if len([pick for pick in tr.stats.picks if pick.phase == 'S']) == 0:
                    tr.stats.picks += default_S_pick
# -----------------------------------------------------------------------------



# FILE PARSING ----------------------------------------------------------------
def __read_dataless__(path):
    if path is None:
        return None

    logging.info('Reading dataless...')
    dataless = dict()
    if os.path.isdir(path):
        listing = os.listdir(path)
        for filename in listing:
            fullpath = os.path.join(path, filename)
            try:
                sp = Parser(fullpath)
                dataless[filename] = sp
            except IOError:
                continue
        #TODO: manage the case in which "path" is a file name
    logging.info('Reading dataless: done')
    return dataless


def __read_paz__(path):
    '''
    Reads a directory with paz files
    or a single file.
    Limitations:
    (1) directory must contain *only* paz files
    (2) paz file can optionally have ".pz" or ".paz" suffixes
    (3) paz file name (without prefix and suffix) *has* to have
        the trace_id (NET.STA.LOC.CHAN) of the corresponding trace
        in the last part of his name (e.g. 20110208_1600.NOW.IV.CRAC.00.EHZ.paz)
    '''
    if path is None:
        return None

    from obspy.sac.sacio import attach_paz
    from obspy.core import Trace
    logging.info('Reading PAZ...')
    paz = dict()
    if os.path.isdir(path):
        listing = os.listdir(path)
        #check if files have a common prefix: we will strip it later
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
    logging.info('Reading PAZ: done')
    return paz


def __parse_hypocenter__(hypo_file):
    hypo = AttribDict()
    hypo.latitude = None
    hypo.longitude = None
    hypo.depth = None
    hypo.origin_time = None
    hypo.evid = None

    if hypo_file is None:
        return None

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
        except IOError, message:
            logging.error(message)
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
        hypo.evid = evid.replace('.phs','').replace('.h','').replace('.hyp','')

    else: #FIXME: put a condition here!
        ev = hypo_file #FIXME: improve this!
        hypo.latitude = ev.latitude
        hypo.longitude = ev.longitude
        hypo.depth = ev.depth
        hypo.origin_time = ev.utcdate
        hypo.evid = ev.event_id

    return hypo


def __is_hypo_format(fp):
    for line in fp.readlines():
        # remove newline
        line = line.replace('\n','')
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
            fp.seek(0) #rewind file
            return True
        else:
            fp.seek(0) #rewind file
            return False


#TODO: def __is_NLL_format(fp):


def __parse_picks__(config):
    pick_file = config.options.pick_file
    if pick_file is None:
        return None

    try:
        with open(pick_file) as fp:
            picks = []
            if __is_hypo_format(fp):
                # Corinth phase file format is hypo
                for line in fp.readlines():
                    # remove newline
                    line = line.replace('\n','')
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
                    pick.station = __correct_station_name__(pick.station, config.traceids)
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
                    except:
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
                    # pick2.time has the same date, hour and minutes than pick.time
                    # We therefore make first a copy of pick.time...
                    pick2.time = UTCDateTime(pick.time)
                    # ...then set seconds and miscorseconds to 0...
                    pick2.time.second = 0
                    pick2.time.microsecond = 0
                    # ...and finally add stime
                    pick2.time += float(stime)

                    picks.append(pick2)
            else:
                raise IOError('%s: Not a phase file' % pick_file)
    except IOError, message:
        logging.error(message)
        ssp_exit(1)
    return picks


def __build_filelist__(path, filelist, tmpdir):
    if os.path.isdir(path):
        listing = os.listdir(path)
        for filename in listing:
            fullpath = os.path.join(path, filename)
            __build_filelist__(fullpath, filelist, tmpdir)
    else:
        try:
            open(path)
        except IOError, message:
            logging.error(message)
            return
        if tarfile.is_tarfile(path) and tmpdir!=None:
            tar = tarfile.open(path)
            tar.extractall(path=tmpdir)
            tar.close()
        else:
            filelist.append(path)
# -----------------------------------------------------------------------------


# PATH DISCOVERY --------------------------------------------------------------
# We try to guess the path of the hypo and pick file from the data dir
# This applies (for the moment) only to the Corinth format
def __set_hypo_file_path__(config):
    if config.options.hypo_file != None:
        return
    # try with the basename of the datadir
    if os.path.isdir(config.args[0]):
        hypo_file = config.args[0] + '.phs.h'
        try:
            open(hypo_file)
            config.options.hypo_file = hypo_file
        except:
            pass
    return


def __set_pick_file_path__(config):
    if config.options.pick_file != None:
        return
    # try with the basename of the datadir
    if os.path.isdir(config.args[0]):
        pick_file = config.args[0] + '.phs'
        try:
            open(pick_file)
            config.options.pick_file = pick_file
        except:
            pass
    return
# -----------------------------------------------------------------------------


# Public interface:
def read_traces(config):
    '''
    Read traces, store waveforms and metadata.
    '''
    # read dataless
    dataless = __read_dataless__(config.dataless)
    # read PAZ (normally this is an alternative to dataless)
    paz = __read_paz__(config.paz)

    # parse hypocenter file
    __set_hypo_file_path__(config)
    hypo = __parse_hypocenter__(config.options.hypo_file)
    # parse pick file
    __set_pick_file_path__(config)
    picks = __parse_picks__(config)

    # finally, read traces
    # traces can be defined in a pickle catalog (Antilles format)...
    if config.pickle_catalog:
        sys.path.append(config.pickle_classpath)
        with open(config.pickle_catalog, 'rb') as fp:
            catalog = pickle.load(fp)
        event = [ ev for ev in catalog if ev.event_id == config.options.evid ][0]
        hypo = __parse_hypocenter__(event)
        st = Stream()
        for trace in event.traces:
            if config.options.station is not None:
                if not trace.station == config.options.station:
                    continue
            try:
                tmpst = read(trace.trace_file, fsize=False)
            except Exception, error:
                print error
                continue
            for trace in tmpst.traces:
                st.append(trace)
                trace.stats.format = config.trace_format
                __add_paz_and_coords__(trace, dataless, paz)
                __add_hypocenter__(trace, hypo)
                __add_picks__(trace, picks) #FIXME: actually add picks!
                #__add_instrtype__(trace)
                trace.stats.instrtype = 'acc' #FIXME
        #ssp_exit()
    # ...or in standard files (CRL and ISNet)
    else:
        logging.info('Reading traces...')
        # phase 1: build a file list
        # ph 1.1: create a temporary dir and run '_build_filelist()'
        #         to move files to it and extract all tar archives
        tmpdir = tempfile.mkdtemp()
        filelist = []
        for arg in config.args:
            __build_filelist__(arg, filelist, tmpdir)
        # ph 1.2: rerun '_build_filelist()' in tmpdir to add to the
        #         filelist all the extraceted files
        listing = os.listdir(tmpdir)
        for filename in listing:
            fullpath = os.path.join(tmpdir, filename)
            __build_filelist__(fullpath, filelist, None)

        # phase 2: build a stream object from the file list
        st = Stream()
        for filename in filelist:
            try:
                tmpst = read(filename, fsize=False)
            except:
                logging.error('%s: Unable to read file as a trace: skipping' % filename)
                continue
            for trace in tmpst.traces:
                st.append(trace)
                trace.stats.format = config.trace_format
                __correct_traceid__(trace, config.traceids)
                __add_paz_and_coords__(trace, dataless, paz)
                __add_instrtype__(trace)
                __add_hypocenter__(trace, hypo)
                __add_picks__(trace, picks)

        shutil.rmtree(tmpdir)

    logging.info('Reading traces: done')
    if len(st.traces) == 0:
        logging.info('No trace loaded')
        ssp_exit()

    __complete_picks__(st)
    st.sort()
    #ssp_exit()
    return st
