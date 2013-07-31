# -*- coding: utf-8 -*-
# ssp_read_traces.py
#
# Read traces for source_spec
# All the functions whose name is between "__" are intended to be private
# (c) 2012 Claudio Satriano <satriano@ipgp.fr>
# (c) 2013 Claudio Satriano <satriano@ipgp.fr>,
#          Emanuela Matrullo <matrullo@geologie.ens.fr>
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
    def __init__(self):
        self.station  = None
        self.flag     = None
        self.phase    = None
        self.polarity = None
        self.quality  = None
        self.time     = None


# TRACE MANIPULATION ----------------------------------------------------------
def __correct_traceid__(trace, traceids):
    if traceids == None:
        return
    try:
        trids = load_source('traceids', traceids)
        traceid = trids.__correct_traceid_dict__[trace.getId()]
        net, sta, loc, chan = traceid.split('.')
        trace.stats.network = net
        trace.stats.station = sta
        trace.stats.location = loc
        trace.stats.channel = chan
    except KeyError:
        pass

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
                coords.elevation /= 1000
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
    if trace.stats.paz == None:
        if paz_dict != None:
            try:
                paz = paz_dict[trace.id]
                trace.stats.paz = paz
            except KeyError:
                pass
    # If we're still out of luck,
    # we try to build the sensitivity from the
    # user2 and user3 header fields (ISNet format)
    if trace.stats.paz == None and trace.stats.format == 'ISNet':
        try:
            # instrument constants
            u2 = trace.stats.sac.user2
            u3 = trace.stats.sac.user3
            if u2==-12345 or u3==-12345: raise AttributeError
            paz = AttribDict()
            paz.sensitivity = u3/u2
            paz.poles = []
            paz.zeros = []
            paz.gain = 1
            trace.stats.paz = paz
        except AttributeError:
            pass
    # Still no paz? Antilles or IPOC format!
    if trace.stats.paz == None\
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
    if trace.stats.coords == None:
        try:
            stla = trace.stats.sac.stla
            stlo = trace.stats.sac.stlo
            stel = trace.stats.sac.stel
            if stla==-12345 or stlo==-12345:
                raise AttributeError
            if stel==-12345:
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

    # First, try to get the instrtype from channel name
    chan = trace.stats.channel
    try: band_code = chan[0]
    except IndexError: band_code = None
    try: instr_code = chan[1]
    except IndexError: instr_code = None
    if instr_code == 'H' or instr_code == 'L':
        if band_code == 'E':
            instrtype = 'shortp'
        if band_code == 'H':
            instrtype = 'broadb'
    if instr_code == 'N' or instr_code =='L':
        instrtype = 'acc'

    # If, not possible, let's see if there is an instrument
    # name in "kinst" (ISNet format)
    if instrtype == None:
        try:
            instr = trace.stats.sac.kinst
            if 'CMG-5T' in instr:
                instrtype = 'acc'
            if 'TRILLIUM' in instr:
                instrtype = 'broadb'
            if 'S13J' in instr:
                instrtype = 'shortp'
            if 'KS2000ED' in instr:
                instrtype = 'shortp'
        except AttributeError:
            pass
    trace.stats.instrtype = instrtype

def __add_hypocenter__(trace, hypo):
    if hypo == None:
        # Try to get hypocenter information from the SAC header
        try:
            evla = trace.stats.sac.evla
            evlo = trace.stats.sac.evlo
            evdp = trace.stats.sac.evdp
            tori = trace.stats.sac.o
            begin = trace.stats.sac.b
            if evla == -12345 or evlo == -12345\
               or evdp == -12345 or begin == -12345:
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
    stat_picks=[]
    station = trace.stats.station

    if picks == None:
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
            pick.station  = station
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
# -----------------------------------------------------------------------------



# FILE PARSING ----------------------------------------------------------------
def __read_dataless__(path):
    if path == None: return None

    logging.info('Reading dataless...')
    dataless = dict()
    if os.path.isdir(path):
        listing = os.listdir(path)
        for filename in listing:
            fullpath='%s/%s' % (path, filename)
            try:
                sp = Parser(fullpath)
                dataless[filename] = sp
            except IOError: continue
        #TODO: manage the case in which "path" is a file name
    logging.info('Reading dataless: done')
    return dataless

def __read_paz__(path):
    '''
    Reads a directory with paz files
    (TODO: read a single file)
    Limitations:
    (1) directory must contain *only* paz files
    (2) all the paz files must have the same prefix (if any)
    (3) paz file can optionally have th ".pz" suffix
    (4) paz file name (without prefix and suffix) *has* to be
        the trace_id (NET.STA.LOC.CHAN) of the corresponding trace
    '''
    if path == None: return None

    from obspy.sac.sacio import attach_paz
    from obspy.core import Trace
    logging.info('Reading PAZ...')
    paz = dict()
    if os.path.isdir(path):
        listing = os.listdir(path)
        #check if files have a common prefix: we will strip it later
        prefix = os.path.commonprefix(listing)
        for filename in listing:
            fullpath='%s/%s' % (path, filename)
            try:
                # This is a horrible hack!
                # Since attach_paz needs a trace,
                # we create a trace and then, later,
                # we just retrieve the paz object
                # from the trace ;)
                tr = Trace()
                attach_paz(tr, fullpath)
                bname = os.path.basename(filename)
                #strip .pz suffix, if there
                bname = re.sub('.pz$', '', bname)
                #and strip any common prefix
                #we assume that the string which remains
                #is the trace_id
                trace_id = re.sub('^' + prefix, '', bname)
                paz[trace_id] = tr.stats.paz.copy()
            except IOError:
                continue
    #TODO: manage the case in which "path" is a file name
    logging.info('Reading PAZ: done')
    return paz

def __parse_hypocenter__(hypo_file):
    hypo = AttribDict()
    hypo.latitude = None
    hypo.longitude = None
    hypo.depth = None
    hypo.origin_time = None
    hypo.evid = None

    if hypo_file == None: return None

    if isinstance(hypo_file, str):
        try: fp = open(hypo_file)
        except: return None

        # Corinth hypocenter file format:
        # TODO: check file format
        line = fp.readline()
        # Skip the first line if it contains characters in the first 10 digits:
        if any(c.isalpha() for c in line[0:10]):
            line = fp.readline()
        fp.close()
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
        stripped_line = line.replace(' ','')
        if stripped_line == '10' or stripped_line == '': continue
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

def __parse_picks__(pick_file):
    if pick_file == None: return None

    try: fp = open(pick_file)
    except: return None

    picks = []

    if __is_hypo_format(fp):
        # Corinth phase file format is hypo
        for line in fp.readlines():
            # remove newline
            line = line.replace('\n','')
            # skip separator and empty lines
            stripped_line = line.replace(' ','')
            if stripped_line == '10' or stripped_line == '':
                continue
            # Check if it is a pick line
            # 6th character should be alpha (phase name: P or S)
            # other character should be digits (date/time)
            if not (line[5].isalpha() and
                             line[9].isdigit() and
                             line[20].isdigit()):
                continue

            #pick = __new_pick__()
            pick = Pick()
            pick.station  = line[0:4]
            pick.flag     = line[4:5]
            pick.phase    = line[5:6]
            pick.polarity = line[6:7]
            try:
                pick.quality = int(line[7:8])
            except ValueError:
                # If we cannot read pick quality,
                # we give the pick the lowest quality
                pick.quality = 4
            timestr       = line[9:24]
            dt = datetime.strptime(timestr, '%y%m%d%H%M%S.%f')
            pick.time = UTCDateTime(dt)

            picks.append(pick)

            try: stime = line[31:36]
            except: continue
            if stime.replace(' ','') == '': continue

            #pick2 = __new_pick__()
            pick2 = Pick()
            pick2.station  = pick.station
            pick2.flag     = line[36:37]
            pick2.phase    = line[37:38]
            pick2.polarity = line[38:39]
            try:
                pick2.quality = int(line[39:40])
            except ValueError:
                # If we cannot read pick quality,
                # we give the pick the lowest quality
                pick2.quality = 4
            pick2.time     = pick.time + float(stime)

            picks.append(pick2)

            fp.close()
            return picks
        else:
            raise IOError('%s: Not a phase file' % pick_file)



def __build_filelist__(path, filelist, tmpdir):
    if os.path.isdir(path):
        listing = os.listdir(path)
        for filename in listing:
            fullpath='%s/%s' % (path, filename)
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
    # read dataless
    dataless = __read_dataless__(config.dataless)
    # read PAZ (normally this is an alternative to dataless)
    paz = __read_paz__(config.paz)

    # parse hypocenter file
    __set_hypo_file_path__(config)
    hypo = __parse_hypocenter__(config.options.hypo_file)
    # parse pick file
    __set_pick_file_path__(config)
    picks = __parse_picks__(config.options.pick_file)

    # finally, read traces
    # traces can be defined in a pickle catalog (Antilles format)...
    if config.pickle_catalog:
        sys.path.append(config.pickle_classpath)
        fp = open(config.pickle_catalog, 'rb')
        catalog = pickle.load(fp)
        fp.close()
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
            fullpath='%s/%s' % (tmpdir, filename)
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

    st.sort()
    return st
