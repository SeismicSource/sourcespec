# -*- coding: utf-8 -*-
# ssp_read_traces.py
#
# Read traces for source_spectrum
# All the functions whose name is between "__" are intended to be private
# (c) 2012 Claudio Satriano <satriano@ipgp.fr>
from __future__ import division
import sys
import os
import shutil
import tarfile
import tempfile
from datetime import datetime
from obspy.core import Stream, read, UTCDateTime
from obspy.core.util import AttribDict
from obspy.xseed import Parser
from obspy.xseed.utils import SEEDParserException

# TRACE MANIPULATION ----------------------------------------------------------
__correct_traceid_dict__={
	'CL.AGE  00..E' : 'CL.AGE.00.EHE',
	'CL.AGE  00..N' : 'CL.AGE.00.EHN',
	'CL.AGE  00..Z' : 'CL.AGE.00.EHZ',
	'CL.AIO  00..E' : 'CL.AIO.00.EHE',
	'CL.AIO  00..N' : 'CL.AIO.00.EHN',
	'CL.AIO  00..Z' : 'CL.AIO.00.EHZ',
	'CL.ALI  00..E' : 'CL.ALI.00.EHE',
	'CL.ALI  00..N' : 'CL.ALI.00.EHN',
	'CL.ALI  00..Z' : 'CL.ALI.00.EHZ',
	'CL.DIM  00..E' : 'CL.DIM.00.EHE',
	'CL.DIM  00..N' : 'CL.DIM.00.EHN',
	'CL.DIM  00..Z' : 'CL.DIM.00.EHZ',
	'CL.KOU  00..E' : 'CL.KOU.00.EHE',
	'CL.KOU  00..N' : 'CL.KOU.00.EHN',
	'CL.KOU  00..Z' : 'CL.KOU.00.EHZ',
	'CL.PAN  00..E' : 'CL.PAN.00.EHE',
	'CL.PAN  00..N' : 'CL.PAN.00.EHN',
	'CL.PAN  00..Z' : 'CL.PAN.00.EHZ',
	'CL.PSA  00..E' : 'CL.PSA.00.EHE',
	'CL.PSA  00..N' : 'CL.PSA.00.EHN',
	'CL.PSA  00..Z' : 'CL.PSA.00.EHZ',
	'CL.PYR  00..E' : 'CL.PYR.00.EHE',
	'CL.PYR  00..N' : 'CL.PYR.00.EHN',
	'CL.PYR  00..Z' : 'CL.PYR.00.EHZ',
	'CL.TEM  00..E' : 'CL.TEM.00.EHE',
	'CL.TEM  00..N' : 'CL.TEM.00.EHN',
	'CL.TEM  00..Z' : 'CL.TEM.00.EHZ'
}

def __correct_traceid__(trace):
	try:
		traceid = __correct_traceid_dict__[trace.getId()]
		net, sta, loc, chan = traceid.split('.')
		trace.stats.network = net
		trace.stats.station = sta
		trace.stats.location = loc
		trace.stats.channel = chan
	except KeyError:
		pass

def __add_paz_and_coords__(trace, dataless):
	trace.stats.paz = None
	trace.stats.coords = None
	# We first look into the dataless dictionary, if available
	if dataless != None:
		for sp in dataless.values():
			try:
				traceid = trace.getId()
				time = trace.stats.starttime
				paz = sp.getPAZ(traceid, time)
				coords = AttribDict(sp.getCoordinates(traceid, time))
				# elevation is in meters in the dataless
				coords.elevation /= 1000
				trace.stats.paz = paz
				trace.stats.coords = coords
			except SEEDParserException:
				pass
	# If we couldn't find any PAZ in the dataless dictionary,
	# we try to build the sensitivity from the
	# user2 and user3 header fields (ISNet format)
	if trace.stats.paz == None:
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
	# Same thing for the station coordinates
	if trace.stats.coords == None:
		try: 
			stla = trace.stats.sac.stla
			stlo = trace.stats.sac.stlo
			stel = trace.stats.sac.stel
			# elevation is in meters in SAC header
			stel /= 1000
		 	if stla==-12345 or stlo==-12345 or stel==-12345:
				raise AttributeError
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
	try: instr_code = chan[1]
	except IndexError: instr_code = None
	if instr_code == 'H': instrtype = 'vel'
	if instr_code == 'L': instrtype = 'vel'
	if instr_code == 'N': instrtype = 'acc'

	# If, not possible, let's see if there is an instrument
	# name in "kinst" (ISNet format)
	if instrtype == None:
		try:
			instr = trace.stats.sac.kinst
			if 'CMG-5T' in instr:
				instrtype = 'acc'
			if 'TRILLIUM' in instr or 'S13J' in instr:
				instrtype = 'vel'
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
			   or evdp == -12345 or tori == -12345\
			   or begin == -12345:
				raise AttributeError
			hypo = AttribDict()
			hypo.latitude = evla
			hypo.longitude = evlo
			hypo.depth = evdp
			hypo.origin_time = trace.stats.starttime + tori - begin
			hypo.evid = None
		except AttributeError:
			pass
	trace.stats.hypo = hypo

def __add_picks__(trace, picks):
	# TODO: try to get picks from SAC header
	if picks == None:
		trace.stats.picks = []
		return
	
	stat_picks=[]
	station = trace.stats.station
	for pick in picks:
		if pick.station == station:
			stat_picks.append(pick)
	trace.stats.picks = stat_picks
# -----------------------------------------------------------------------------



# FILE PARSING ----------------------------------------------------------------
def __read_dataless__(path):
	if path == None: return None

	sys.stdout.write('Reading dataless...')
	dataless=dict()
	if os.path.isdir(path):
		listing = os.listdir(path)
		for filename in listing:
			fullpath='%s/%s' % (path, filename)
			try:
				sp = Parser(fullpath)
				dataless[filename] = sp
			except IOError: continue
	sys.stdout.write(' done.\n')
	return dataless

def __parse_hypocenter__(hypo_file):
	hypo = AttribDict()
	hypo.latitude = None
	hypo.longitude = None
	hypo.depth = None
	hypo.origin_time = None
	hypo.evid = None

	if hypo_file == None: return None

	try: fp = open(hypo_file)
	except: return None

	# Corinth hypocenter file format:
	# TODO: check file format
	line = fp.readline()
	fp.close()
	timestr = line[0:17]
	dt = datetime.strptime(timestr, '%y%m%d %H %M%S.%f')
	hypo.origin_time = UTCDateTime(dt)

	lat = float(line[17:20])
	lat_deg = float(line[20:26])
	hypo.latitude = lat + lat_deg/60
	lon = float(line[26:30])
	lon_deg = float(line[30:36])
	hypo.longitude = lon + lon_deg/60
	hypo.depth = float(line[36:42])
	hypo.evid = line[72:82]

	return hypo

def __new_pick__():
	pick = AttribDict()
	pick.station  = None
	pick.flag     = None
	pick.phase    = None
	pick.polarity = None
	pick.quality  = None
	pick.time     = None
	return pick

def __parse_picks__(pick_file):
	if pick_file == None: return None

	try: fp = open(pick_file)
	except: return None

	picks = []

	# Corinth hypocenter file format:
	# TODO: check file format
	for line in fp.readlines():
		# remove newline
		line = line.replace('\n','')
		# skip separator and empty lines
		stripped_line = line.replace(' ','')
		if stripped_line == '10' or stripped_line == '': continue

		pick = __new_pick__()
		pick.station  = line[0:4]
		pick.flag     = line[4:5]
		pick.phase    = line[5:6]
		pick.polarity = line[6:7]
		pick.quality  = int(line[7:8])
		timestr       = line[9:24]
		dt = datetime.strptime(timestr, '%y%m%d%H%M%S.%f')
		pick.time = UTCDateTime(dt)

		picks.append(pick)

		try: stime = line[31:36]
		except: continue
		if stime.replace(' ','') == '': continue

		pick2 = __new_pick__()
		pick2.station  = pick.station
		pick2.flag     = line[36:37]
		pick2.phase    = line[37:38]
		pick2.polarity = line[38:39]
		pick2.quality  = int(line[39:40])
		pick2.time     = pick.time + float(stime)
		#print 'aa' + stime + 'aa'

		picks.append(pick2)

	fp.close()
	return picks


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
			sys.stderr.write('%s\n' % message)
			return
		if tarfile.is_tarfile(path) and tmpdir!=None:
			tar = tarfile.open(path)
			tar.extractall(path=tmpdir)
			tar.close()
		else:
			filelist.append(path)
# -----------------------------------------------------------------------------



# Public interface:
def read_traces(args, options):
	# read dataless	
	dataless = __read_dataless__(options.dataless)
	# parse hypocenter file
	hypo = __parse_hypocenter__(options.hypo_file)
	# parse pick file
	picks = __parse_picks__(options.pick_file)

	# finally, read traces
	sys.stdout.write('Reading traces...')
	# phase 1: build a file list
	# ph 1.1: create a temporary dir and run '_build_filelist()'
	#         to move files to it and extract all tar archives
	tmpdir = tempfile.mkdtemp()
	filelist = []
	for arg in args:
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
			tmpst = read(filename)
		except:
			sys.stderr.write('%s: Unable to read file as a trace: skipping\n' % filename)
			continue
		for trace in tmpst.traces:
			st.append(trace)
			__correct_traceid__(trace)
			__add_paz_and_coords__(trace, dataless)
			__add_instrtype__(trace)
			__add_hypocenter__(trace, hypo)
			__add_picks__(trace, picks)

	shutil.rmtree(tmpdir)
	sys.stdout.write(' done.\n')
	return st
