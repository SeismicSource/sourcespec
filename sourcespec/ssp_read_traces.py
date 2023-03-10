# -*- coding: utf-8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Read traces in multiple formats of data and metadata.

:copyright:
    2012 Claudio Satriano <satriano@ipgp.fr>

    2013-2014 Claudio Satriano <satriano@ipgp.fr>,
              Emanuela Matrullo <matrullo@geologie.ens.fr>

    2015-2023 Claudio Satriano <satriano@ipgp.fr>,
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
import contextlib
from obspy import read
from obspy.core import Stream
from obspy.core.util import AttribDict
from obspy.io.sac import attach_paz
from sourcespec.ssp_setup import (
    ssp_exit, instr_codes_vel, instr_codes_acc, traceid_map)
from sourcespec.ssp_util import get_vel
from sourcespec.ssp_read_station_metadata import (
    read_station_metadata, PAZ)
from sourcespec.ssp_read_event_metadata import (
    parse_qml, parse_hypo_file, parse_hypo71_picks, Pick)
logger = logging.getLogger(__name__.split('.')[-1])


# TRACE MANIPULATION ----------------------------------------------------------
def _correct_traceid(trace):
    if traceid_map is None:
        return
    with contextlib.suppress(KeyError):
        traceid = traceid_map[trace.get_id()]
        net, sta, loc, chan = traceid.split('.')
        trace.stats.network = net
        trace.stats.station = sta
        trace.stats.location = loc
        trace.stats.channel = chan


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
        except Exception as e:
            raise TypeError(f'{trace.id}: trace must be in SAC format') from e
    try:
        sensitivity = eval(inp, {}, namespace)
    except NameError as msg:
        hdr_field = str(msg).split()[1]
        logger.error(f'SAC header field {hdr_field} does not exist')
        ssp_exit(1)
    return sensitivity


def _add_inventory(trace, inventory, config):
    """Add inventory to trace."""
    net, sta, loc, chan = trace.id.split('.')
    inv = (
        inventory.select(
            network=net, station=sta, location=loc, channel=chan)
        or inventory.select(
            network='XX', station='GENERIC', location='XX', channel='XXX')
    )
    if 'XX.GENERIC.XX.XXX' in inv.get_contents()['channels']:
        inv = inv.copy()
        inv.networks[0].code = net
        inv.networks[0].stations[0].code = sta
        inv.networks[0].stations[0].channels[0].code = chan
        inv.networks[0].stations[0].channels[0].location_code = loc
    # If a "sensitivity" config option is provided, override the Inventory
    # object with a new one constructed from the sensitivity value
    if config.sensitivity is not None:
        # save coordinates from the inventory, if available
        coords = None
        with contextlib.suppress(Exception):
            coords = inv.get_coordinates(trace.id, trace.stats.starttime)
        paz = PAZ()
        paz.seedID = trace.id
        paz.sensitivity = _compute_sensitivity(trace, config)
        paz.poles = []
        paz.zeros = []
        if inv:
            logger.warning(
                f'Overriding response for {trace.id} with constant '
                f'sensitivity {paz.sensitivity}')
        inv = paz.to_inventory()
        # restore coordinates, if available
        if coords is not None:
            chan = inv.networks[0].stations[0].channels[0]
            chan.latitude = coords['latitude']
            chan.longitude = coords['longitude']
            chan.elevation = coords['elevation']
            chan.depth = coords['local_depth']
    trace.stats.inventory = inv


def _add_coords(trace):
    """Add coordinates to trace."""
    # If we already know that traceid is skipped, raise a silent exception
    if trace.id in _add_coords.skipped:
        raise RuntimeError()
    coords = None
    with contextlib.suppress(Exception):
        coords = AttribDict(
            trace.stats.inventory.get_coordinates(
                trace.id, trace.stats.starttime))
        coords = (
            None if (
                coords.longitude == 0 and coords.latitude == 0 and
                coords.local_depth == 123456 and coords.elevation == 123456)
            else coords)
    # If we still don't have trace coordinates,
    # we try to get them from SAC header
    if coords is None:
        try:
            stla = trace.stats.sac.stla
            stlo = trace.stats.sac.stlo
            stel = trace.stats.sac.get('stel', 0.)
            coords = AttribDict(latitude=stla, longitude=stlo, elevation=stel)
            logger.info(
                f'{trace.id}: station coordinates read from SAC header')
        except Exception as e:
            _add_coords.skipped.append(trace.id)
            raise RuntimeError(
                f'{trace.id}: could not find coords for trace: skipping trace'
            ) from e
    if coords.latitude == coords.longitude == 0:
        logger.warning(
            f'{trace.id}: trace has latitude and longitude equal to zero!')
    # elevation is in meters in StationXML or SAC header
    coords.elevation /= 1e3
    trace.stats.coords = coords
# list to keep track of skipped traces
_add_coords.skipped = []  #noqa


def _add_paz(trace):
    """Add PAZ to trace."""
    traceid = trace.id
    time = trace.stats.starttime
    try:
        with warnings.catch_warnings(record=True) as warns:
            # get_sacpz() can issue warnings on more than one PAZ found,
            # so let's catch those warnings and log them properly
            # warnings.filterwarnings('ignore', message='Found more than')
            # warnings.filterwarnings('ignore', message='More than')
            inventory = trace.stats.inventory
            sacpz = inventory.get_response(traceid, time).get_sacpz()
            for w in warns:
                msg = str(w.message)
                logger.warning(f'{traceid}: {msg} Time: {time}')
        attach_paz(trace, io.StringIO(sacpz))
    except Exception as msg:
        logger.warning(f'{traceid}: {msg} Time: {time}')


def _add_instrtype(trace):
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
    trace.stats.info = f'{trace.id} {trace.stats.instrtype}'


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
def _hypo_vel(hypo, config):
    hypo.vp = get_vel(hypo.longitude, hypo.latitude, hypo.depth, 'P', config)
    hypo.vs = get_vel(hypo.longitude, hypo.latitude, hypo.depth, 'S', config)
    msg = 'Vp_hypo: {:.2f} km/s, Vs_hypo: {:.2f} km/s'.format(
        hypo.vp, hypo.vs)
    logger.info(msg)


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
    # read station metadata into an ObsPy ``Inventory`` object
    inventory = read_station_metadata(config.station_metadata)

    hypo = picks = None
    # parse hypocenter file
    if config.options.hypo_file is not None:
        hypo, picks = parse_hypo_file(config.options.hypo_file)
    # parse pick file
    if config.options.pick_file is not None:
        picks = parse_hypo71_picks(config)
    # parse QML file
    if config.options.qml_file is not None:
        hypo, picks = parse_qml(config.options.qml_file, config.options.evid)

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
            # only use the station specified by the command line option
            # "--station", if any
            if (config.options.station is not None and
                    trace.stats.station != config.options.station):
                continue
            _correct_traceid(trace)
            try:
                _add_inventory(trace, inventory, config)
                _add_coords(trace)
                _add_paz(trace)
                _add_instrtype(trace)
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
