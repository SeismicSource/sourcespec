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
import logging
import shutil
import tarfile
import tempfile
import contextlib
from obspy import read
from obspy.core import Stream
from obspy.core.util import AttribDict
from sourcespec.ssp_setup import (
    ssp_exit, instr_codes_vel, instr_codes_acc, traceid_map)
from sourcespec.ssp_util import get_vel
from sourcespec.ssp_read_station_metadata import (
    read_station_metadata, PAZ)
from sourcespec.ssp_read_event_metadata import (
    parse_qml, parse_hypo_file, parse_hypo71_picks)
from sourcespec.ssp_read_sac_header import (
    compute_sensitivity_from_SAC,
    get_instrument_from_SAC, get_station_coordinates_from_SAC,
    get_hypocenter_from_SAC, get_picks_from_SAC)
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


def _add_instrtype(trace):
    """Add instrtype to trace."""
    instrtype = None
    band_code = None
    instr_code = None
    trace.stats.instrtype = None
    # First, try to get the instrtype from channel name
    chan = trace.stats.channel
    if len(chan) > 2:
        band_code = chan[0]
        instr_code = chan[1]
    if instr_code in instr_codes_vel:
        # SEED standard band codes from higher to lower sampling rate
        # https://ds.iris.edu/ds/nodes/dmc/data/formats/seed-channel-naming/
        if band_code in ['G', 'D', 'E', 'S']:
            instrtype = 'shortp'
        if band_code in ['F', 'C', 'H', 'B']:
            instrtype = 'broadb'
    if instr_code in instr_codes_acc:
        instrtype = 'acc'
    if instrtype is None:
        # Let's see if there is an instrument name in SAC header (ISNet format)
        # In this case, we define band and instrument codes a posteriori
        instrtype, band_code, instr_code = get_instrument_from_SAC(trace)
        orientation = trace.stats.channel[-1]
        trace.stats.channel = ''.join((band_code, instr_code, orientation))
    trace.stats.instrtype = instrtype
    trace.stats.info = f'{trace.id} {trace.stats.instrtype}'


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
        paz.sensitivity = compute_sensitivity_from_SAC(trace, config)
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


def _check_instrtype(trace):
    """Check if instrument type is consistent with units in inventory."""
    inv = trace.stats.inventory
    if not inv:
        raise RuntimeError(
            f'{trace.id}: cannot get instrtype from inventory: '
            'inventory is empty: skipping trace')
    instrtype = trace.stats.instrtype
    new_instrtype = None
    try:
        units = inv.get_response(trace.id, trace.stats.starttime).\
            instrument_sensitivity.input_units
    except Exception as e:
        # inventory attached to trace has only one channel
        chan = inv[0][0][0]
        start_date = chan.start_date
        end_date = chan.end_date
        raise RuntimeError(
            f'{trace.id}: cannot get units from inventory.\n'
            f'> {e.__class__.__name__}: {e}\n'
            f'> Channel start/end date: {start_date} {end_date}\n'
            f'> Trace start time: {trace.stats.starttime}\n'
            '> Skipping trace'
        ) from e
    trace.stats.units = units
    if units.lower() == 'm' and trace.stats.instrtype != 'disp':
        new_instrtype = 'disp'
    if units.lower() == 'm/s' and instrtype not in ['shortp', 'broadb']:
        new_instrtype = 'broadb'
    if units.lower() == 'm/s**2' and instrtype != 'acc':
        new_instrtype = 'acc'
    if new_instrtype is not None:
        logger.warning(
            f'{trace.id}: instrument response units are "{units}" but '
            f'instrument type is "{instrtype}". Changing instrument type '
            f'to "{new_instrtype}"')
        trace.stats.instrtype = new_instrtype


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
    if coords is None:
        # If we still don't have trace coordinates,
        # we try to get them from SAC header
        coords = get_station_coordinates_from_SAC(trace)
    if coords is None:
        # Give up!
        _add_coords.skipped.append(trace.id)
        raise RuntimeError(
            f'{trace.id}: could not find coords for trace: skipping trace')
    if coords.latitude == coords.longitude == 0:
        logger.warning(
            f'{trace.id}: trace has latitude and longitude equal to zero!')
    # elevation is in meters in StationXML or SAC header
    coords.elevation /= 1e3
    trace.stats.coords = coords
# list to keep track of skipped traces
_add_coords.skipped = []  #noqa


def _add_hypocenter(trace, hypo=None):
    """Add hypocenter information to trace."""
    if hypo is None:
        # Try to get hypocenter information from the SAC header
        try:
            hypo = get_hypocenter_from_SAC(trace)
        except Exception:
            return
    trace.stats.hypo = hypo


def _add_picks(trace, picks=None):
    """Add picks to trace."""
    if picks is None:
        picks = []
    trace_picks = []
    with contextlib.suppress(Exception):
        trace_picks = get_picks_from_SAC(trace)
    for pick in picks:
        if pick.station == trace.stats.station:
            trace_picks.append(pick)
    trace.stats.picks = trace_picks
    # Create empty dicts for arrivals, travel_times and takeoff angles.
    # They will be used later.
    trace.stats.arrivals = {}
    trace.stats.travel_times = {}
    trace.stats.takeoff_angles = {}


def _complete_picks(st):
    """Add component-specific picks to all components."""
    for station in {tr.stats.station for tr in st}:
        st_sel = st.select(station=station)
        # 'code' is band+instrument code
        for code in {tr.stats.channel[:-1] for tr in st_sel}:
            st_sel2 = st_sel.select(channel=f'{code}?')
            all_picks = [p for tr in st_sel2 for p in tr.stats.picks]
            all_P_picks = [p for p in all_picks if p.phase == 'P']
            all_S_picks = [p for p in all_picks if p.phase == 'S']
            # Select default P and S picks as the first in list (or empty list)
            default_P_pick = all_P_picks[:1]
            default_S_pick = all_S_picks[:1]
            for tr in st_sel2:
                # Attribute default picks to components without picks
                if not [p for p in tr.stats.picks if p.phase == 'P']:
                    tr.stats.picks += default_P_pick
                if not [p for p in tr.stats.picks if p.phase == 'S']:
                    tr.stats.picks += default_S_pick
# -----------------------------------------------------------------------------


# FILE PARSING ----------------------------------------------------------------
def _hypo_vel(hypo, config):
    hypo.vp = get_vel(hypo.longitude, hypo.latitude, hypo.depth, 'P', config)
    hypo.vs = get_vel(hypo.longitude, hypo.latitude, hypo.depth, 'S', config)
    logger.info(f'Vp_hypo: {hypo.vp:.2f} km/s, Vs_hypo: {hypo.vs:.2f} km/s')


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


def _read_trace_files(config, inventory, hypo, picks, green=False):
    """
    Read trace files from a given path. Complete trace metadata and
    return a stream object.
    """
    # phase 1: build a file list
    # ph 1.1: create a temporary dir and run '_build_filelist()'
    #         to move files to it and extract all tar archives
    tmpdir = tempfile.mkdtemp()
    filelist = []
    config_trace_path =\
        config.options.green_path if green else config.options.trace_path
    for trace_path in config_trace_path:
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
            logger.warning(
                f'{filename}: Unable to read file as a trace: skipping')
            continue
        for trace in tmpst.traces:
            orientation = trace.stats.channel[-1]
            if orientation not in orientation_codes:
                logger.warning(
                    f'{trace.id}: Unknown channel orientation: '
                    f'"{orientation}": skipping trace'
                )
                continue
            # only use the station specified by the command line option
            # "--station", if any
            if (config.options.station is not None and
                    trace.stats.station != config.options.station):
                continue
            _correct_traceid(trace)
            try:
                _add_instrtype(trace)
                _add_inventory(trace, inventory, config)
                _check_instrtype(trace)
                _add_coords(trace)
                _add_hypocenter(trace, hypo)
                _add_picks(trace, picks)
            except Exception as err:
                for line in str(err).splitlines():
                    logger.warning(line)
                continue
            st.append(trace)
    shutil.rmtree(tmpdir)
    return st
# -----------------------------------------------------------------------------


# Public interface:
def read_traces(config):
    """Read traces, store waveforms and metadata."""
    # read station metadata into an ObsPy ``Inventory`` object
    inventory = read_station_metadata(config.station_metadata)

    picks, picksG = [], []
    hypo, hypoG = None, None
    # parse hypocenter file
    if config.options.hypo_file is not None:
        hypo, picks = parse_hypo_file(config.options.hypo_file)
    # parse pick file
    if config.options.pick_file is not None:
        picks = parse_hypo71_picks(config)
    # parse QML file
    if config.options.qml_file is not None:
        hypo, picks = parse_qml(config.options.qml_file, config.options.evid)

    if config.options.greenqml_file is not None:
        hypoG, picksG = parse_qml(
            config.options.greenqml_file, config.options.evid)

    # finally, read trace files
    logger.info('Reading traces...')
    st = _read_trace_files(config, inventory, hypo, picks)
    logger.info('Reading traces: done')
    if len(st) == 0:
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
    # add vp and vs to hypo
    _hypo_vel(hypo, config)
    # add hypo to config file
    config.hypo = hypo
    # same for hypoG, if available
    if hypoG is not None:
        _hypo_vel(hypoG, config)
        config.hypoG = hypoG
    else:
        config.hypoG = None

    # if green function traces are available, read traces
    if config.options.green_path is not None:
        logger.info('Green function traces available')
        # read green function traces and adds it to a stream
        st_green = _read_trace_files(
            config, inventory, hypoG, picksG, green=True)
        logger.info('Reading green function traces: done')
        if len(st_green) == 0:
            logger.info('No green function trace loaded')
            ssp_exit()
        _complete_picks(st_green)
        # add green function to the main stream object
        # Traces of the main event and of the green function
        # are identified by evid metadata in the stream
        st = st + st_green

    st.sort()
    return st
