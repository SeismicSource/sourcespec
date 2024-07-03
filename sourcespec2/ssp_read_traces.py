# -*- coding: utf-8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Read traces in multiple formats of data and metadata.

:copyright:
    2012 Claudio Satriano <satriano@ipgp.fr>

    2013-2014 Claudio Satriano <satriano@ipgp.fr>,
              Emanuela Matrullo <matrullo@geologie.ens.fr>

    2015-2025 Claudio Satriano <satriano@ipgp.fr>,
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
import zipfile
import tempfile
import contextlib
from obspy import read
from obspy.core import Stream
from obspy.core.util import AttribDict
from .ssp_setup import (
    ssp_exit, INSTR_CODES_VEL, INSTR_CODES_ACC, TRACEID_MAP)
from .ssp_util import MediumProperties
from .ssp_read_station_metadata import (
    read_station_metadata, PAZ)
from .ssp_read_event_metadata import (
    parse_qml, parse_hypo_file, parse_hypo71_picks)
from .ssp_read_sac_header import (
    compute_sensitivity_from_SAC,
    get_instrument_from_SAC, get_station_coordinates_from_SAC,
    get_event_from_SAC, get_picks_from_SAC)
logger = logging.getLogger(__name__.rsplit('.', maxsplit=1)[-1])


# TRACE MANIPULATION ----------------------------------------------------------
def _correct_traceid(trace):
    if TRACEID_MAP is None:
        return
    with contextlib.suppress(KeyError):
        traceid = TRACEID_MAP[trace.get_id()]
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
    if instr_code in INSTR_CODES_VEL:
        # SEED standard band codes from higher to lower sampling rate
        # https://ds.iris.edu/ds/nodes/dmc/data/formats/seed-channel-naming/
        if band_code in ['G', 'D', 'E', 'S']:
            instrtype = 'shortp'
        if band_code in ['F', 'C', 'H', 'B']:
            instrtype = 'broadb'
    if instr_code in INSTR_CODES_ACC:
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
        # Inventory.get_coordinates() raises a generic Exception
        # if coordinates are not found
        coords = trace.stats.inventory.get_coordinates(
                    trace.id, trace.stats.starttime)
    if coords is not None:
        # Build an AttribDict and make sure that coordinates are floats
        coords = AttribDict({
            key: float(value) for key, value in coords.items()})
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
_add_coords.skipped = []  # noqa


def _add_event(trace, ssp_event=None):
    """Add ssp_event object to trace."""
    if ssp_event is None:
        # Try to get hypocenter information from the SAC header
        try:
            ssp_event = get_event_from_SAC(trace)
        except Exception:
            return
    trace.stats.event = ssp_event


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
    medium_properties = MediumProperties(
        hypo.longitude, hypo.latitude, hypo.depth.value_in_km, config)
    hypo.vp = medium_properties.get(mproperty='vp', where='source')
    hypo.vs = medium_properties.get(mproperty='vs', where='source')
    hypo.rho = medium_properties.get(mproperty='rho', where='source')
    depth_string = medium_properties.to_string(
        'source depth', hypo.depth.value_in_km)
    vp_string = medium_properties.to_string('vp_source', hypo.vp)
    vs_string = medium_properties.to_string('vs_source', hypo.vs)
    rho_string = medium_properties.to_string('rho_source', hypo.rho)
    logger.info(f'{depth_string}, {vp_string}, {vs_string}, {rho_string}')


def _build_filelist(path, filelist, tmpdir):
    if os.path.isdir(path):
        listing = os.listdir(path)
        for filename in listing:
            fullpath = os.path.join(path, filename)
            _build_filelist(fullpath, filelist, tmpdir)
    else:
        try:
            # pylint: disable=unspecified-encoding consider-using-with
            open(path)
        except IOError as err:
            logger.error(err)
            return
        if tarfile.is_tarfile(path) and tmpdir is not None:
            with tarfile.open(path) as tar:
                try:
                    tar.extractall(path=tmpdir)
                except Exception as msg:
                    logger.warning(
                        f'{path}: Unable to fully extract tar archive: {msg}')
        elif zipfile.is_zipfile(path) and tmpdir is not None:
            with zipfile.ZipFile(path) as zipf:
                try:
                    zipf.extractall(path=tmpdir)
                except Exception as msg:
                    logger.warning(
                        f'{path}: Unable to fully extract zip archive: {msg}')
        else:
            filelist.append(path)


def _read_trace_files(config, inventory, ssp_event, picks):
    """
    Read trace files from a given path. Complete trace metadata and
    return a stream object.
    """
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
                _add_event(trace, ssp_event)
                _add_picks(trace, picks)
            except Exception as err:
                for line in str(err).splitlines():
                    logger.warning(line)
                continue
            st.append(trace)
    shutil.rmtree(tmpdir)
    return st
# -----------------------------------------------------------------------------


def _log_event_info(ssp_event):
    for line in str(ssp_event).splitlines():
        logger.info(line)
    logger.info('---------------------------------------------------')


# Public interface:
def read_traces(config):
    """Read traces, store waveforms and metadata."""
    # read station metadata into an ObsPy ``Inventory`` object
    inventory = read_station_metadata(config.station_metadata)

    picks = []
    ssp_event = None
    # parse hypocenter file
    if config.options.hypo_file is not None:
        ssp_event, picks, file_format = parse_hypo_file(
            config.options.hypo_file, config.options.evid)
        config.hypo_file_format = file_format
    # parse pick file
    if config.options.pick_file is not None:
        picks = parse_hypo71_picks(config)
    # parse QML file
    if config.options.qml_file is not None:
        ssp_event, picks = parse_qml(config)
    if ssp_event is not None:
        _log_event_info(ssp_event)

    # finally, read trace files
    logger.info('Reading traces...')
    st = _read_trace_files(config, inventory, ssp_event, picks)
    logger.info('Reading traces: done')
    logger.info('---------------------------------------------------')
    if len(st) == 0:
        logger.error('No trace loaded')
        ssp_exit(1)
    _complete_picks(st)

    # if ssp_event is still None, get it from first trace
    if ssp_event is None:
        try:
            ssp_event = st[0].stats.event
            _log_event_info(ssp_event)
        except AttributeError:
            logger.error('No hypocenter information found.')
            sys.stderr.write(
                '\n'
                'Use "-q" or "-H" options to provide hypocenter information\n'
                'or add hypocenter information to the SAC file header\n'
                '(if you use the SAC format).\n'
            )
            ssp_exit(1)
    # add velocity info to hypocenter
    try:
        _hypo_vel(ssp_event.hypocenter, config)
    except Exception as e:
        logger.error(
            f'Unable to compute velocity at hypocenter: {e}\n')
        ssp_exit(1)
    if config.options.evname is not None:
        # add evname from command line, if any, overriding the one in ssp_event
        ssp_event.name = config.options.evname
    else:
        # add evname from ssp_event, if any, to config file
        config.options.evname = ssp_event.name
    # add event to config file
    config.event = ssp_event

    st.sort()
    return st
