# -*- coding: utf-8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Augment traces with station and event metadata.

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
import re
import logging
import contextlib
from obspy.core.util import AttribDict
from ..setup import config
from .station_metadata import PAZ
from .sac_header import (
    compute_sensitivity_from_SAC,
    get_instrument_from_SAC, get_station_coordinates_from_SAC,
    get_picks_from_SAC)
logger = logging.getLogger(__name__.rsplit('.', maxsplit=1)[-1])


def _skip_traces_from_config(traceid):
    """
    Skip traces with unknown channel orientation or ignored from config file.

    :param traceid: Trace ID.
    :type traceid: str

    :raises: RuntimeError if traceid is ignored from config file.
    """
    network, station, location, channel = traceid.split('.')
    orientation_codes = config.vertical_channel_codes +\
        config.horizontal_channel_codes_1 +\
        config.horizontal_channel_codes_2
    orientation = channel[-1]
    if orientation not in orientation_codes:
        raise RuntimeError(
            f'{traceid}: Unknown channel orientation: '
            f'"{orientation}": skipping trace'
        )
    # build a list of all possible ids, from station only
    # to full net.sta.loc.chan
    ss = [
        station,
        '.'.join((network, station)),
        '.'.join((network, station, location)),
        '.'.join((network, station, location, channel)),
    ]
    if config.use_traceids is not None:
        # - combine all use_traceids in a single regex
        # - escape the dots, otherwise they are interpreted as any character
        # - add a dot before the first asterisk, to avoid a pattern error
        combined = (
            "(" + ")|(".join(config.use_traceids) + ")"
        ).replace('.', r'\.').replace('(*', '(.*')
        if not any(re.match(combined, s) for s in ss):
            raise RuntimeError(f'{traceid}: ignored from config file')
    if config.ignore_traceids is not None:
        # - combine all ignore_traceids in a single regex
        # - escape the dots, otherwise they are interpreted as any character
        # - add a dot before the first asterisk, to avoid a pattern error
        combined = (
            "(" + ")|(".join(config.ignore_traceids) + ")"
        ).replace('.', r'\.').replace('(*', '(.*')
        if any(re.match(combined, s) for s in ss):
            raise RuntimeError(f'{traceid}: ignored from config file')


def _select_components(st):
    """
    Select requested components from stream

    :param st: ObsPy Stream object
    :type st: :class:`obspy.core.stream.Stream`

    :return: ObsPy Stream object
    :rtype: :class:`obspy.core.stream.Stream`
    """
    traces_to_keep = []
    for trace in st:
        try:
            _skip_traces_from_config(trace.id)
        except RuntimeError as e:
            logger.warning(str(e))
            continue
        # TODO: should we also filter by station here?
        # only use the station specified by the command line option
        # "--station", if any
        #if (config.options.station is not None and
        #        trace.stats.station != config.options.station):
        #    continue
        traces_to_keep.append(trace)
    # in-place update of st
    st.traces[:] = traces_to_keep[:]


def _correct_traceid(trace):
    """
    Correct traceid from config.TRACEID_MAP, if available.

    :param trace: ObsPy Trace object
    :type trace: :class:`obspy.core.trace.Trace`
    """
    if config.TRACEID_MAP is None:
        return
    with contextlib.suppress(KeyError):
        traceid = config.TRACEID_MAP[trace.get_id()]
        net, sta, loc, chan = traceid.split('.')
        trace.stats.network = net
        trace.stats.station = sta
        trace.stats.location = loc
        trace.stats.channel = chan


def _add_instrtype(trace):
    """
    Add instrtype to trace.

    :param trace: ObsPy Trace object
    :type trace: :class:`obspy.core.trace.Trace`
    """
    instrtype = None
    band_code = None
    instr_code = None
    trace.stats.instrtype = None
    # First, try to get the instrtype from channel name
    chan = trace.stats.channel
    if len(chan) > 2:
        band_code = chan[0]
        instr_code = chan[1]
    if instr_code in config.INSTR_CODES_VEL:
        # SEED standard band codes from higher to lower sampling rate
        # https://ds.iris.edu/ds/nodes/dmc/data/formats/seed-channel-naming/
        if band_code in ['G', 'D', 'E', 'S']:
            instrtype = 'shortp'
        if band_code in ['F', 'C', 'H', 'B']:
            instrtype = 'broadb'
    if instr_code in config.INSTR_CODES_ACC:
        instrtype = 'acc'
    if instrtype is None:
        # Let's see if there is an instrument name in SAC header (ISNet format)
        # In this case, we define band and instrument codes a posteriori
        instrtype, band_code, instr_code = get_instrument_from_SAC(trace)
        orientation = trace.stats.channel[-1]
        trace.stats.channel = ''.join((band_code, instr_code, orientation))
    trace.stats.instrtype = instrtype
    trace.stats.info = f'{trace.id} {trace.stats.instrtype}'


def _add_inventory(trace, inventory):
    """
    Add inventory to trace.

    :param trace: ObsPy Trace object
    :type trace: :class:`obspy.core.trace.Trace`
    :param inventory: ObsPy Inventory object
    :type inventory: :class:`obspy.core.inventory.Inventory`
    """
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
        paz.sensitivity = compute_sensitivity_from_SAC(trace)
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
    """
    Check if instrument type is consistent with units in inventory.

    :param trace: ObsPy Trace object
    :type trace: :class:`obspy.core.trace.Trace`
    """
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
    """
    Add coordinates to trace.

    :param trace: ObsPy Trace object
    :type trace: :class:`obspy.core.trace.Trace`
    """
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


def _add_event(trace, ssp_event):
    """
    Add ssp_event object to trace.

    :param trace: ObsPy Trace object
    :type trace: :class:`obspy.core.trace.Trace`
    :param ssp_event: SSPEvent object
    :type ssp_event: :class:`sourcespec.ssp_event.SSPEvent`
    """
    trace.stats.event = ssp_event


def _add_picks(trace, picks=None):
    """
    Add picks to trace.

    :param trace: ObsPy Trace object
    :type trace: :class:`obspy.core.trace.Trace`
    :param picks: list of picks (default: None)
    :type picks: list of :class:`sourcespec.ssp_event.Pick`
    """
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
    """
    Add component-specific picks to all components in a stream.

    :param st: ObsPy Stream object
    :type st: :class:`obspy.core.stream.Stream`
    """
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


def _augment_trace(trace, inventory, ssp_event, picks):
    """
    Augment trace with station and event metadata.

    :param trace: ObsPy Trace object
    :type trace: :class:`obspy.core.trace.Trace`
    :param inventory: ObsPy Inventory object
    :type inventory: :class:`obspy.core.inventory.Inventory`
    :param ssp_event: SSPEvent object
    :type ssp_event: :class:`sourcespec.ssp_event.SSPEvent`
    :param picks: list of picks
    :type picks: list of :class:`sourcespec.ssp_event.Pick`
    """
    _add_instrtype(trace)
    _add_inventory(trace, inventory)
    _check_instrtype(trace)
    _add_coords(trace)
    _add_event(trace, ssp_event)
    _add_picks(trace, picks)
    trace.stats.ignore = False


def augment_traces(st, inventory, ssp_event, picks):
    """
    Add all required information to trace headers.

    Only the traces that satisfy the conditions in the config file are kept.
    Problematic traces are also removed.

    :param st: Traces to be augmented
    :type st: :class:`obspy.core.stream.Stream`
    :param inventory: Station metadata. If it is None or an empty Inventory
        object, the code will try to read the station metadata from the
        trace headers (only SAC format is supported).
    :type inventory: :class:`obspy.core.inventory.Inventory`
    :param ssp_event: Event information
    :type ssp_event: :class:`sourcespec.ssp_event.SSPEvent`
    :param picks: list of picks
    :type picks: list of :class:`sourcespec.ssp_event.Pick`
    """
    # First, select the components based on the config options
    _select_components(st)
    # Then, augment the traces and remove the problematic ones
    traces_to_keep = []
    for trace in st:
        _correct_traceid(trace)
        try:
            _augment_trace(trace, inventory, ssp_event, picks)
        except Exception as err:
            for line in str(err).splitlines():
                logger.warning(line)
            continue
        traces_to_keep.append(trace)
    # in-place update of st
    st.traces[:] = traces_to_keep[:]
    _complete_picks(st)
