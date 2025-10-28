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
import logging
import contextlib
from obspy.core.util import AttribDict
from .instrument_type import get_instrument_type
from ..setup import config
logger = logging.getLogger(__name__.rsplit('.', maxsplit=1)[-1])


def _add_instrtype(trace):
    """
    Add instrtype to trace.

    :param trace: ObsPy Trace object
    :type trace: :class:`obspy.core.trace.Trace`
    """
    instrtype = get_instrument_type(trace)
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
    trace.stats.inventory = inv


def _check_instrtype(trace):
    """
    Check if instrument type is consistent with units in inventory.

    :param trace: ObsPy Trace object
    :type trace: :class:`obspy.core.trace.Trace`
    """
    # If traces are already corrected, instrument response may not be present
    # In this case, we just write the correct units to the trace headers
    instrtype = trace.stats.instrtype
    if config.correct_instrumental_response is False:
        # Override instrtype from config.trace_units, which may be
        # different from native instrtype inferred from channel code
        if config.trace_units:
            config.trace_units = config.trace_units.lower()
        if config.trace_units == 'vel':
            if instrtype not in ('shortp', 'broadb'):
                instrtype = 'broadb'
        elif config.trace_units in ('acc', 'disp'):
            instrtype = config.trace_units
        trace.stats.instrtype = instrtype
        if instrtype == 'disp':
            trace.stats.units = 'm'
        elif instrtype in ('shortp', 'broadb'):
            trace.stats.units = 'm/s'
        elif instrtype == 'acc':
            trace.stats.units = 'm/s**2'
        else:
            raise RuntimeError(
                f'{trace.id}: cannot get units for instrument type '
                f'{instrtype}.\n'
                f'> Trace start time: {trace.stats.starttime}\n'
                '> Skipping trace'
            )
    else:
        inv = trace.stats.inventory
        if not inv:
            raise RuntimeError(
                f'{trace.id}: cannot get instrtype from inventory: '
                'inventory is empty: skipping trace')
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
    if coords is None:
        # Give up!
        _add_coords.skipped.append(trace.id)
        raise RuntimeError(
            f'{trace.id}: could not find coords for trace: skipping trace')
    # Build an AttribDict and make sure that coordinates are floats
    coords = AttribDict({
        key: float(value) for key, value in coords.items()})
    coords = (
        None if (
            coords.longitude == 0 and coords.latitude == 0 and
            coords.local_depth == 123456 and coords.elevation == 123456)
        else coords)
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
    trace_picks = [
        pick for pick in picks if pick.station == trace.stats.station]
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


def augment_traces(stream, inventory, ssp_event, picks):
    """
    Add all required information to trace headers.

    Trace with no or incomplete metadata are skipped.

    :param stream: Traces to be augmented
    :type stream: :class:`obspy.core.stream.Stream`
    :param inventory: Station metadata
    :type inventory: :class:`obspy.core.inventory.Inventory`
    :param ssp_event: Event information
    :type ssp_event: :class:`sourcespec.ssp_event.SSPEvent`
    :param picks: list of picks
    :type picks: list of :class:`sourcespec.ssp_event.Pick`
    """
    traces_to_keep = []
    for trace in stream:
        try:
            _augment_trace(trace, inventory, ssp_event, picks)
        except Exception as err:
            for line in str(err).splitlines():
                logger.warning(line)
            continue
        traces_to_keep.append(trace)
    # in-place update of st
    stream.traces[:] = traces_to_keep[:]
    _complete_picks(stream)
