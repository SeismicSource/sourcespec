# -*- coding: utf-8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Read station metadata in StationXML, dataless SEED, SEED RESP,
PAZ (SAC polezero format), ASDF

:copyright:
    2012-2024 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
import os
import logging
import contextlib
from obspy import read_inventory
from obspy.core.inventory import Inventory
from ..setup import config
from .station_metadata_parsers import (
    PAZ,
    read_paz_file, parse_asdf_inventory,
    compute_sensitivity_from_SAC,
    get_station_coordinates_from_SAC
)
logger = logging.getLogger(__name__.rsplit('.', maxsplit=1)[-1])


def _update_coords(obj, coords):
    """
    Update coordinates the given object.

    :param obj: ObsPy Station or Channel object
    :type obj: :class:`obspy.core.inventory.station.Station`
        or :class:`obspy.core.inventory.channel.Channel`
    :param coords: Dictionary containing coordinates
    :type coords: dict
    """
    obj.latitude = coords['latitude']
    obj.longitude = coords['longitude']
    obj.elevation = coords['elevation']
    with contextlib.suppress(KeyError):
        obj.depth = coords['local_depth']


def _update_inventory_with_trace_coords(trace, inventory):
    """
    Update inventory with station coordinates from a trace.

    The inventory is updated in place.

    :param trace: ObsPy trace object
    :type trace: :class:`obspy.core.trace.Trace`
    :param inventory: ObsPy Inventory object to update
    :type inventory: :class:`obspy.core

    .. note::
        Currently only SAC files are supported.
    """
    if not hasattr(trace.stats, 'sac'):
        return
    # get coordinates from SAC headers
    sac_coords = get_station_coordinates_from_SAC(trace)
    if sac_coords is None:
        return
    for net in inventory:
        if net.code != trace.stats.network:
            continue
        for sta in net:
            if sta.code != trace.stats.station:
                continue
            for chan in sta:
                if chan.location_code != trace.stats.location:
                    continue
                if chan.code != trace.stats.channel:
                    continue
                _update_coords(chan, sac_coords)
            _update_coords(sta, sac_coords)


def _update_inventory_with_trace_sensitivity(trace, inventory):
    """
    Update inventory with sensitivity from a trace.

    :param trace: ObsPy trace object
    :type trace: :class:`obspy.core.trace.Trace`
    :param inventory: ObsPy Inventory object to update
    :type inventory: :class:`obspy.core.inventory.inventory.Inventory`

    :return: inventory
    :rtype: :class:`~obspy.core.inventory.inventory.Inventory`

    .. note::
        Currently only SAC files are supported.
    """
    if not hasattr(trace.stats, 'sac'):
        return inventory
    update_inventory = inventory.copy()
    trace_inv = update_inventory.select(
        network=trace.stats.network,
        station=trace.stats.station,
        location=trace.stats.location,
        channel=trace.stats.channel,
        starttime=trace.stats.starttime
    )
    # save coordinates from the inventory, if available
    coords = None
    with contextlib.suppress(Exception):
        coords = trace_inv.get_coordinates(
            trace.id, trace.stats.starttime)
    paz = PAZ()
    paz.seedID = trace.id
    paz.sensitivity = compute_sensitivity_from_SAC(trace)
    paz.poles = []
    paz.zeros = []
    if trace_inv:
        logger.warning(
            f'Overriding response for {trace.id} with constant '
            f'sensitivity {paz.sensitivity}')
    # override the trace_inv object with a new one constructed from the
    # sensitivity value
    trace_inv = paz.to_inventory()
    # restore coordinates, if available
    if coords is not None:
        for sta in trace_inv[0].stations:
            _update_coords(sta, coords)
            for chan in sta.channels:
                _update_coords(chan, coords)
    else:
        _update_inventory_with_trace_coords(trace, trace_inv)
    # update inventory
    update_inventory = update_inventory.remove(
        network=trace.stats.network,
        station=trace.stats.station,
        location=trace.stats.location,
        channel=trace.stats.channel
    ) + trace_inv
    return update_inventory


def _update_inventory_from_stream(stream, inventory):
    """
    Update inventory with station metadata from a stream.

    :param stream: ObsPy Stream object containing station metadata
    :type stream: :class:`obspy.core.stream.Stream`
    :param inventory: ObsPy Inventory object to update
    :type inventory: :class:`obspy.core.inventory.inventory.Inventory`

    :return: inventory
    :rtype: :class:`~obspy.core.inventory.inventory.Inventory`

    .. note::
        Currently only SAC files are supported.
    """
    update_inventory = inventory.copy()
    for trace in stream:
        # ignore traces that are not in SAC format
        if not hasattr(trace.stats, 'sac'):
            continue
        _update_inventory_with_trace_coords(trace, update_inventory)
        # If a "sensitivity" config option is provided, override the Inventory
        # object with a new one constructed from the sensitivity value
        if config.sensitivity is not None:
            update_inventory = _update_inventory_with_trace_sensitivity(
                trace, update_inventory)
    return update_inventory


def _read_asdf_inventory():
    """
    Read station metadata from ASDF file specified in the configuration.

    :return: inventory
    :rtype: :class:`~obspy.core.inventory.inventory.Inventory`
    """
    inventory = Inventory()
    asdf_path = getattr(config.options, 'asdf_path', None)
    if not asdf_path:
        return inventory
    for asdf_file in asdf_path:
        logger.info(f'Reading station metadata from ASDF file: {asdf_file}')
        inventory += parse_asdf_inventory(asdf_file)
    return inventory


def _read_station_metadata_from_files():
    """
    Read station metadata into an ObsPy ``Inventory`` object.

    :return: inventory
    :rtype: :class:`~obspy.core.inventory.inventory.Inventory`

    .. note::
        - Station metadata can be in StationXML, dataless SEED, SEED RESP,
          PAZ (SAC polezero format) format.
        - An empty inventory is returned if the station metadata path is
          not set in the configuration or if no valid files are found
    """
    inventory = Inventory()
    metadata_path = config.station_metadata
    if not metadata_path:
        return inventory
    if os.path.isdir(metadata_path):
        filelist = [
            os.path.join(metadata_path, file)
            for file in os.listdir(metadata_path)
        ]
    else:
        filelist = [metadata_path, ]
    for file in sorted(filelist):
        if os.path.isdir(file):
            # we do not enter into subdirs of "path"
            continue
        logger.info(f'Reading station metadata from file: {file}')
        try:
            inventory += read_inventory(file)
        except Exception:
            msg1 = f'Unable to parse file "{file}" as Inventory'
            try:
                inventory += read_paz_file(file)
            except Exception as msg2:
                logger.warning(msg1)
                logger.warning(msg2)
                continue
    return inventory


def read_station_metadata(stream=None):
    """
    Read station metadata.

    :param stream: ObsPy Stream object containing station metadata (optional)
    :type stream: :class:`obspy.core.stream.Stream`

    :return: inventory
    :rtype: :class:`~obspy.core.inventory.inventory.Inventory`
    """
    logger.info('Reading station metadata...')
    inventory = _read_station_metadata_from_files() + _read_asdf_inventory()
    if stream is not None:
        inventory = _update_inventory_from_stream(stream, inventory)
    nstations = len(inventory.get_contents()['stations'])
    logger.info(f'Reading station metadata: {nstations} stations read')
    logger.info('---------------------------------------------------')
    return inventory
