# -*- coding: utf-8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Read station metadata in StationXML, dataless SEED, SEED RESP,
PAZ (SAC polezero format), ASDF

:copyright:
    2012-2025 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
import os
import logging
from obspy import read_inventory
from obspy.core.inventory import Inventory
from ..setup import config
from .station_metadata_parsers import read_paz_file, parse_asdf_inventory
logger = logging.getLogger(__name__.rsplit('.', maxsplit=1)[-1])


def _read_asdf_inventory():
    """
    Read station metadata from ASDF file specified in the configuration.

    :return: inventory
    :rtype: :class:`~obspy.core.inventory.inventory.Inventory`
    """
    asdf_file = config.options.asdf_file
    if not asdf_file:
        return Inventory()
    logger.info(f'Reading station metadata from ASDF file: {asdf_file}')
    return parse_asdf_inventory(asdf_file)


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


def read_station_metadata():
    """
    Read station metadata.

    :return: inventory
    :rtype: :class:`~obspy.core.inventory.inventory.Inventory`
    """
    logger.info('Reading station metadata...')
    inventory = _read_station_metadata_from_files() + _read_asdf_inventory()
    nstations = len(inventory.get_contents()['stations'])
    logger.info(f'Reading station metadata: {nstations} stations read')
    logger.info('---------------------------------------------------')
    return inventory
