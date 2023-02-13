# -*- coding: utf-8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Read station metadata in StationXML, dataless SEED, SEED RESP,
PAZ (SAC polezero format).

:copyright:
    2012-2023 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
import os
import re
import logging
from obspy import read_inventory
from obspy.core.inventory import Inventory, Network, Station, Channel, Response
logger = logging.getLogger(__name__.split('.')[-1])


def _parse_paz_file(file):
    lines = iter(open(file, 'r'))
    zeros = []
    poles = []
    constant = None
    linenumber = 0
    try:
        for line in lines:
            linenumber += 1
            word = line.split()
            if not word:
                continue
            if word[0] == 'ZEROS':
                nzeros = int(word[1])
                for _ in range(nzeros):
                    linenumber += 1
                    _zero = complex(*map(float, next(lines).split()))
                    zeros.append(_zero)
            if word[0] == 'POLES':
                npoles = int(word[1])
                for _ in range(npoles):
                    linenumber += 1
                    _pole = complex(*map(float, next(lines).split()))
                    poles.append(_pole)
            if word[0] == 'CONSTANT':
                constant = float(word[1])
    except Exception:
        msg = 'Unable to parse file "{}" as PAZ. '.format(file)
        msg += 'Parse error at line {}'.format(linenumber)
        raise TypeError(msg)
    if constant is None:
        msg = 'Unable to parse file "{}" as PAZ. '.format(file)
        msg += 'Cannot find a "CONSTANT" value'
        raise TypeError(msg)
    return zeros, poles, constant


def _read_paz_file(file):
    """
    Read a paz file into an ``Inventory``object.

    - paz file must have ".pz" or ".paz" suffix (or no suffix)
    - paz file name (without prefix and suffix) can have
      the trace_id (NET.STA.LOC.CHAN) of the corresponding trace in the last
      part of his name (e.g., 20110208_1600.NOW.IV.CRAC.00.EHZ.paz),
      otherwhise it will be treaten as a generic paz.
    """
    bname = os.path.basename(file)
    # strip .pz suffix, if there
    bname = re.sub('.pz$', '', bname)
    # strip .paz suffix, if there
    bname = re.sub('.paz$', '', bname)
    # we assume that the last four fields of bname
    # (separated by '.') are the trace_id
    trace_id = '.'.join(bname.split('.')[-4:])
    try:
        # check if trace_id is an actual trace ID
        net, sta, loc, chan = trace_id.split('.')
    except ValueError:
        # otherwhise, let's use this PAZ for a generic trace ID
        net = 'XX'
        sta = 'GENERIC'
        loc = 'XX'
        chan = 'XXX'
    zeros, poles, constant = _parse_paz_file(file)
    resp = Response().from_paz(
        zeros, poles, stage_gain=1, input_units='M/S', output_units='COUNTS')
    resp.instrument_sensitivity.value = constant
    channel = Channel(
        code=chan, location_code=loc, response=resp,
        latitude=0, longitude=0, elevation=123456, depth=123456)
    station = Station(
        code=sta, channels=[channel, ],
        latitude=0, longitude=0, elevation=123456)
    network = Network(code=net, stations=[station, ])
    inv = Inventory(networks=[network, ])
    return inv


def read_station_metadata(path):
    """
    Read station metadata into an ObsPy ``Inventory`` object.
    """
    if path is None:
        return None
    logger.info('Reading station metadata...')
    inventory = Inventory()
    if os.path.isdir(path):
        filelist = [os.path.join(path, file) for file in os.listdir(path)]
    else:
        filelist = [path, ]
    for file in sorted(filelist):
        if os.path.isdir(file):
            # we do not enter into subdirs of "path"
            continue
        logger.info('Reading station metadata from file: {}'.format(file))
        try:
            inventory += read_inventory(file)
        except Exception:
            msg1 = 'Unable to parse file "{}" as Inventory'.format(file)
            try:
                inventory += _read_paz_file(file)
            except Exception as msg2:
                logger.warning(msg1)
                logger.warning(msg2)
                continue
    if not inventory:
        inventory = None
    logger.info('Reading station metadata: done')
    return inventory
