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
from sourcespec.ssp_setup import instr_codes_vel, instr_codes_acc
logger = logging.getLogger(__name__.split('.')[-1])


def _parse_poles_and_zeros_from_paz_file(lines, what='poles'):
    """
    Parse poles or zeros from a PAZ file.

    :param lines: lines of the PAZ file
    :type lines: list of str
    :param what: 'poles' or 'zeros'
    :type what: str

    :return: list of poles or zeros
    :rtype: list of complex
    """
    if what not in ['poles', 'zeros']:
        msg = f'Invalid value for "what": {what}'
        raise ValueError(msg)
    # Change to iterator, so that we can use next()
    lines = iter(lines)
    poles_or_zeros = []
    linenum = 0
    for line in lines:
        linenum += 1
        try:
            word = line.split()
            if not word:
                continue
            if word[0].lower() == what:
                nvalues = int(word[1])
                for _ in range(nvalues):
                    linenum += 1
                    value = complex(*map(float, next(lines).split()))
                    poles_or_zeros.append(value)
        except Exception as e:
            msg = f'Parse error at line {linenum}: {e}'
            raise TypeError(msg) from e
    return poles_or_zeros


def _parse_constant_from_paz_file(lines):
    """
    Parse constant from a PAZ file.

    :param lines: lines of the PAZ file
    :type lines: list of str

    :return: constant
    :rtype: float
    """
    constant = None
    for linenum, line in enumerate(lines):
        try:
            word = line.split()
            if not word:
                continue
            if word[0] == 'CONSTANT':
                constant = float(word[1])
        except Exception as e:
            msg = f'Parse error at line {linenum}: {e}'
            raise TypeError(msg) from e
    if constant is None:
        msg = 'Cannot find a "CONSTANT" value'
        raise TypeError(msg)
    return constant


def _parse_paz_file(file):
    """
    Parse a PAZ file.

    :param file: path to the PAZ file
    :type file: str

    :return: zeros, poles, constant
    :rtype: list of complex, list of complex, float
    """
    try:
        lines = open(file, 'r').readlines()
        poles = _parse_poles_and_zeros_from_paz_file(lines, what='poles')
        zeros = _parse_poles_and_zeros_from_paz_file(lines, what='zeros')
        constant = _parse_constant_from_paz_file(lines)
    except Exception as e:
        msg = f'Unable to parse file "{file}" as PAZ: {e}'
        raise TypeError(msg) from e
    return zeros, poles, constant


def _find_input_units(trace_id):
    """
    Find the input units of a trace, based on its trace ID.

    :param traceid: trace ID
    :type traceid: str
    """
    chan = trace_id.split('.')[-1]
    if len(chan) < 3:
        raise ValueError(f'Cannot find input units for trace ID "{trace_id}"')
    band_code = chan[0]
    instr_code = chan[1]
    input_units = None
    if instr_code in instr_codes_vel:
        # SEED standard band codes from higher to lower sampling rate
        # https://ds.iris.edu/ds/nodes/dmc/data/formats/seed-channel-naming/
        if band_code in ['B', 'C', 'D', 'E', 'F', 'G', 'H', 'S']:
            input_units = 'M/S'
    elif instr_code in instr_codes_acc:
        input_units = 'M/S**2'
    if input_units is None:
        raise ValueError(f'Cannot find input units for trace ID "{trace_id}"')
    return input_units


def _read_paz_file(file):
    """
    Read a paz file into an ``Inventory``object.

    :note:
    - paz file must have ".pz" or ".paz" suffix (or no suffix)
    - paz file name (without prefix and suffix) can have
      the trace_id (NET.STA.LOC.CHAN) of the corresponding trace in the last
      part of his name (e.g., 20110208_1600.NOW.IV.CRAC.00.EHZ.paz),
      otherwhise it will be treaten as a generic paz.

    :param file: path to the PAZ file
    :type file: str

    :return: inventory
    :rtype: :class:`~obspy.core.inventory.inventory.Inventory`
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
        input_units = _find_input_units(trace_id)
    except ValueError:
        # otherwhise, let's use this PAZ for a generic trace ID
        net = 'XX'
        sta = 'GENERIC'
        loc = 'XX'
        chan = 'XXX'
        # We use M/S as default input units, if "trace_units" is specified
        # in the config file, we will change this later
        input_units = 'M/S'
        logger.info(f'Using generic trace ID for PAZ file {file}')
    zeros, poles, constant = _parse_paz_file(file)
    resp = Response().from_paz(
        zeros, poles, stage_gain=1,
        input_units=input_units, output_units='COUNTS')
    resp.instrument_sensitivity.value = constant
    channel = Channel(
        code=chan, location_code=loc, response=resp,
        latitude=0, longitude=0, elevation=123456, depth=123456)
    station = Station(
        code=sta, channels=[channel, ],
        latitude=0, longitude=0, elevation=123456)
    network = Network(code=net, stations=[station, ])
    return Inventory(networks=[network, ])


def read_station_metadata(path):
    """
    Read station metadata into an ObsPy ``Inventory`` object.

    :param path: path to the station metadata file or directory
    :type path: str

    :return: inventory
    :rtype: :class:`~obspy.core.inventory.inventory.Inventory`
    """
    inventory = Inventory()
    if path is None:
        return inventory
    logger.info('Reading station metadata...')
    if os.path.isdir(path):
        filelist = [os.path.join(path, file) for file in os.listdir(path)]
    else:
        filelist = [path, ]
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
                inventory += _read_paz_file(file)
            except Exception as msg2:
                logger.warning(msg1)
                logger.warning(msg2)
                continue
    logger.info('Reading station metadata: done')
    return inventory
