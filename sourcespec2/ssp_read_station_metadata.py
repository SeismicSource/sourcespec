# -*- coding: utf-8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Read station metadata in StationXML, dataless SEED, SEED RESP,
PAZ (SAC polezero format).

:copyright:
    2012-2025 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
import os
import re
import logging
from obspy import read_inventory
from obspy.core.inventory import Inventory, Network, Station, Channel, Response
from .config import config
logger = logging.getLogger(__name__.rsplit('.', maxsplit=1)[-1])


class PAZ():
    """Instrument response defined through poles and zeros."""
    zeros = []
    poles = []
    sensitivity = 1.
    _seedID = None
    network = None
    station = None
    location = None
    channel = None
    input_units = None
    linenum = None

    def __init__(self, file=None):
        """
        Init PAZ object.

        :param file: path to the PAZ file
        :type file: str
        """
        if file is not None:
            self._read(file)

    def __str__(self):
        return (
            f'PAZ: {self.seedID}'
            f'  zeros: {self.zeros}'
            f'  poles: {self.poles}'
            f'  sensitivity: {self.sensitivity}'
            f'  input_units: {self.input_units}'
        )

    @property
    def seedID(self):
        """Return the seedID."""
        return self._seedID

    @seedID.setter
    def seedID(self, seedID):
        try:
            self.network, self.station, self.location, self.channel =\
                seedID.split('.')
        except ValueError as e:
            raise ValueError(
                f'Invalid seedID "{seedID}". '
                'SeedID must be in the form NET.STA.LOC.CHAN'
            ) from e
        self._seedID = seedID
        self._guess_input_units()

    def _guess_input_units(self):
        """
        Guess the input units from the seedID.
        """
        if len(self.channel) < 3:
            return
        instr_code = self.channel[1]
        self.input_units = None
        if instr_code in config.INSTR_CODES_VEL:
            band_code = self.channel[0]
            # SEED standard band codes for velocity channels
            # https://ds.iris.edu/ds/nodes/dmc/data/formats/seed-channel-naming
            if band_code in ['B', 'C', 'D', 'E', 'F', 'G', 'H', 'S']:
                self.input_units = 'M/S'
        elif instr_code in config.INSTR_CODES_ACC:
            self.input_units = 'M/S**2'

    def _read(self, file):
        """Read a PAZ file."""
        # We cannot use the "with" statement here because we need to keep
        # self.lines alive as an iterator
        # pylint: disable=consider-using-with
        fp = open(file, 'r', encoding='ascii')  # sourcery skip
        # pylint: enable=consider-using-with
        self.lines = enumerate(fp, start=1)
        while True:
            try:
                self._parse_paz_file_lines()
            except StopIteration:
                break
            except Exception as e:
                raise TypeError(
                    f'Unable to parse file "{file}" as PAZ file. '
                    f'Error at line {self.linenum}: {e}'
                ) from e
        fp.close()

    def _parse_paz_file_lines(self):
        """Parse a line or a set of lines of a PAZ file."""
        self.linenum, line = next(self.lines)
        word = line.split()
        if not word:
            return
        what = word[0].lower()
        if what in ['poles', 'zeros']:
            nvalues = int(word[1])
            poles_zeros = []
            for _ in range(nvalues):
                self.linenum, line = next(self.lines)
                value = complex(*map(float, line.split()))
                poles_zeros.append(value)
            setattr(self, what, poles_zeros)
        elif what == 'constant':
            self.sensitivity = float(word[1])

    def to_inventory(self):
        """
        Convert PAZ object to an Inventory object.
        """
        resp = Response().from_paz(
            self.zeros, self.poles, stage_gain=self.sensitivity,
            input_units=self.input_units, output_units='COUNTS')
        resp.instrument_sensitivity.value = self.sensitivity
        channel = Channel(
            code=self.channel, location_code=self.location, response=resp,
            latitude=0, longitude=0, elevation=123456, depth=123456)
        station = Station(
            code=self.station, channels=[channel, ],
            latitude=0, longitude=0, elevation=123456)
        network = Network(
            code=self.network, stations=[station, ])
        return Inventory(networks=[network, ])


def _read_paz_file(file):
    """
    Read a paz file into an ``Inventory`` object.

    :note:
    - paz file must have ".pz" or ".paz" suffix (or no suffix)
    - paz file name (without prefix and suffix) can have
      the trace_id (NET.STA.LOC.CHAN) of the corresponding trace in the last
      part of his name (e.g., 20110208_1600.NOW.IV.CRAC.00.EHZ.paz),
      otherwise it will be treaten as a generic paz.

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
    paz = PAZ(file)
    try:
        paz.seedID = trace_id
    except ValueError:
        paz.seedID = 'XX.GENERIC.XX.XXX'
        logger.info(f'Using generic trace ID for PAZ file {file}')
    if paz.input_units is None:
        paz.input_units = 'M/S'
        logger.warning(
            f'Cannot find input units for ID "{paz.seedID}". '
            'Defaulting to M/S')
    return paz.to_inventory()


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
    logger.info('---------------------------------------------------')
    return inventory
