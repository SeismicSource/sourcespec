# -*- coding: utf-8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Read station metadata from ASDF files.

:copyright:
    2012-2024 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
import logging
from obspy.core.inventory import Inventory
from ...setup import ssp_exit
logger = logging.getLogger(__name__.rsplit('.', maxsplit=1)[-1])


def parse_asdf_inventory(asdf_file):
    """
    Read station metadata from ASDF file

    :param asdf_file: full path to ASDF file
    :type asdf_file: str

    :return: inventory
    :rtype: :class:`~obspy.core.inventory.inventory.Inventory`
    """
    try:
        # pylint: disable=import-outside-toplevel
        # pyasdf is not a hard dependency, so we import it here
        # and check for ImportError
        import pyasdf
    except ImportError:
        ssp_exit(
            'Error importing pyasdf. '
            'See https://seismicdata.github.io/pyasdf/ for installation '
            'instructions.'
        )
    inventory = Inventory()
    try:
        ds = pyasdf.ASDFDataSet(asdf_file, mode='r')
    except OSError:
        logger.warning(f'Unable to read ASDF file: {asdf_file}')
        return inventory
    for nw_stat_code in ds.waveforms.list():
        if 'StationXML' in ds.waveforms[nw_stat_code]:
            station_inv = ds.waveforms[nw_stat_code].StationXML
            inventory += station_inv
    ds._close()
    return inventory
