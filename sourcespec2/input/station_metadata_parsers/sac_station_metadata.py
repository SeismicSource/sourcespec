# -*- coding: utf-8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Read station metadata from SAC file headers.

:copyright:
    2023-2024 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
import re
import logging
import contextlib
from obspy.core.util import AttribDict
from ...setup import config, ssp_exit
logger = logging.getLogger(__name__.rsplit('.', maxsplit=1)[-1])


def compute_sensitivity_from_SAC(trace):
    """
    Compute sensitivity from SAC header fields.

    :param trace: ObsPy trace object
    :type trace: :class:`obspy.core.trace.Trace`

    :return: Sensitivity
    :rtype: float
    """
    # Securize the string before calling eval()
    # see https://stackoverflow.com/a/25437733/2021880
    inp = re.sub(r'\.(?![0-9])', '', config.sensitivity)
    # Check if string contains letters, meaning that
    # it must contain SAC header fields and trace must be in SAC format
    namespace = None
    if re.search(r'[a-zA-Z]', inp):
        try:
            namespace = trace.stats.sac
        except Exception as e:
            raise TypeError(f'{trace.id}: trace must be in SAC format') from e
    try:
        sensitivity = eval(inp, {}, namespace)  # pylint: disable=eval-used
    except NameError as msg:
        hdr_field = str(msg).split()[1]
        logger.error(f'SAC header field {hdr_field} does not exist')
        ssp_exit(1)
    return sensitivity


# handpicked list of instruments, instrtypes and band/instr codes,
# mainly for ISNet compatibility
_instruments = {
    'CMG-5T': {'instrtype': 'acc', 'band_code': 'H', 'instr_code': 'N'},
    'CMG-40T': {'instrtype': 'broadb', 'band_code': 'H', 'instr_code': 'H'},
    'TRILLIUM': {'instrtype': 'broadb', 'band_code': 'H', 'instr_code': 'H'},
    'S13J': {'instrtype': 'shortp', 'band_code': 'S', 'instr_code': 'H'},
    'KS2000ED': {'instrtype': 'shortp', 'band_code': 'S', 'instr_code': 'H'},
}


def get_instrument_from_SAC(trace):
    """
    Get instrument information from SAC header.

    :param trace: ObsPy trace object
    :type trace: :class:`obspy.core

    :return: Instrument type, band code, instrument code
    :rtype: tuple
    """
    try:
        codes = _instruments[trace.stats.sac.kinst]
        instrtype = codes['instrtype']
        band_code = codes['band_code']
        instr_code = codes['instr_code']
    except AttributeError as e:
        raise RuntimeError(
            f'{trace.id}: cannot find instrtype for trace: '
            'missing SAC header field "kinst": skipping trace'
        ) from e
    except KeyError as e:
        raise RuntimeError(
            f'{trace.id}: cannot find instrtype for trace: '
            f'unknown instrument "{trace.stats.sac.kinst}": '
            'skipping trace'
        ) from e
    return instrtype, band_code, instr_code


def get_station_coordinates_from_SAC(trace):
    """
    Get station coordinates from SAC header.

    :param trace: ObsPy trace object
    :type trace: :class:`obspy.core

    :return: Station coordinates
    :rtype: :class:`AttribDict` or None
    """
    with contextlib.suppress(Exception):
        stla = trace.stats.sac.stla
        stlo = trace.stats.sac.stlo
        stel = trace.stats.sac.get('stel', 0.)
        coords = AttribDict(latitude=stla, longitude=stlo, elevation=stel)
        msg = f'{trace.id}: station coordinates read from SAC header'
        if msg not in get_station_coordinates_from_SAC.msgs:
            logger.info(msg)
            get_station_coordinates_from_SAC.msgs.append(msg)
        return coords
    return None
get_station_coordinates_from_SAC.msgs = []  # noqa
