# -*- coding: utf-8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Get the instrument type from the channel name of a trace.

:copyright:
    2023-2024 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
from ..setup import config


def get_instrument_type(trace):
    """
    Get the instrument type from the channel name of the trace.

    :param trace: ObsPy Trace object
    :type trace: :class:`obspy.core.trace.Trace`

    :return: Instrument type
    :rtype: str
    """
    instrtype = None
    _band_code = None
    _instr_code = None
    # First, try to get the instrtype from channel name
    chan = trace.stats.channel
    if len(chan) > 2:
        _band_code = chan[0]
        _instr_code = chan[1]
    if _instr_code in config.INSTR_CODES_VEL:
        # SEED standard band codes from higher to lower sampling rate
        # https://ds.iris.edu/ds/nodes/dmc/data/formats/seed-channel-naming/
        if _band_code in ['G', 'D', 'E', 'S']:
            instrtype = 'shortp'
        if _band_code in ['F', 'C', 'H', 'B']:
            instrtype = 'broadb'
    if _instr_code in config.INSTR_CODES_ACC:
        instrtype = 'acc'
    return instrtype
