# -*- coding: utf-8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Read metadata from SAC file headers.

:copyright:
    2023-2025 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
import logging
import contextlib
from ...ssp_event import SSPEvent
from ...ssp_pick import SSPPick
logger = logging.getLogger(__name__.rsplit('.', maxsplit=1)[-1])


def read_event_from_SAC(trace):
    """
    Read event information from SAC header.

    :param trace: ObsPy trace object
    :type trace: :class:`obspy.core

    :return: Event information
    :rtype: :class:`ssp_event.SSPEvent`

    :raises: ValueError if trace is not a SAC trace
    :raises: RuntimeError if hypocenter information is not found in header
    """
    try:
        sac_hdr = trace.stats.sac
    except AttributeError as e:
        raise ValueError(f'{trace.id}: not a SAC trace') from e
    try:
        evla = sac_hdr['evla']
        evlo = sac_hdr['evlo']
        evdp = sac_hdr['evdp']
        begin = sac_hdr['b']
    except KeyError as e:
        raise RuntimeError(
            f'{trace.id}: cannot find hypocenter information in SAC header: '
            f'{e}'
        ) from e
    starttime = trace.stats.starttime
    try:
        tori = sac_hdr['o']
        origin_time = starttime + tori - begin
        # make a copy of origin_time and round it to the nearest second
        if origin_time.microsecond >= 500000:
            evid_time = (origin_time + 1).replace(microsecond=0)
        else:
            evid_time = origin_time.replace(microsecond=0)
    except Exception:
        origin_time = None
        # make a copy of starttime and round it to the nearest minute
        if starttime.second >= 30:
            evid_time = (starttime + 60).replace(second=0, microsecond=0)
        else:
            evid_time = starttime.replace(second=0, microsecond=0)
    # Check if kevnm is not empty and does not contain spaces
    # (if it has spaces, then kevnm is probably not an evid)
    kevnm = sac_hdr.get('kevnm', '').strip()
    if kevnm and ' ' not in kevnm:
        evid = kevnm
    else:
        # create evid from origin time
        evid = evid_time.strftime('%Y%m%d_%H%M%S')
    ssp_event = SSPEvent()
    ssp_event.event_id = evid
    ssp_event.hypocenter.latitude.value_in_deg = evla
    ssp_event.hypocenter.longitude.value_in_deg = evlo
    ssp_event.hypocenter.depth.value = evdp
    ssp_event.hypocenter.depth.units = 'km'
    ssp_event.hypocenter.origin_time = origin_time
    return ssp_event


def read_picks_from_SAC(trace):
    """
    Read picks from SAC header.

    :param trace: ObsPy trace object
    :type trace: :class:`obspy.core

    :return: List of picks
    :rtype: list of :class:`ssp_pick.SSPPick`

    :raises: ValueError if trace is not a SAC trace
    :raises: RuntimeError if pick information is not found in header
    """
    try:
        sac_hdr = trace.stats.sac
    except AttributeError as e:
        raise ValueError(f'{trace.id}: not a SAC trace') from e
    trace_picks = []
    pick_fields = (
        'a', 't0', 't1', 't2', 't3', 't4', 't5', 't6', 't7', 't8', 't9')
    for field in pick_fields:
        try:
            time = sac_hdr[field]
            begin = sac_hdr['b']
        except KeyError:
            continue
        pick = SSPPick()
        pick.station = trace.stats.station
        pick.time = trace.stats.starttime + time - begin
        # we will try to get the phase from the label later
        pick.phase = 'X'
        # now look at labels (ka, kt0, ...)
        label = ''
        with contextlib.suppress(Exception):
            label = sac_hdr[f'k{field}'].strip()
        if len(label) == 4:
            # label is something like 'IPU0' or 'ES 2'
            pick.flag = label[0]
            pick.phase = label[1]
            pick.polarity = label[2]
            pick.quality = label[3]
        if pick.phase.upper() not in 'PS' and len(label) > 0:
            # we assume that label starts with 'P' or 'S'
            pick.phase = label[0]
        else:
            # no label, use default phase for field
            default_phases = {'a': 'P', 't0': 'S'}
            pick.phase = default_phases.get(field, 'X')
        trace_picks.append(pick)
    return trace_picks
