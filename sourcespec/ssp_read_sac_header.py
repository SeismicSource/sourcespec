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
import re
import logging
import contextlib
from obspy.core.util import AttribDict
from sourcespec.ssp_setup import ssp_exit
from sourcespec.ssp_event import SSPEvent
from sourcespec.ssp_pick import SSPPick
logger = logging.getLogger(__name__.rsplit('.', maxsplit=1)[-1])


def compute_sensitivity_from_SAC(trace, config):
    """Compute sensitivity from SAC header fields."""
    # Securize the string before calling eval()
    # see https://stackoverflow.com/a/25437733/2021880
    inp = re.sub(r'\.(?![0-9])', '', config.sensitivity)
    # Try returning a float conversion of the string
    with contextlib.suppress(ValueError):
        return float(inp)
    # If the above fails, try evaluating the string as a combination of
    # SAC header fields
    try:
        namespace = trace.stats.sac
    except Exception as e:
        raise TypeError(
            f'{trace.id}: non-numerical values in the "sensitivity" config '
            'option are only allowed for SAC traces.'
        ) from e
    try:
        sensitivity = eval(inp, {}, namespace)  # pylint: disable=eval-used
    except NameError as msg:
        hdr_field = str(msg).split()[1]
        logger.error(
            f'SAC header field {hdr_field} in the "sensitivity" config '
            'option does not exist')
        ssp_exit(1)
    return sensitivity


# handpicked list of instruments, instrtypes and band/instr codes,
# mainly for ISNet compatibility
instruments = {
    'CMG-5T': {'instrtype': 'acc', 'band_code': 'H', 'instr_code': 'N'},
    'CMG-40T': {'instrtype': 'broadb', 'band_code': 'H', 'instr_code': 'H'},
    'TRILLIUM': {'instrtype': 'broadb', 'band_code': 'H', 'instr_code': 'H'},
    'S13J': {'instrtype': 'shortp', 'band_code': 'S', 'instr_code': 'H'},
    'KS2000ED': {'instrtype': 'shortp', 'band_code': 'S', 'instr_code': 'H'},
}


def get_instrument_from_SAC(trace):
    """Get instrument information from SAC header."""
    try:
        codes = instruments[trace.stats.sac.kinst]
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
    """Get station coordinates from SAC header."""
    with contextlib.suppress(Exception):
        stla = trace.stats.sac.stla
        stlo = trace.stats.sac.stlo
        stel = trace.stats.sac.get('stel', 0.)
        coords = AttribDict(latitude=stla, longitude=stlo, elevation=stel)
        logger.info(
            f'{trace.id}: station coordinates read from SAC header')
        return coords
    return None


def get_event_from_SAC(trace):
    """Get event information from SAC header."""
    try:
        sac_hdr = trace.stats.sac
    except AttributeError as e:
        raise RuntimeError(
            f'{trace.id}: not a SAC trace: cannot get hypocenter from header'
        ) from e
    evla = sac_hdr['evla']
    evlo = sac_hdr['evlo']
    evdp = sac_hdr['evdp']
    begin = sac_hdr['b']
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


def get_picks_from_SAC(trace):
    """Get picks from SAC header."""
    try:
        sac_hdr = trace.stats.sac
    except AttributeError as e:
        raise RuntimeError(
            f'{trace.id}: not a SAC trace: cannot get picks from header'
        ) from e
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
