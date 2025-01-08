# -*- coding: utf-8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Parse event metadata and picks from an ObsPy catalog object.

:copyright:
    2012-2025 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
import contextlib
import logging
import re
from ...setup import config, ssp_exit
from ...ssp_event import SSPEvent
from ...ssp_pick import SSPPick
logger = logging.getLogger(__name__.rsplit('.', maxsplit=1)[-1])


def _get_evid_from_resource_id(resource_id):
    """
    Get evid from resource_id.

    :param resource_id: resource_id string
    :type resource_id: str

    :returns: evid string
    :rtype: str
    """
    evid = resource_id
    if '/' in evid:
        evid = resource_id.split('/')[-1]
    if '?' in evid:
        evid = resource_id.split('?')[-1]
    if '&' in evid:
        evid = evid.split('&')[0]
    if '=' in evid:
        evid = evid.split('=')[-1]
    return evid


def _parse_event_metadata(obspy_event):
    """
    Parse event metadata from an ObsPy event object.

    :param obspy_event: ObsPy event object
    :type obspy_event: obspy.core.event.Event

    :return: a tuple of SSPEvent object and ObsPy origin object
    :rtype: (SSPEvent, obspy.core.event.Origin)
    """
    # No need to parse event name from the ObsPy event if event name is given
    # in the command line
    if getattr(config.options, 'evname', None):
        parse_event_name_from_description = False
        event_description_regex = None
    else:
        parse_event_name_from_description = config.qml_event_description
        event_description_regex = config.qml_event_description_regex
    ssp_event = SSPEvent()
    ssp_event.event_id = _get_evid_from_resource_id(
        str(obspy_event.resource_id.id))
    if parse_event_name_from_description:
        try:
            ssp_event.name = str(obspy_event.event_descriptions[0].text)
        except IndexError:
            logger.warning(
                'The event does not contain a description. Cannot parse '
                'event name from description.'
            )
    if ssp_event.name and event_description_regex:
        pattern = re.compile(event_description_regex)
        match = pattern.search(ssp_event.name)
        if match:
            name = match.group()
            # capitalize first letter
            name = name[0].upper() + name[1:]
            ssp_event.name = name
    # See if there is a preferred origin...
    obspy_origin = obspy_event.preferred_origin()
    # ...or just use the first one
    if obspy_origin is None:
        obspy_origin = obspy_event.origins[0]
    ssp_event.hypocenter.longitude.value_in_deg = obspy_origin.longitude
    ssp_event.hypocenter.latitude.value_in_deg = obspy_origin.latitude
    ssp_event.hypocenter.depth.value = obspy_origin.depth
    ssp_event.hypocenter.depth.units = 'm'
    ssp_event.hypocenter.origin_time = obspy_origin.time
    return ssp_event, obspy_origin


def _parse_magnitude_from_obspy_event(obspy_event, ssp_event):
    """
    Parse magnitude from an ObsPy event.

    :param obspy_event: ObsPy event
    :type obspy_event: obspy.core.event.Event
    :param ssp_event: SSPEvent object
    :type ssp_event: SSPEvent
    """
    mag = obspy_event.preferred_magnitude() or obspy_event.magnitudes[0]
    ssp_event.magnitude.value = mag.mag
    ssp_event.magnitude.mag_type = mag.magnitude_type


def _parse_scalar_moment_from_obspy_event(obspy_event, ssp_event):
    """
    Parse scalar moment from an ObsPy event.

    :param obspy_event: ObsPy event
    :type obspy_event: obspy.core.event.Event
    :param ssp_event: SSPEvent object
    :type ssp_event: SSPEvent
    """
    fm = obspy_event.preferred_focal_mechanism()\
        or obspy_event.focal_mechanisms[0]
    ssp_event.scalar_moment.value = fm.moment_tensor.scalar_moment
    ssp_event.scalar_moment.units = 'N-m'


def _parse_moment_tensor_from_obspy_event(obspy_event, ssp_event):
    """
    Parse moment tensor from an ObsPy event.

    :param obspy_event: ObsPy event
    :type obspy_event: obspy.core.event.Event
    :param ssp_event: SSPEvent object
    :type ssp_event: SSPEvent
    """
    fm = obspy_event.preferred_focal_mechanism()\
        or obspy_event.focal_mechanisms[0]
    mt = fm.moment_tensor.tensor
    ssp_event.moment_tensor.m_rr = mt.m_rr
    ssp_event.moment_tensor.m_tt = mt.m_tt
    ssp_event.moment_tensor.m_pp = mt.m_pp
    ssp_event.moment_tensor.m_rt = mt.m_rt
    ssp_event.moment_tensor.m_rp = mt.m_rp
    ssp_event.moment_tensor.m_tp = mt.m_tp
    ssp_event.moment_tensor.units = 'N-m'


def _parse_focal_mechanism_from_obspy_event(obspy_event, ssp_event):
    """
    Parse focal mechanism from an ObsPy event.

    :param obspy_event: ObsPy event
    :type obspy_event: obspy.core.event.Event
    :param ssp_event: SSPEvent object
    :type ssp_event: SSPEvent
    """
    fm = obspy_event.focal_mechanisms[0]
    nodal_plane = fm.nodal_planes.nodal_plane_1
    ssp_event.focal_mechanism.strike = nodal_plane.strike
    ssp_event.focal_mechanism.dip = nodal_plane.dip
    ssp_event.focal_mechanism.rake = nodal_plane.rake


def _parse_picks_from_obspy_event(obspy_event, obspy_origin):
    """
    Parse picks from an ObsPy event.

    :param obspy_event: ObsPy event
    :type obspy_event: obspy.core.event.Event
    :param obspy_origin: ObsPy origin object
    :type obspy_origin: obspy.core.event.Origin

    :return: list of SSPPick objects
    :rtype: list
    """
    picks = []
    for pck in obspy_event.picks:
        pick = SSPPick()
        pick.station = pck.waveform_id.station_code
        pick.network = pck.waveform_id.network_code
        pick.channel = pck.waveform_id.channel_code
        if pck.waveform_id.location_code is not None:
            pick.location = pck.waveform_id.location_code
        else:
            pick.location = ''
        if pck.onset == 'emergent':
            pick.flag = 'E'
        elif pck.onset == 'impulsive':
            pick.flag = 'I'
        pick.phase = None
        if pck.phase_hint:
            pick.phase = pck.phase_hint[:1]
        else:
            # try to get the phase from the arrival object that uses this pick
            pick_id = pck.resource_id.id
            arrivals = [
                arr for arr in obspy_origin.arrivals
                if arr.pick_id.id == pick_id
            ]
            if arrivals:
                pick.phase = arrivals[0].phase[:1]
        if pick.phase is None:
            # ignore picks with no phase hint
            continue
        if pck.polarity == 'negative':
            pick.polarity = 'D'
        elif pck.polarity == 'positive':
            pick.polarity = 'U'
        pick.time = pck.time
        picks.append(pick)
    return picks


def _get_event_from_obspy_catalog(
        obspy_catalog, event_id=None, file_name=''):
    """
    Get an event from an ObsPy catalog object.

    :param obspy_catalog: ObsPy catalog object
    :type obspy_catalog: instance of :class:`obspy.core.event.Catalog`
    :param event_id: event id
    :type event_id: str
    :param file_name: name of the file containing the catalog
    :type file_name: str

    :return: ObsPy event object
    :rtype: obspy.core.event.Event
    """
    if event_id is not None:
        _obspy_events = [
            ev for ev in obspy_catalog if event_id in str(ev.resource_id)
        ]
        try:
            obspy_event = _obspy_events[0]
        except IndexError as e:
            raise ValueError(
                f'Event {event_id} not found in {file_name}') from e
    else:
        obspy_event = obspy_catalog[0]
        if len(obspy_catalog) > 1:
            logger.warning(
                f'Found {len(obspy_catalog)} events in {file_name}. '
                'Using the first one.')
    return obspy_event


def _parse_obspy_event(obspy_event):
    """
    Parse event metadata and picks from an ObsPy event object.

    :param obspy_event: ObsPy event object to parse
    :type obspy_event: instance of :class:`obspy.core.event.Event`

    :return: a tuple of (SSPEvent, picks)
    :rtype: tuple
    """
    try:
        ssp_event, obspy_origin = _parse_event_metadata(obspy_event)
        picks = _parse_picks_from_obspy_event(obspy_event, obspy_origin)
    except Exception as err:
        ssp_exit(err)
    log_messages = []
    with contextlib.suppress(Exception):
        _parse_moment_tensor_from_obspy_event(obspy_event, ssp_event)
        # Compute focal mechanism, scalar moment and magnitude.
        # They will be overwritten later on, if they are found in the
        # ObsPy event.
        ssp_event.focal_mechanism.from_moment_tensor(ssp_event.moment_tensor)
        ssp_event.scalar_moment.from_moment_tensor(ssp_event.moment_tensor)
        ssp_event.magnitude.from_scalar_moment(ssp_event.scalar_moment)
    with contextlib.suppress(Exception):
        _parse_scalar_moment_from_obspy_event(obspy_event, ssp_event)
        # Compute magnitude from scalar moment. It will be overwritten later
        # on, if it is found in the ObsPy event.
        ssp_event.magnitude.from_scalar_moment(ssp_event.scalar_moment)
    with contextlib.suppress(Exception):
        _parse_magnitude_from_obspy_event(obspy_event, ssp_event)
    with contextlib.suppress(Exception):
        _parse_focal_mechanism_from_obspy_event(obspy_event, ssp_event)
    for msg in log_messages:
        logger.info(msg)
    return ssp_event, picks


def parse_obspy_catalog(obspy_catalog, event_id=None, file_name=''):
    """
    Parse event metadata and picks from an ObsPy catalog object.

    :param obspy_catalog: ObsPy catalog object to parse
    :type obspy_catalog: instance of :class:`obspy.core.event.Catalog`
    :param event_id: event id to extract from the catalog. If None, the value
        of config.options.evid is used. If this is also None,
        the first event in the catalog is used.
    :type event_id: str
    :param file_name: name of the file containing the catalog
    :type file_name: str

    :return: a tuple of (SSPEvent, picks)
    :rtype: tuple
    """
    event_id = event_id or getattr(config.options, 'evid', None)
    try:
        obspy_event = _get_event_from_obspy_catalog(
            obspy_catalog, event_id, file_name)
    except Exception as err:
        ssp_exit(err)
    else:
        ssp_event, picks = _parse_obspy_event(obspy_event)
    return ssp_event, picks
