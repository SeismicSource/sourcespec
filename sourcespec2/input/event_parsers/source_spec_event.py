# -*- coding: utf-8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Parse event metadata from a SourceSpec event file.

:copyright:
    2012-2025 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
import warnings
import logging
import xml.etree.ElementTree as ET
import yaml
from ...ssp_event import SSPEvent
logger = logging.getLogger(__name__.rsplit('.', maxsplit=1)[-1])


def parse_source_spec_event_file(event_file, event_id=None):
    """
    Parse a SourceSpec Event File, which is a YAML file.

    :param event_file: path to SourceSpec event file
    :type event_file: str
    :param evid: event id
    :type evid: str

    :return: SSPEvent object
    :rtype: SSPEvent

    .. note::
        The returned picks list is empty, for consistency with other parsers.
    """
    try:
        with open(event_file, encoding='utf-8') as fp:
            events = yaml.safe_load(fp)
    except Exception as e:
        raise TypeError('Not a valid YAML file.') from e
    # XML is valid YAML, but obviously not a SourceSpec Event File
    try:
        root = ET.fromstring(events)
    except Exception:
        root = None
    if root:
        raise TypeError(
            'The file is an XML file, not a YAML file. '
            'Try reading it with the "-q" option (QuakeML format).')
    # raise TypeError if events is not a list
    if not isinstance(events, list):
        raise TypeError(
            'This is a valid YAML file, but it does not contain the key: '
            '"- event_id:", preceded by a dash (-).')
    # make sure that all the event_id are a string
    for ev in events:
        ev['event_id'] = str(ev['event_id'])
    if event_id is not None:
        _events = [ev for ev in events if ev.get('event_id') == event_id]
        try:
            event = _events[0]
        except IndexError as e:
            raise ValueError(
                f'Event {event_id} not found in {event_file}') from e
    else:
        event = events[0]
        if len(events) > 1:
            logger.warning(
                f'Found {len(events)} events in {event_file}. '
                'Using the first one.')
    try:
        with warnings.catch_warnings(record=True) as w:
            ssp_event = SSPEvent(event)
            if len(w) > 0:
                logger.warning(f'Warnings while parsing {event_file}:')
            for warning in w:
                logger.warning(warning.message)
    except Exception as e:
        raise TypeError(
            'This is a valid YAML file, but the following error occurred: '
            f'{e}.'
        ) from e
    # empty picks list, for consistency with other parsers
    picks = []
    return ssp_event, picks
