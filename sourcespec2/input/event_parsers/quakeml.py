# -*- coding: utf-8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Parse event metadata and picks from a QuakeML file.

:copyright:
    2012-2025 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
import logging
from obspy import read_events
from .obspy_catalog import parse_obspy_catalog
logger = logging.getLogger(__name__.rsplit('.', maxsplit=1)[-1])


def parse_qml_event_picks(qml_file, event_id=None):
    """
    Parse event metadata and picks from a QuakeML file.

    :param qml_file: Path to the QuakeML file.
    :type qml_file: str
    :param event_id: event id
    :type event_id: str

    :return: a tuple of (SSPEvent, picks)
    :rtype: tuple
    """
    ssp_event = None
    picks = []
    if qml_file is None:
        return ssp_event, picks
    try:
        obspy_catalog = read_events(qml_file)
        return parse_obspy_catalog(obspy_catalog, event_id, qml_file)
    except Exception as err:
        logger.warning(err)
        return ssp_event, picks
