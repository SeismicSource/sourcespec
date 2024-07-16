# -*- coding: utf-8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Parse event metadata and picks from a QuakeML file.

:copyright:
    2012-2024 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
import logging
from obspy import read_events
from ...setup import config, ssp_exit
from .obspy_catalog import parse_obspy_catalog
logger = logging.getLogger(__name__.rsplit('.', maxsplit=1)[-1])


def parse_qml_file():
    """
    Parse event metadata and picks from a QuakeML file.

    :return: a tuple of (SSPEvent, picks)
    :rtype: tuple
    """
    qml_file = config.options.qml_file
    if qml_file is None:
        ssp_event = None
        picks = []
        return ssp_event, picks
    try:
        obspy_catalog = read_events(qml_file)
        return parse_obspy_catalog(obspy_catalog, file_name=qml_file)
    except Exception as err:
        logger.error(err)
        ssp_exit(1)
