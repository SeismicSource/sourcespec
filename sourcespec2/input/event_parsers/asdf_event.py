# -*- coding: utf-8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Parse event metadata and picks from ASDF file.

:copyright:
    2012-2026 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
import logging
from ...setup import ssp_exit
from .obspy_catalog import parse_obspy_catalog
logger = logging.getLogger(__name__.rsplit('.', maxsplit=1)[-1])


def parse_asdf_event_picks(asdf_file, event_id=None):
    """
    Parse event metadata and picks from ASDF file

    :param asdf_file: full path to ASDF file
    :type asdf_file: str
    :param event_id: event id
    :type event_id: str

    :return: a tuple of (SSPEvent, picks)
    :rtype: tuple
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
    try:
        obspy_catalog = pyasdf.ASDFDataSet(asdf_file, mode='r').events
        return parse_obspy_catalog(obspy_catalog, event_id, asdf_file)
    except Exception as err:
        logger.warning(err)
        ssp_event = None
        picks = []
        return ssp_event, picks
