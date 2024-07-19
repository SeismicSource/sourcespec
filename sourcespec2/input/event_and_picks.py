# -*- coding: utf-8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Read event and phase picks.

:copyright:
    2012 Claudio Satriano <satriano@ipgp.fr>

    2013-2014 Claudio Satriano <satriano@ipgp.fr>,
              Emanuela Matrullo <matrullo@geologie.ens.fr>

    2015-2024 Claudio Satriano <satriano@ipgp.fr>,
              Sophie Lambotte <sophie.lambotte@unistra.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
import sys
import logging

from ..setup import config, ssp_exit
from .event_parsers import (
    parse_hypo71_hypocenter,
    parse_hypo71_picks,
    parse_hypo2000_file,
    parse_source_spec_event_file,
    parse_qml_event_picks,
    parse_asdf_event_picks,
    override_event_depth
)
logger = logging.getLogger(__name__.rsplit('.', maxsplit=1)[-1])


# pylint: disable=inconsistent-return-statements
def _parse_hypo_file(hypo_file, event_id=None):
    """
    Parse a SourceSpec Event File, hypo71 or hypo2000 hypocenter file.

    :param hypo_file: Path to the hypocenter file.
    :type hypo_file: str
    :param event_id: Event ID.
    :type event_id: str

    :return: A tuple of (SSPEvent, picks, format).
    :rtype: tuple
    """
    err_msgs = []
    parsers = {
        'ssp_event_file': parse_source_spec_event_file,
        'hypo71': parse_hypo71_hypocenter,
        'hypo2000': parse_hypo2000_file,
    }
    format_strings = {
        'ssp_event_file': 'SourceSpec Event File',
        'hypo71': 'hypo71 hypocenter file',
        'hypo2000': 'hypo2000 hypocenter file',
    }
    for file_format, parser in parsers.items():
        try:
            ssp_event, picks = parser(hypo_file, event_id)
            return ssp_event, picks, file_format
        except Exception as err:
            format_str = format_strings[file_format]
            msg = f'{hypo_file}: Not a {format_str}'
            err_msgs.append(msg)
            msg = f'Parsing error: {err}'
            err_msgs.append(msg)
    # If we arrive here, the file was not recognized as valid
    for msg in err_msgs:
        logger.error(msg)
    ssp_exit(1)


def _log_event_info(ssp_event):
    """
    Log event information.

    :param ssp_event: SSPEvent object
    :type ssp_event: :class:`sourcespec.ssp_event.SSPEvent`
    """
    for line in str(ssp_event).splitlines():
        logger.info(line)
    logger.info('---------------------------------------------------')


def read_event_and_picks(trace1=None):
    """
    Read event and phase picks

    :param trace1: ObsPy Trace object containing event info (optional)
    :type trace1: :class:`obspy.core.stream.Stream`

    :return: (ssp_event, picks)
    :rtype: tuple of
        :class:`sourcespec.ssp_event.SSPEvent`,
        list of :class:`sourcespec.ssp_event.Pick`
    """
    picks = []
    ssp_event = None
    depth_override = getattr(config.options, 'depth', None)
    # parse hypocenter file
    if config.options.hypo_file is not None:
        ssp_event, picks, file_format = _parse_hypo_file(
            config.options.hypo_file, config.options.evid)
        config.hypo_file_format = file_format
    # parse pick file
    if config.options.pick_file is not None:
        picks = parse_hypo71_picks()
    # parse QML file
    if config.options.qml_file is not None:
        ssp_event, picks = parse_qml_event_picks(config.options.qml_file)
    # parse ASDF file
    if config.options.asdf_file is not None:
        ssp_event, picks = parse_asdf_event_picks(config.options.asdf_file)
    if ssp_event is not None:
        override_event_depth(ssp_event, depth_override)
        _log_event_info(ssp_event)

    # if ssp_event is still None, get it from first trace
    if ssp_event is None and trace1 is not None:
        try:
            ssp_event = trace1.stats.event
            _log_event_info(ssp_event)
        except AttributeError:
            logger.error('No hypocenter information found.')
            sys.stderr.write(
                '\n'
                'Use "-q" or "-H" options to provide hypocenter information\n'
                'or add hypocenter information to the SAC file header\n'
                '(if you use the SAC format).\n'
            )
            ssp_exit(1)
    # TODO: log also if trace1 is None?

    return (ssp_event, picks)
