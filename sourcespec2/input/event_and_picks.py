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
    parse_source_spec_event_file,
    parse_hypo71_hypocenter, parse_hypo71_picks, parse_hypo2000_file,
    parse_qml_event_picks, parse_asdf_event_picks,
    read_event_from_SAC, read_picks_from_SAC
)
logger = logging.getLogger(__name__.rsplit('.', maxsplit=1)[-1])


def _read_event_and_picks_from_stream(stream):
    """
    Read event and phase picks from a stream.

    :param stream: ObsPy Stream object containing event info and phase picks
    :type stream: :class:`obspy.core.stream.Stream`

    :return: (ssp_event, picks)
    :rtype: tuple of
        :class:`sourcespec.ssp_event.SSPEvent`,
        list of :class:`sourcespec.ssp_event.Pick`

    .. note::
        Currently only SAC files are supported.
    """
    ssp_event = None
    picks = []
    for trace in stream:
        if ssp_event is None:
            try:
                ssp_event = read_event_from_SAC(trace)
            except RuntimeError as err:
                _read_event_and_picks_from_stream.event_warnings.append(err)
        try:
            picks += read_picks_from_SAC(trace)
        except RuntimeError as err:
            _read_event_and_picks_from_stream.picks_warnings.append(err)
    return ssp_event, picks
_read_event_and_picks_from_stream.event_warnings = []  # noqa
_read_event_and_picks_from_stream.picks_warnings = []  # noqa


def _read_event_and_picks_from_ASDF_file_list(filelist):
    """
    Read event and phase picks from a list of ASDF files.

    :param filelist: List of ASDF files
    :type filelist: list of str

    :return: (ssp_event, picks, event_source)
    :rtype: tuple of
        :class:`sourcespec.ssp_event.SSPEvent`,
        list of :class:`sourcespec.ssp_event.Pick`,
        str
    """
    ssp_event = None
    picks = []
    event_source = None
    for asdf_file in filelist:
        _ssp_event, _picks = parse_asdf_event_picks(asdf_file)
        if ssp_event is None:
            ssp_event = _ssp_event
            event_source = asdf_file
        picks += _picks
    return ssp_event, picks, event_source


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
        logger.warning(msg)
    picks = []
    ssp_event = None
    file_format = None
    return ssp_event, picks, file_format


def _replace_event(old_event, new_event, old_source, new_source):
    """
    Replace old event with new event.

    :param old_event: Old SSPEvent object
    :type old_event: :class:`sourcespec.ssp_event.SSPEvent`
    :param new_event: New SSPEvent object
    :type new_event: :class:`sourcespec.ssp_event.SSPEvent`
    :param old_source: Old source
    :type old_source: str
    :param new_source: New source
    :type new_source: str

    :return: New SSPEvent object, new source
    :rtype: tuple
    """
    if new_event is None:
        return old_event, old_source
    if old_event is not None:
        logger.warning(
            f'Replacing event information found in {old_source} '
            f'with information from {new_source}'
        )
    return new_event, new_source


def _replace_picks(old_picks, new_picks, old_source, new_source):
    """
    Replace old picks with new picks.

    :param old_picks: Old picks
    :type old_picks: list of :class:`sourcespec.ssp_event.Pick`
    :param new_picks: New picks
    :type new_picks: list of :class:`sourcespec.ssp_event.Pick`
    :param old_source: Old source
    :type old_source: str
    :param new_source: New source
    :type new_source: str

    :return: New picks, new source
    :rtype: tuple
    """
    if not new_picks:
        return old_picks, old_source
    if old_picks:
        logger.warning(
            f'Replacing picks found in {old_source} '
            f'with picks from {new_source}'
        )
    return new_picks, new_source


def _log_event_and_pick_info(ssp_event, picks, event_source, picks_source):
    """
    Log event information.

    :param ssp_event: SSPEvent object
    :type ssp_event: :class:`sourcespec.ssp_event.SSPEvent`
    :param picks: List of Pick objects
    :type picks: list of :class:`sourcespec.ssp_event.Pick`
    :param event_source: Event source
    :type event_source: str
    :param picks_source: Picks source
    """
    if ssp_event is None:
        # only log warnings if no event information was found
        for warning in _read_event_and_picks_from_stream.event_warnings:
            logger.warning(warning)
    else:
        logger.info(f'Event information read from: {event_source}')
        for line in str(ssp_event).splitlines():
            logger.info(line)
    if not picks:
        # only log warnings if no pick information was found
        for warning in _read_event_and_picks_from_stream.picks_warnings:
            logger.warning(warning)
    else:
        logger.info(f'Pick information read from: {picks_source}')
        logger.info(f'{len(picks)} picks read')
    logger.info('---------------------------------------------------')


def read_event_and_picks(stream=None):
    """
    Read event and phase picks

    :param stream: ObsPy Stream object containing event info and phase picks
        (optional)
    :type stream: :class:`obspy.core.stream.Stream`

    :return: (ssp_event, picks)
    :rtype: tuple of
        :class:`sourcespec.ssp_event.SSPEvent`,
        list of :class:`sourcespec.ssp_event.Pick`

    .. note::
        The function reads event and phase picks from the following sources,
        from the less prioritary to the most prioritary:
        - SAC trace headers
        - ASDF file
        - QML file
        - Hypocenter file
        - Pick file
    """
    picks = []
    ssp_event = None
    event_source = None
    picks_source = None
    # first, try to read event and picks from stream (SAC trace headers)
    if stream is not None:
        ssp_event, picks = _read_event_and_picks_from_stream(stream)
        event_source = 'traces'
        picks_source = 'traces'
    asdf_path = getattr(config.options, 'asdf_path', None)
    hypo_file = getattr(config.options, 'hypo_file', None)
    pick_file = getattr(config.options, 'pick_file', None)
    qml_file = getattr(config.options, 'qml_file', None)
    # parse ASDF file, possibly replacing event and picks
    if asdf_path is not None:
        _new_ssp_event, _new_picks, _new_event_source =\
            _read_event_and_picks_from_ASDF_file_list(asdf_path)
        ssp_event, event_source = _replace_event(
            ssp_event, _new_ssp_event, event_source, _new_event_source
        )
        picks, picks_source = _replace_picks(
            picks, _new_picks, picks_source, 'ASDF files'
        )
    # parse QML file, possibly replacing event and picks
    if qml_file is not None:
        _new_ssp_event, _new_picks = parse_qml_event_picks(qml_file)
        ssp_event, event_source = _replace_event(
            ssp_event, _new_ssp_event, event_source, qml_file
        )
        picks, picks_source = _replace_picks(
            picks, _new_picks, picks_source, qml_file
        )
    # parse hypocenter file, possibly replacing event and picks
    if hypo_file is not None:
        _new_ssp_event, _new_picks, file_format = _parse_hypo_file(
            hypo_file, config.options.evid)
        # this is needed when writing the output file in hypo71 format
        config.hypo_file_format = file_format
        ssp_event, event_source = _replace_event(
            ssp_event, _new_ssp_event, event_source, hypo_file
        )
        picks, picks_source = _replace_picks(
            picks, _new_picks, picks_source, hypo_file
        )
    # parse pick file, possibly replacing picks
    if pick_file is not None:
        _new_picks = parse_hypo71_picks()
        picks, picks_source = _replace_picks(
            picks, _new_picks, picks_source, pick_file
        )
    _log_event_and_pick_info(ssp_event, picks, event_source, picks_source)

    if ssp_event is None:
        logger.error('No hypocenter information found.')
        sys.stderr.write(
            '\n'
            'Use "-q" or "-H" options to provide hypocenter information\n'
            'or add hypocenter information to the SAC file header\n'
            '(if you use the SAC format).\n'
        )
        ssp_exit(1)

    return ssp_event, picks
