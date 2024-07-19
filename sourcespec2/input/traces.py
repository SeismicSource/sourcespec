# -*- coding: utf-8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Read traces in multiple formats.

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
import os
import logging
import shutil
import tarfile
import zipfile
import tempfile
from obspy import read
from obspy.core import Stream
from ..setup import config, ssp_exit
from .sac_header import get_event_from_SAC
from .trace_parsers import parse_asdf_traces
logger = logging.getLogger(__name__.rsplit('.', maxsplit=1)[-1])


def _build_filelist(path, filelist, tmpdir):
    """
    Build a list of files to read.

    :param path: Path to a file or directory
    :type path: str
    :param filelist: List of files to read
    :type filelist: list
    :param tmpdir: Temporary directory
    :type tmpdir: str
    """
    if os.path.isdir(path):
        listing = os.listdir(path)
        for filename in listing:
            fullpath = os.path.join(path, filename)
            _build_filelist(fullpath, filelist, tmpdir)
    else:
        try:
            # pylint: disable=unspecified-encoding consider-using-with
            open(path)
        except IOError as err:
            logger.error(err)
            return
        if tarfile.is_tarfile(path) and tmpdir is not None:
            with tarfile.open(path) as tar:
                try:
                    tar.extractall(path=tmpdir)
                except Exception as msg:
                    logger.warning(
                        f'{path}: Unable to fully extract tar archive: {msg}')
        elif zipfile.is_zipfile(path) and tmpdir is not None:
            with zipfile.ZipFile(path) as zipf:
                try:
                    zipf.extractall(path=tmpdir)
                except Exception as msg:
                    logger.warning(
                        f'{path}: Unable to fully extract zip archive: {msg}')
        else:
            filelist.append(path)


def _filter_by_station(input_stream):
    """
    Select traces for a given station, if specified in the configuration.

    :param input_stream: Input stream to filter
    :type tmpst: :class:`obspy.core.stream.Stream`

    :return: Stream object
    :rtype: :class:`obspy.core.stream.Stream`
    """
    if config.options.station is None:
        return input_stream
    return Stream([
        trace for trace in input_stream.traces
        if trace.stats.station == config.options.station
    ])


def _read_asdf_traces():
    """
    Read traces from ASDF file specified in the configuration.

    :return: ObsPy Stream object
    :rtype: :class:`obspy.core.stream.Stream`
    """
    asdf_file = config.options.asdf_file
    if not asdf_file:
        return Stream()
    return _filter_by_station(parse_asdf_traces(asdf_file))


def _read_trace_files():
    """
    Read trace files from the path specified in the configuration.

    :return: ObsPy Stream object
    :rtype: :class:`obspy.core.stream.Stream`
    """
    if config.options.trace_path is None:
        return Stream()
    # phase 1: build a file list
    # ph 1.1: create a temporary dir and run '_build_filelist()'
    #         to move files to it and extract all tar archives
    tmpdir = tempfile.mkdtemp()
    filelist = []
    for trace_path in config.options.trace_path:
        _build_filelist(trace_path, filelist, tmpdir)
    # ph 1.2: rerun '_build_filelist()' in tmpdir to add to the
    #         filelist all the extraceted files
    listing = os.listdir(tmpdir)
    for filename in listing:
        fullpath = os.path.join(tmpdir, filename)
        _build_filelist(fullpath, filelist, None)
    # phase 2: build a stream object from the file list
    st = Stream()
    for filename in sorted(filelist):
        try:
            st += _filter_by_station(read(filename, fsize=False))
        except (TypeError, FileNotFoundError):
            logger.warning(
                f'{filename}: Unable to read file as a trace: skipping')
            continue
    shutil.rmtree(tmpdir)
    return st


def _read_event_from_traces(stream):
    """
    Read event information from trace headers.
    The event information is stored in the trace.stats.event attribute.

    Currently supports only the SAC header.

    :param stream: ObsPy Stream object
    :type stream: :class:`obspy.core.stream.Stream`
    """
    for trace in stream:
        try:
            trace.stats.event = get_event_from_SAC(trace)
        except RuntimeError:
            continue


def read_traces():
    """
    Read traces from the files or paths specified in the configuration.

    :return: Traces
    :rtype: :class:`obspy.core.stream.Stream`
    """
    logger.info('Reading traces...')
    stream = _read_asdf_traces() + _read_trace_files()
    _read_event_from_traces(stream)
    logger.info('Reading traces: done')
    logger.info('---------------------------------------------------')
    if len(stream) == 0:
        logger.error('No trace loaded')
        ssp_exit(1)
    return stream
