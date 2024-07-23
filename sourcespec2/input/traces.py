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
import re
import logging
import warnings
import contextlib
import shutil
import tarfile
import zipfile
import tempfile
from obspy import read
from obspy.core import Stream
from ..setup import config, ssp_exit
from .trace_parsers import parse_asdf_traces
from .station_metadata_parsers import get_instrument_from_SAC
from .instrument_type import get_instrument_type
logger = logging.getLogger(__name__.rsplit('.', maxsplit=1)[-1])
# Silence specific ObsPy SAC sample spacing warning (see obspy#3408)
warnings.filterwarnings(
    'ignore',
    message=r'.*Sample spacing read from SAC file',
    category=UserWarning,
    module=r'obspy\.io\.sac\.util'
)


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
    if getattr(config.options, 'station', None) is None:
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
    stream = Stream()
    asdf_path = getattr(config.options, 'asdf_path', None)
    if not asdf_path:
        return stream
    asdf_tag = getattr(config.options, 'asdf_tag', None)
    for asdf_file in asdf_path:
        logger.info(f'Reading traces from ASDF file: {asdf_file}')
        stream += _filter_by_station(
            parse_asdf_traces(asdf_file, tag=asdf_tag, read_headers=True)
        )
    return stream


def _read_trace_files():
    """
    Read trace files from the path specified in the configuration.

    :return: ObsPy Stream object
    :rtype: :class:`obspy.core.stream.Stream`
    """
    if getattr(config.options, 'trace_path', None) is None:
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


def _correct_traceids(stream):
    """
    Correct traceids from config.TRACEID_MAP, if available.

    :param stream: ObsPy Stream object
    :type stream: :class:`obspy.core.stream.Stream`
    """
    if config.TRACEID_MAP is None:
        return
    for trace in stream:
        with contextlib.suppress(KeyError):
            traceid = config.TRACEID_MAP[trace.get_id()]
            net, sta, loc, chan = traceid.split('.')
            trace.stats.network = net
            trace.stats.station = sta
            trace.stats.location = loc
            trace.stats.channel = chan


def _update_non_standard_trace_ids(stream):
    """
    Update non-standard trace IDs with a standard SEED ID obtained from the
    instrument type.

    :param stream: Stream object
    :type stream: :class:`obspy.core

    .. note::
        Currently only SAC files are supported.
    """
    traces_to_skip = []
    for trace in stream:
        if not hasattr(trace.stats, 'sac'):
            continue
        instrtype = get_instrument_type(trace)
        if instrtype is not None:
            continue
        try:
            instrtype, band_code, instr_code = get_instrument_from_SAC(trace)
        except RuntimeError as e:
            logger.warning(e)
            traces_to_skip.append(trace)
            continue
        old_id = trace.id
        orientation = trace.stats.channel[-1]
        trace.stats.channel = ''.join((band_code, instr_code, orientation))
        msg = f'{old_id}: non-standard trace ID updated to {trace.id}'
        if msg not in _update_non_standard_trace_ids.msgs:
            logger.info(msg)
            _update_non_standard_trace_ids.msgs.append(msg)
    stream.traces = [trace for trace in stream if trace not in traces_to_skip]
_update_non_standard_trace_ids.msgs = []  # noqa


def _glob_to_regex(pattern):
    """
    Convert a glob-style pattern to a regex pattern.
    """
    # Escape regex-special characters except for ? and *
    pattern = pattern.strip()
    pattern = re.escape(pattern)            # Escapes all regex chars
    pattern = pattern.replace(r'\?', '.')   # Convert escaped ? to .
    pattern = pattern.replace(r'\*', '.*')  # Convert escaped * to .*
    return pattern


def _should_keep_trace(traceid):
    """
    Check if trace should be kept.

    :param traceid: Trace ID.
    :type traceid: str

    :raises: RuntimeError if traceid is to be skipped.
    """
    network, station, location, channel = traceid.split('.')
    orientation_codes = config.vertical_channel_codes +\
        config.horizontal_channel_codes_1 +\
        config.horizontal_channel_codes_2
    orientation = channel[-1]
    if orientation not in orientation_codes:
        raise RuntimeError(
            f'{traceid}: Unknown channel orientation: "{orientation}"'
        )
    network, station, location, channel = traceid.split('.')

    # Build all possible IDs: station → full net.sta.loc.chan
    ss = [
        station,
        '.'.join((network, station)),
        '.'.join((network, station, location)),
        '.'.join((network, station, location, channel)),
    ]

    def _matches(patterns, strings):
        """Return True if any string matches any glob pattern."""
        regex = r'^(' + '|'.join(_glob_to_regex(p) for p in patterns) + r')$'
        return any(re.match(regex, s) for s in strings)

    if (
        config.use_traceids is not None and
        not _matches(config.use_traceids, ss)
    ):
        raise RuntimeError(f'{traceid}: not selected from config file')

    if (
        config.ignore_traceids is not None and
        _matches(config.ignore_traceids, ss)
    ):
        raise RuntimeError(f'{traceid}: ignored from config file')


def _select_requested_components(stream):
    """
    Select requested components from stream

    :param stream: ObsPy Stream object
    :type stream: :class:`obspy.core.stream.Stream`

    :return: ObsPy Stream object
    :rtype: :class:`obspy.core.stream.Stream`
    """
    traces_to_keep = []
    for trace in stream:
        try:
            _should_keep_trace(trace.id)
        except RuntimeError as e:
            logger.warning(str(e))
            continue
        traces_to_keep.append(trace)
    # in-place update of st
    stream.traces[:] = traces_to_keep[:]


def read_traces():
    """
    Read traces from the files or paths specified in the configuration.

    :return: Traces
    :rtype: :class:`obspy.core.stream.Stream`
    """
    logger.info('Reading traces...')
    stream = _read_asdf_traces() + _read_trace_files()
    _correct_traceids(stream)
    _update_non_standard_trace_ids(stream)
    _select_requested_components(stream)
    ntraces = len(stream)
    logger.info(f'Reading traces: {ntraces} traces loaded')
    logger.info('---------------------------------------------------')
    if not ntraces:
        logger.error('No trace loaded')
        ssp_exit(1)
    return stream
