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
import warnings
import shutil
import tarfile
import zipfile
import tempfile
import json
from obspy import read
from obspy.core import Stream
from ..setup import config, ssp_exit
from .sac_header import get_event_from_SAC
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


def _parse_asdf_trace_headers(ds, st, nw_stat_codes, tag):
    """
    Parse ASDF trace headers.

    :param ds: ASDF dataset
    :type ds: :class:`pyasdf.ASDFDataSet`
    :param st: ObsPy Stream object
    :type st: :class:`obspy.core.stream.Stream`
    :param nw_stat_codes: Network and station codes
    :type nw_stat_codes: list
    :param tag: waveform tag in ASDF file
    :type tag: str
    """
    # pylint: disable=import-outside-toplevel
    from pyasdf.utils import AuxiliaryDataContainer
    for nw_stat_code, tr in zip(nw_stat_codes, st.traces):
        if nw_stat_code in ds.auxiliary_data['TraceHeaders']:
            auxiliary_root = ds.auxiliary_data['TraceHeaders'][nw_stat_code]
        else:
            # ESM
            _nw_stat_code = nw_stat_code.replace('.', '_')
            auxiliary_root = ds.auxiliary_data['TraceHeaders'][_nw_stat_code]
        if not tag or tag not in auxiliary_root:
            continue
        header = None
        tr_id = tr.id.replace('.', '_')
        if isinstance(auxiliary_root[tag], AuxiliaryDataContainer):
            # ESM
            header = auxiliary_root[tag].parameters
        elif tr_id in auxiliary_root[tag]:
            header = auxiliary_root[tag][tr_id].parameters
        if not header:
            continue
        for key, val in header.items():
            try:
                val = json.loads(val)
            except json.JSONDecodeError:
                if isinstance(val, type(b'')):
                    # ESM
                    val = val.decode('ascii')
                if key == 'processing':
                    val = [val]
            # Try preserving original _format
            if (
                key == '_format'
                and val != 'ASDF'
                and 'original_format' not in header.items()
            ):
                key = 'original_format'
            # Write to trace header
            if key not in tr.stats:
                tr.stats[key] = val


def parse_asdf_traces(asdf_file, tag=None, read_headers=False):
    """
    Read all traces from ASDF file with given waveform tag
    If tag is not specified, the first available one will be taken

    :param asdf_file: full path to ASDF file
    :type asdf_file: str
    :param tag: waveform tag in ASDF file
    :type tag: str
    :param read_headers: flag to control reading of (non-standard)
        trace headers
    :type read_headers: bool

    :return: ObsPy Stream object
    :rtype: :class:`obspy.core.stream.Stream`
    """
    try:
        # pylint: disable=import-outside-toplevel
        # pyasdf is not a hard dependency, so we import it here
        # and check for ImportError
        import pyasdf
    except ImportError:
        logger.error(
            'Error importing pyasdf. '
            'See https://seismicdata.github.io/pyasdf/ for installation '
            'instructions.'
        )
        ssp_exit(1)
    try:
        ds = pyasdf.ASDFDataSet(asdf_file, mode='r')
    except Exception as err:
        logger.error(err)
        ssp_exit(1)
    # Read waveform data
    st = Stream()
    nw_stat_codes = []
    for nw_stat_code in ds.waveforms.list():
        wf_tags = ds.waveforms[nw_stat_code].get_waveform_tags()
        # If tag is not specified, take first available tag
        if not tag:
            # Maybe this should be logged
            tag = wf_tags[0]
        if tag in wf_tags:
            station_st = ds.waveforms[nw_stat_code][tag]
            st.extend(station_st)
            nw_stat_codes.extend([nw_stat_code] * len(station_st))
    # Try reading trace headers if present
    if 'TraceHeaders' in ds.auxiliary_data:
        header_key = 'TraceHeaders'
    elif 'Headers' in ds.auxiliary_data:
        header_key = 'Headers'
    else:
        header_key = None
    if read_headers and header_key:
        _parse_asdf_trace_headers(ds, st, nw_stat_codes, tag)
    ds._close()
    return st


def _read_trace_files():
    """
    Read trace files from the path specified in the configuration file.

    :return: ObsPy Stream object
    :rtype: :class:`obspy.core.stream.Stream`
    """
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
            if os.path.splitext(filename)[-1].lower() in ('.asdf', '.h5'):
                tmpst = parse_asdf_traces(filename)
            else:
                tmpst = read(filename, fsize=False)
        except Exception:
            logger.warning(
                f'{filename}: Unable to read file as a trace: skipping')
            continue
        for trace in tmpst.traces:
            # only use the station specified by the command line option
            # "--station", if any
            if (config.options.station is not None and
                    trace.stats.station != config.options.station):
                continue
            st.append(trace)
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
    Read trace files

    :return: Traces
    :rtype: :class:`obspy.core.stream.Stream`
    """
    logger.info('Reading traces...')
    stream = _read_trace_files()
    _read_event_from_traces(stream)
    logger.info('Reading traces: done')
    logger.info('---------------------------------------------------')
    if len(stream) == 0:
        logger.error('No trace loaded')
        ssp_exit(1)
    return stream
