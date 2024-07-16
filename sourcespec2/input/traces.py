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
from obspy import read
from obspy.core import Stream
from ..setup import config, ssp_exit
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


def read_traces():
    """
    Read trace files

    :return: Traces
    :rtype: :class:`obspy.core.stream.Stream`
    """
    logger.info('Reading traces...')
    st = _read_trace_files()
    logger.info('Reading traces: done')
    logger.info('---------------------------------------------------')
    if len(st) == 0:
        logger.error('No trace loaded')
        ssp_exit(1)
    return st
