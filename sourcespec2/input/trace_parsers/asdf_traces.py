# -*- coding: utf-8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Read traces from ASDF files.

:copyright:
    2012-2024 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
import logging
import json
from obspy.core import Stream
from ...setup import ssp_exit
logger = logging.getLogger(__name__.rsplit('.', maxsplit=1)[-1])


def _parse_asdf_trace_headers(ds, stream, nw_stat_codes, tag):
    """
    Parse ASDF trace headers.

    :param ds: ASDF dataset
    :type ds: :class:`pyasdf.ASDFDataSet`
    :param stream: ObsPy Stream object
    :type stream: :class:`obspy.core.stream.Stream`
    :param nw_stat_codes: Network and station codes
    :type nw_stat_codes: list
    :param tag: waveform tag in ASDF file
    :type tag: str
    """
    # pylint: disable=import-outside-toplevel
    from pyasdf.utils import AuxiliaryDataContainer
    for nw_stat_code, tr in zip(nw_stat_codes, stream.traces):
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
    stream = Stream()
    try:
        ds = pyasdf.ASDFDataSet(asdf_file, mode='r')
    except OSError:
        logger.warning(f'Unable to read ASDF file: {asdf_file}')
        return stream
    if not ds.waveforms:
        return stream
    nw_stat_codes = []
    for nw_stat_code in ds.waveforms.list():
        wf_tags = ds.waveforms[nw_stat_code].get_waveform_tags()
        # If tag is not specified, take first available tag
        if not tag:
            # Maybe this should be logged
            tag = wf_tags[0]
        if tag in wf_tags:
            station_st = ds.waveforms[nw_stat_code][tag]
            stream.extend(station_st)
            nw_stat_codes.extend([nw_stat_code] * len(station_st))
    # Try reading trace headers if present
    if 'TraceHeaders' in ds.auxiliary_data:
        header_key = 'TraceHeaders'
    elif 'Headers' in ds.auxiliary_data:
        header_key = 'Headers'
    else:
        header_key = None
    if read_headers and header_key:
        _parse_asdf_trace_headers(ds, stream, nw_stat_codes, tag)
    ds._close()
    return stream