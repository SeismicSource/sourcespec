# -*- coding: utf-8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
QuakeML output for source_spec.

:copyright:
    2016-2026 Claudio Satriano <satriano@ipgp.fr>

:license:
    CeCILL Free Software License Agreement, Version 2.1
    (http://www.cecill.info/index.en.html)
"""
import os
import logging
import warnings
from obspy import read_events, UTCDateTime
from obspy.core import AttribDict
from obspy.core.event import (CreationInfo, FocalMechanism, Magnitude,
                              MomentTensor, QuantityError, ResourceIdentifier,
                              StationMagnitude, StationMagnitudeContribution,
                              WaveformStreamID)
from .config import config
from ._version import get_versions
logger = logging.getLogger(__name__.rsplit('.', maxsplit=1)[-1])


def _to_camel_case(snake_str):
    # source: http://stackoverflow.com/q/19053707
    components = snake_str.split('_')
    # We capitalize the first letter of each component except the first one
    # with the 'title' method and join them together.
    return components[0] + ''.join(x.title() for x in components[1:])


class SSPExtra(AttribDict):
    """Container for custom tags."""

    def __setattr__(self, key, value):
        """Set class attribute."""
        key = _to_camel_case(key)
        return super(AttribDict, self).__setattr__(key, value)

    def __getattr__(self, key, default=None):
        """Get class attribute."""
        key = _to_camel_case(key)
        return super(AttribDict, self).__getattr__(key, default)


class SSPContainerTag(AttribDict):
    """Container for nested custom tags."""

    def __init__(self):
        super().__init__()
        self.namespace = 'https://sourcespec.seismicsource.org'
        self.value = SSPExtra()


class SSPTag(AttribDict):
    """Custom tag object."""

    def __init__(self, value=None):
        super().__init__()
        self.namespace = 'https://sourcespec.seismicsource.org'
        self.value = value


def write_qml(sspec_output):
    """
    Write QuakeML output.

    :param sspec_output: Output from spectral inversion.
    :type sspec_output: :class:`~sourcespec.ssp_data_types.SourceSpecOutput`
    """
    if not config.options.qml_file:
        config.qml_file_out = None
        return
    qml_file = config.options.qml_file
    cat = read_events(qml_file)
    evid = config.event.event_id
    try:
        ev = [e for e in cat if evid in str(e.resource_id)][0]
    except Exception:
        logger.warning(
            f'Unable to find evid "{evid}" in QuakeML file. '
            'QuakeML output will not be written.')
        config.qml_file_out = None
        return

    origin = ev.preferred_origin()
    if origin is None:
        origin = ev.origins[0]
    origin_id = origin.resource_id
    origin_id_strip = origin_id.id.split('/')[-1]
    origin_id_strip = origin_id_strip.replace(
        config.smi_strip_from_origin_id, '')

    # Common parameters
    ssp_version = get_versions()['version']
    method_id = f'{config.smi_base}/sourcespec/{ssp_version}'
    cr_info = CreationInfo()
    if config.agency_short_name is not None:
        cr_info.agency_id = config.agency_short_name
    if config.author_name is not None:
        cr_info.author = config.author_name
    cr_info.creation_time = UTCDateTime()

    summary_parameters = sspec_output.reference_summary_parameters()
    stationpar = sspec_output.station_parameters

    # Magnitude
    Mw_summary = summary_parameters['Mw']
    mag = Magnitude()
    _id = config.smi_magnitude_template.replace('$SMI_BASE', config.smi_base)
    _id = _id.replace('$ORIGIN_ID', origin_id_strip)
    mag.resource_id = ResourceIdentifier(id=_id)
    mag.method_id = ResourceIdentifier(id=method_id)
    mag.origin_id = origin_id
    mag.magnitude_type = 'Mw'
    mag.mag = Mw_summary.value
    mag_err = QuantityError()
    if Mw_summary.uncertainty is not None:
        mag_err.uncertainty = Mw_summary.uncertainty
    elif Mw_summary.lower_uncertainty is not None:
        mag_err.lower_uncertainty = Mw_summary.lower_uncertainty
        mag_err.upper_uncertainty = Mw_summary.upper_uncertainty
    mag_err.confidence_level = Mw_summary.confidence_level
    mag.mag_errors = mag_err
    mag.station_count = len(list(stationpar.keys()))
    mag.evaluation_mode = 'automatic'
    mag.creation_info = cr_info

    # Seismic moment -- It has to be stored in a MomentTensor object
    # which, in turn, is part of a FocalMechanism object
    Mo_summary = summary_parameters['Mo']
    mt = MomentTensor()
    _id = config.smi_moment_tensor_template.replace(
        '$SMI_BASE', config.smi_base)
    _id = _id.replace('$ORIGIN_ID', origin_id_strip)
    mt.resource_id = ResourceIdentifier(id=_id)
    mt.derived_origin_id = origin_id
    mt.moment_magnitude_id = mag.resource_id
    mt.scalar_moment = Mo_summary.value
    mt_err = QuantityError()
    mt_err.lower_uncertainty = Mo_summary.lower_uncertainty
    mt_err.upper_uncertainty = Mo_summary.upper_uncertainty
    mt_err.confidence_level = Mo_summary.confidence_level
    mt.scalar_moment_errors = mt_err
    mt.method_id = method_id
    mt.creation_info = cr_info
    # And here is the FocalMechanism object
    fm = FocalMechanism()
    _id = config.smi_focal_mechanism_template.replace(
        '$SMI_BASE', config.smi_base)
    _id = _id.replace('$ORIGIN_ID', origin_id_strip)
    fm.resource_id = ResourceIdentifier(id=_id)
    fm.triggering_origin_id = origin_id
    fm.method_id = ResourceIdentifier(id=method_id)
    fm.moment_tensor = mt
    fm.creation_info = cr_info
    ev.focal_mechanisms.append(fm)

    # Station magnitudes
    for statId in sorted(stationpar.keys()):
        par = stationpar[statId]
        if par.ignored:
            continue
        st_mag = StationMagnitude()
        _id = config.smi_station_magnitude_template.replace(
            '$SMI_MAGNITUDE_TEMPLATE', config.smi_magnitude_template)
        _id = _id.replace('$ORIGIN_ID', origin_id_strip)
        _id = _id.replace('$SMI_BASE', config.smi_base)
        _id = _id.replace('$WAVEFORM_ID', statId)
        st_mag.resource_id = ResourceIdentifier(id=_id)
        st_mag.origin_id = origin_id
        st_mag.mag = par.Mw.value
        st_mag.station_magnitude_type = 'Mw'
        st_mag.method_id = mag.method_id
        st_mag.creation_info = cr_info
        st_mag.waveform_id = WaveformStreamID(seed_string=statId)
        st_mag.extra = SSPExtra()
        st_mag.extra.moment = SSPTag(par.Mo.value)
        st_mag.extra.corner_frequency = SSPTag(par.fc.value)
        st_mag.extra.t_star = SSPTag(par.t_star.value)
        ev.station_magnitudes.append(st_mag)
        st_mag_contrib = StationMagnitudeContribution()
        st_mag_contrib.station_magnitude_id = st_mag.resource_id
        mag.station_magnitude_contributions.append(st_mag_contrib)
    ev.magnitudes.append(mag)

    # Write other summary parameters as custom tags
    ev.extra = SSPExtra()
    ev.extra.corner_frequency = _summary_parameter_tag(
        summary_parameters['fc'])
    ev.extra.t_star = _summary_parameter_tag(summary_parameters['t_star'])
    ev.extra.source_radius = _summary_parameter_tag(
        summary_parameters['radius'])
    ev.extra.static_stress_drop = _summary_parameter_tag(
        summary_parameters['ssd'])
    ev.extra.radiated_energy = _summary_parameter_tag(
        summary_parameters['Er'])
    ev.extra.apparent_stress = _summary_parameter_tag(
        summary_parameters['sigma_a'])
    if config.set_preferred_magnitude:
        ev.preferred_magnitude_id = mag.resource_id.id

    qml_file_out = os.path.join(config.options.outdir, f'{evid}.xml')
    with warnings.catch_warnings(record=True) as warns:
        ev.write(qml_file_out, format='QUAKEML')
        for w in warns:
            message = str(w.message)
            # Ignore a couple of harmless warnings
            if 'trimmed mean' in message or 'ITAPER' in message:
                continue
            logger.warning(f'Warning while writing QuakeML: {message}')
    logger.info(f'QuakeML file written to: {qml_file_out}')
    config.qml_file_out = qml_file_out


def _summary_parameter_tag(param):
    tag = SSPContainerTag()
    tag.value.value = SSPTag(param.value)
    if param.uncertainty is not None:
        tag.value.uncertainty = SSPTag(param.uncertainty)
    elif param.lower_uncertainty is not None:
        tag.value.lower_uncertainty = SSPTag(param.lower_uncertainty)
        tag.value.upper_uncertainty = SSPTag(param.upper_uncertainty)
    tag.value.confidence_level = SSPTag(param.confidence_level)
    return tag
