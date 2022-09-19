# -*- coding: utf-8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
QuakeML output for source_spec.

:copyright:
    2016-2022 Claudio Satriano <satriano@ipgp.fr>

:license:
    CeCILL Free Software License Agreement, Version 2.1
    (http://www.cecill.info/index.en.html)
"""
import os
import logging
from obspy import read_events, UTCDateTime
from obspy.core import AttribDict
from obspy.core.event import (CreationInfo, FocalMechanism, Magnitude,
                              MomentTensor, QuantityError, ResourceIdentifier,
                              StationMagnitude, StationMagnitudeContribution,
                              WaveformStreamID)
from sourcespec._version import get_versions


def _to_camel_case(snake_str):
    # source: http://stackoverflow.com/q/19053707
    components = snake_str.split('_')
    # We capitalize the first letter of each component except the first one
    # with the 'title' method and join them together.
    return components[0] + "".join(x.title() for x in components[1:])


class SSPExtra(AttribDict):
    """Container for custom tags."""

    def __setattr__(self, key, value):
        """Set class attribute."""
        key = _to_camel_case(key)
        return super(AttribDict, self).__setattr__(key, value)

    def __getattr__(self, key):
        """Get class attribute."""
        key = _to_camel_case(key)
        return self.__dict__[key]


class SSPContainerTag(AttribDict):
    """Container for nested custom tags."""

    def __init__(self):
        self.namespace = 'https://sourcespec.seismicsource.org'
        self.value = SSPExtra()


class SSPTag(AttribDict):
    """Custom tag object."""

    def __init__(self, value=None):
        self.namespace = 'https://sourcespec.seismicsource.org'
        self.value = value


def write_qml(config, sspec_output):
    if not config.options.qml_file:
        config.qml_file_out = None
        return
    qml_file = config.options.qml_file
    cat = read_events(qml_file)
    evid = config.hypo.evid
    try:
        ev = [e for e in cat if evid in str(e.resource_id)][0]
    except Exception:
        logging.warning('Unable to find evid "{}" in QuakeML file. '
                        'QuakeML output will not be written.'.format(evid))

    origin = ev.preferred_origin()
    if origin is None:
        origin = ev.origins[0]
    origin_id = origin.resource_id
    origin_id_strip = origin_id.id.split('/')[-1]
    origin_id_strip = origin_id_strip.replace(
        config.smi_strip_from_origin_id, '')

    # Common parameters
    ssp_version = get_versions()['version']
    method_id = config.smi_base + '/sourcespec/' + ssp_version
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
    mag.station_count = len([_s for _s in stationpar.keys()])
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
    fc_summary = summary_parameters['fc']
    ev.extra = SSPExtra()
    ev.extra.corner_frequency = SSPContainerTag()
    ev.extra.corner_frequency.value.value = SSPTag(fc_summary.value)
    ev.extra.corner_frequency.value.lower_uncertainty =\
        SSPTag(fc_summary.lower_uncertainty)
    ev.extra.corner_frequency.value.upper_uncertainty =\
        SSPTag(fc_summary.upper_uncertainty)
    ev.extra.corner_frequency.value.confidence_level =\
        SSPTag(fc_summary.confidence_level)
    t_star_summary = summary_parameters['t_star']
    ev.extra.t_star = SSPContainerTag()
    ev.extra.t_star.value.value = SSPTag(t_star_summary.value)
    if t_star_summary.uncertainty is not None:
        ev.extra.t_star.value.uncertainty =\
            SSPTag(t_star_summary.uncertainty)
    elif t_star_summary.lower_uncertainty is not None:
        ev.extra.t_star.value.lower_uncertainty =\
            SSPTag(t_star_summary.lower_uncertainty)
        ev.extra.t_star.value.upper_uncertainty =\
            SSPTag(t_star_summary.upper_uncertainty)
    ev.extra.t_star.value.confidence_level =\
        SSPTag(t_star_summary.confidence_level)
    radius_summary = summary_parameters['radius']
    ev.extra.source_radius = SSPContainerTag()
    ev.extra.source_radius.value.value = SSPTag(radius_summary.value)
    ev.extra.source_radius.value.lower_uncertainty =\
        SSPTag(radius_summary.lower_uncertainty)
    ev.extra.source_radius.value.upper_uncertainty =\
        SSPTag(radius_summary.upper_uncertainty)
    ev.extra.source_radius.value.confidence_level =\
        SSPTag(radius_summary.confidence_level)
    bsd_summary = summary_parameters['bsd']
    ev.extra.stress_drop = SSPContainerTag()
    ev.extra.stress_drop.value.value = SSPTag(bsd_summary.value)
    ev.extra.stress_drop.value.lower_uncertainty =\
        SSPTag(bsd_summary.lower_uncertainty)
    ev.extra.stress_drop.value.upper_uncertainty =\
        SSPTag(bsd_summary.upper_uncertainty)
    ev.extra.stress_drop.value.confidence_level =\
        SSPTag(bsd_summary.confidence_level)

    if config.set_preferred_magnitude:
        ev.preferred_magnitude_id = mag.resource_id.id

    qml_file_out = os.path.join(config.options.outdir, evid + '.xml')
    ev.write(qml_file_out, format='QUAKEML')
    logging.info('QuakeML file written to: ' + qml_file_out)
    config.qml_file_out = qml_file_out
