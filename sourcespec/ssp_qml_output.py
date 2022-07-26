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


def write_qml(config, sourcepar):
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

    means = sourcepar.means_weight
    errors = sourcepar.errors_weight
    stationpar = sourcepar.station_parameters

    # Magnitude
    mag = Magnitude()
    _id = config.smi_magnitude_template.replace('$SMI_BASE', config.smi_base)
    _id = _id.replace('$ORIGIN_ID', origin_id_strip)
    mag.resource_id = ResourceIdentifier(id=_id)
    mag.method_id = ResourceIdentifier(id=method_id)
    mag.origin_id = origin_id
    mag.magnitude_type = 'Mw'
    mag.mag = means['Mw']
    mag_err = QuantityError()
    mag_err.uncertainty = errors['Mw']
    mag_err.confidence_level = 68.2
    mag.mag_errors = mag_err
    mag.station_count = len([_s for _s in stationpar.keys()])
    mag.evaluation_mode = 'automatic'
    mag.creation_info = cr_info

    # Seismic moment -- It has to be stored in a MomentTensor object
    # which, in turn, is part of a FocalMechanism object
    mt = MomentTensor()
    _id = config.smi_moment_tensor_template.replace(
        '$SMI_BASE', config.smi_base)
    _id = _id.replace('$ORIGIN_ID', origin_id_strip)
    mt.resource_id = ResourceIdentifier(id=_id)
    mt.derived_origin_id = origin_id
    mt.moment_magnitude_id = mag.resource_id
    mt.scalar_moment = means['Mo']
    mt_err = QuantityError()
    mt_err.lower_uncertainty = errors['Mo'][0]
    mt_err.upper_uncertainty = errors['Mo'][1]
    mt_err.confidence_level = 68.2
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
        seed_id = statId.split()[0]
        _id = config.smi_station_magnitude_template.replace(
            '$SMI_MAGNITUDE_TEMPLATE', config.smi_magnitude_template)
        _id = _id.replace('$ORIGIN_ID', origin_id_strip)
        _id = _id.replace('$SMI_BASE', config.smi_base)
        _id = _id.replace('$WAVEFORM_ID', seed_id)
        st_mag.resource_id = ResourceIdentifier(id=_id)
        st_mag.origin_id = origin_id
        st_mag.mag = par['Mw']
        st_mag.station_magnitude_type = 'Mw'
        st_mag.method_id = mag.method_id
        st_mag.creation_info = cr_info
        st_mag.waveform_id = WaveformStreamID(seed_string=seed_id)
        st_mag.extra = SSPExtra()
        st_mag.extra.moment = SSPTag(par['Mo'])
        st_mag.extra.corner_frequency = SSPTag(par['fc'])
        st_mag.extra.t_star = SSPTag(par['t_star'])
        ev.station_magnitudes.append(st_mag)
        st_mag_contrib = StationMagnitudeContribution()
        st_mag_contrib.station_magnitude_id = st_mag.resource_id
        mag.station_magnitude_contributions.append(st_mag_contrib)
    ev.magnitudes.append(mag)

    # Write other average parameters as custom tags
    ev.extra = SSPExtra()
    ev.extra.corner_frequency = SSPContainerTag()
    ev.extra.corner_frequency.value.value = SSPTag(means['fc'])
    ev.extra.corner_frequency.value.lower_uncertainty =\
        SSPTag(errors['fc'][0])
    ev.extra.corner_frequency.value.upper_uncertainty =\
        SSPTag(errors['fc'][1])
    ev.extra.corner_frequency.value.confidence_level = SSPTag(68.2)
    ev.extra.t_star = SSPContainerTag()
    ev.extra.t_star.value.value = SSPTag(means['t_star'])
    ev.extra.t_star.value.uncertainty = SSPTag(errors['t_star'])
    ev.extra.t_star.value.confidence_level = SSPTag(68.2)
    ev.extra.source_radius = SSPContainerTag()
    ev.extra.source_radius.value.value = SSPTag(means['ra'])
    ev.extra.source_radius.value.lower_uncertainty =\
        SSPTag(errors['ra'][0])
    ev.extra.source_radius.value.upper_uncertainty =\
        SSPTag(errors['ra'][1])
    ev.extra.source_radius.value.confidence_level = SSPTag(68.2)
    ev.extra.stress_drop = SSPContainerTag()
    ev.extra.stress_drop.value.value = SSPTag(means['bsd'])
    ev.extra.stress_drop.value.lower_uncertainty =\
        SSPTag(errors['bsd'][0])
    ev.extra.stress_drop.value.upper_uncertainty =\
        SSPTag(errors['bsd'][1])
    ev.extra.stress_drop.value.confidence_level = SSPTag(68.2)

    if config.set_preferred_magnitude:
        ev.preferred_magnitude_id = mag.resource_id.id

    qml_file_out = os.path.join(config.options.outdir, evid + '.xml')
    ev.write(qml_file_out, format='QUAKEML')
    logging.info('QuakeML file written to: ' + qml_file_out)
    config.qml_file_out = qml_file_out
