
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
SourceSpec event class and supporting classes.

:copyright:
    2023 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
import numpy as np
import contextlib
from obspy import UTCDateTime
from obspy.imaging.scripts import mopad


def _float(value):
    return None if value is None else float(value)


def _time(value):
    return None if value is None else UTCDateTime(value)


def _km_to_m(value):
    return None if value is None else value * 1e3


def _m_to_km(value):
    return None if value is None else value / 1e3


def _dyne_cm_to_N_m(value):
    return None if value is None else value * 1e-7


class SSPDepth(object):
    """
    SourceSpec depth class.
    """
    def __init__(self, value=None, units=None):
        self.value = _float(value)
        self.units = units

    def __str__(self):
        try:
            return(f'{self.value:.1f} {self.units}')
        except TypeError:
            return('Incomplete depth data')

    @property
    def value_in_m(self):
        """Return the depth value in m."""
        if self.value is None or self.units is None:
            return None
        if self.units == 'm':
            return self.value
        elif self.units == 'km':
            return _km_to_m(self.value)
        else:
            raise ValueError('depth units must be m or km')

    @property
    def value_in_km(self):
        """Return the depth value in km."""
        if self.value is None or self.units is None:
            return None
        if self.units == 'km':
            return self.value
        elif self.units == 'm':
            return _m_to_km(self.value)
        else:
            raise ValueError('depth units must be m or km')


class SSPHypocenter(object):
    """
    SourceSpec hypocenter class.
    """
    def __init__(self, longitude=None, latitude=None, depth=None,
                 origin_time=None):
        self.longitude = _float(longitude)
        self.latitude = _float(latitude)
        depth = {} if depth is None else depth
        if not hasattr(depth, 'keys'):
            raise ValueError('depth must be a dictionary-like object')
        self.depth = SSPDepth(**depth)
        self.origin_time = _time(origin_time)

    def __str__(self):
        try:
            return(
                f'Longitude: {self.longitude:.4f}°, '
                f'Latitude: {self.latitude:.4f}°, '
                f'Depth: {self.depth.value:.1f} {self.depth.units}, '
                f'Origin time: {self.origin_time}'
            )
        except TypeError:
            return('Incomplete hypocenter data')


class SSPMagnitude(object):
    """
    SourceSpec magnitude class.
    """
    def __init__(self, value=None, type=None):
        self.value = _float(value)
        self.type = type

    def __str__(self):
        try:
            return(
                f'Magnitude: {self.type} {self.value:.1f}'
            )
        except TypeError:
            return('Incomplete magnitude data')

    def from_scalar_moment(self, scalar_moment):
        """
        Initialize the magnitude object from a scalar moment.

        Args:
            scalar_moment (SSPScalarMoment): scalar moment object
        """
        if not isinstance(scalar_moment, SSPScalarMoment):
            raise TypeError('scalar_moment must be an SSPScalarMoment object')
        moment = scalar_moment.value
        if scalar_moment.units == 'dyne-cm':
            moment = _dyne_cm_to_N_m(moment)
        elif scalar_moment.units != 'N-m':
            raise ValueError('scalar_moment units must be N-m or dyne-cm')
        self.value = 2/3 * (np.log10(moment) - 9.1)
        self.type = 'Mw'


class SSPScalarMoment(object):
    """
    SourceSpec scalar moment class.
    """
    def __init__(self, value=None, units=None):
        self.value = _float(value)
        self.units = units
        self.to_N_m()

    def __str__(self):
        try:
            return(
                f'Scalar moment: {self.value:.1e} {self.units}'
            )
        except TypeError:
            return('Incomplete scalar moment data')

    def to_N_m(self):
        """Convert the scalar moment to N-m, if necessary."""
        if self.units is None or self.value is None:
            return
        if self.units not in ['N-m', 'dyne-cm']:
            raise ValueError('units must be N-m or dyne-cm')
        if self.units == 'dyne-cm':
            self.value = _dyne_cm_to_N_m(self.value)
            self.units = 'N-m'

    def from_moment_tensor(self, moment_tensor):
        """
        Initialize the scalar moment object from a moment tensor.

        Args:
            moment_tensor (SSPMomentTensor): moment tensor object
        """
        if not isinstance(moment_tensor, SSPMomentTensor):
            raise TypeError('moment_tensor must be an SSPMomentTensor object')
        _tensor_components = ('m_rr', 'm_tt', 'm_pp', 'm_rt', 'm_rp', 'm_tp')
        _tensor = [getattr(moment_tensor, cmp) for cmp in _tensor_components]
        self.value = mopad.MomentTensor(_tensor).get_moment()
        self.units = moment_tensor.units
        self.to_N_m()


class SSPFocalMechanism(object):
    """
    SourceSpec focal mechanism class.
    """
    def __init__(self, strike=None, dip=None, rake=None):
        self.strike = _float(strike)
        self.dip = _float(dip)
        self.rake = _float(rake)

    def __str__(self):
        try:
            return(
                f'Strike: {self.strike:.1f}°, '
                f'Dip: {self.dip:.1f}°, '
                f'Rake: {self.rake:.1f}°'
            )
        except TypeError:
            return('Incomplete focal mechanism data')

    def from_moment_tensor(self, moment_tensor):
        """
        Initialize the focal mechanism object from a moment tensor.

        Args:
            moment_tensor (SSPMomentTensor): moment tensor object
        """
        if not isinstance(moment_tensor, SSPMomentTensor):
            raise TypeError('moment_tensor must be an SSPMomentTensor object')
        _tensor_components = ('m_rr', 'm_tt', 'm_pp', 'm_rt', 'm_rp', 'm_tp')
        _tensor = [getattr(moment_tensor, cmp) for cmp in _tensor_components]
        _fps1, _fps2 = mopad.MomentTensor(_tensor).get_fps()
        self.strike, self.dip, self.rake = _fps1


class SSPMomentTensor(object):
    """
    SourceSpec moment tensor class.
    """
    def __init__(self, units=None, m_rr=None, m_tt=None, m_pp=None, m_rt=None,
                 m_rp=None, m_tp=None):
        self.units = units
        self.m_rr = _float(m_rr)
        self.m_tt = _float(m_tt)
        self.m_pp = _float(m_pp)
        self.m_rt = _float(m_rt)
        self.m_rp = _float(m_rp)
        self.m_tp = _float(m_tp)
        self.to_N_m()

    def __str__(self):
        try:
            return(
                f'units: {self.units}, '
                f'm_rr: {self.m_rr:.1e}, '
                f'm_tt: {self.m_tt:.1e}, '
                f'm_pp: {self.m_pp:.1e}, '
                f'm_rt: {self.m_rt:.1e}, '
                f'm_rp: {self.m_rp:.1e}, '
                f'm_tp: {self.m_tp:.1e}'
            )
        except TypeError:
            return('Incomplete moment tensor data')

    def to_N_m(self):
        """Convert the moment tensor to N-m, if necessary."""
        if self.units is None:
            return
        if self.units not in ['N-m', 'dyne-cm']:
            raise ValueError('units must be N-m or dyne-cm')
        if self.units == 'dyne-cm':
            for component in ('m_rr', 'm_tt', 'm_pp', 'm_rt', 'm_rp', 'm_tp'):
                setattr(
                    self, component, _dyne_cm_to_N_m(getattr(self, component)))
            self.units = 'N-m'


class SSPEvent(object):
    """
    SourceSpec event class.
    """
    def __init__(self, event_dict=None):
        """
        Initialize the event object.

        Args:
            file (str): path to the event file
        """
        self.event_id = None
        self.name = None
        self.hypocenter = SSPHypocenter()
        self.magnitude = SSPMagnitude()
        self.scalar_moment = SSPScalarMoment()
        self.focal_mechanism = SSPFocalMechanism()
        self.moment_tensor = SSPMomentTensor()
        if event_dict is not None:
            self.from_event_dict(event_dict)

    def __str__(self):
        return(
            f'Event ID: {self.event_id}\n'
            f'Name: {self.name}\n'
            f'Hypocenter:\n  {self.hypocenter}\n'
            f'Magnitude:\n  {self.magnitude}\n'
            f'Scalar moment:\n  {self.scalar_moment}\n'
            f'Focal mechanism:\n  {self.focal_mechanism}\n'
            f'Moment tensor:\n  {self.moment_tensor}'
        )

    def from_event_dict(self, event_dict):
        """
        Initialize the event object from a event dictionary.
        """
        try:
            self.event_id = event_dict['event_id']
        except KeyError as e:
            raise KeyError('event_id is required') from e
        self.name = event_dict.get('name')
        hypo_dict = event_dict.get('hypocenter', {})
        self.hypocenter = SSPHypocenter(**hypo_dict)
        mag_dict = event_dict.get('magnitude', {})
        self.magnitude = SSPMagnitude(**mag_dict)
        moment_dict = event_dict.get('scalar_moment', {})
        self.scalar_moment = SSPScalarMoment(**moment_dict)
        fm_dict = event_dict.get('focal_mechanism', {})
        self.focal_mechanism = SSPFocalMechanism(**fm_dict)
        mt_dict = event_dict.get('moment_tensor', {})
        self.moment_tensor = SSPMomentTensor(**mt_dict)
        # Try recomputing focal mechanism and scalar moment from moment tensor
        # and magnitude from scalar moment. Silently fail if any of these fail.
        with contextlib.suppress(Exception):
            self.focal_mechanism.from_moment_tensor(self.moment_tensor)
        with contextlib.suppress(Exception):
            self.scalar_moment.from_moment_tensor(self.moment_tensor)
        with contextlib.suppress(Exception):
            self.magnitude.from_scalar_moment(self.scalar_moment)
