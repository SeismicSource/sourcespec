# -*- coding: utf-8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Define Spectrum and SpectrumStream classes,
similar to ObsPy's Trace and Stream.

:copyright:
    2012-2024 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
import copy
import fnmatch
import warnings
import math
import h5py
import numpy as np


def signal_fft(signal, delta):
    """
    Compute the complex Fourier transform of a signal.

    :param signal: The signal to transform.
    :param delta: The sampling interval.
    :return: The Fourier transform and the frequency axis.
    """
    npts = len(signal)
    # if npts is even, we make it odd
    # so that we do not have a negative frequency in the last point
    # (see numpy.fft.rfft doc)
    if not npts % 2:
        npts -= 1
    # note that fft has the dimensions of the signal multiplied by time (delta)
    fft = np.fft.rfft(signal, n=npts) * delta
    fftfreq = np.fft.fftfreq(len(signal), d=delta)
    fftfreq = fftfreq[:fft.size]
    return fft, fftfreq


class AttributeDict(dict):
    """A dictionary that allows attribute-style access."""
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def __getattr__(self, name):
        return self[name]

    def __setattr__(self, name, value):
        self[name] = value

    def __dir__(self):
        return self.keys()

    def __copy__(self):
        new_dict = AttributeDict()
        for key, value in self.items():
            new_dict[key] = value
        return new_dict

    def __deepcopy__(self, memo):
        new_dict = AttributeDict()
        for key, value in self.items():
            new_dict[key] = copy.copy(value)
        return new_dict


def _n_significant_digits(x):
    """
    Helper function to compute the number of significant digits of a number.

    - If the number is greater than 1, the number of significant digits is
      zero.
    - If the number is less than 1, the number of significant digits is
      the number of digits after the decimal point.
    - If the number is zero, the number of significant digits is zero.
    """
    try:
        x = math.fabs(x)
    except TypeError as e:
        raise ValueError('x must be a number') from e
    return 0 if x == 0 or x > 1 else -int(math.floor(math.log10(x)))


class Spectrum():
    """
    A class to handle amplitude spectra.

    :param obspy_trace: An ObsPy Trace object to compute the spectrum from.
    """
    def __init__(self, obspy_trace=None):
        """
        Initialize the Spectrum object.

        :param obspy_trace: An ObsPy Trace object to compute the spectrum from.
        """
        self._data = np.array([], dtype=float)
        self._data_logspaced = np.array([], dtype=float)
        self._data_mag = np.array([], dtype=float)
        self._data_mag_logspaced = np.array([], dtype=float)
        self._freq = np.array([], dtype=float)
        self._freq_logspaced = np.array([], dtype=float)
        self.stats = AttributeDict()
        self.stats.delta = 1.
        self.stats.npts = 0
        self.stats.delta_logspaced = 1.
        self.stats.npts_logspaced = 0
        self.stats.station = ''
        self.stats.network = ''
        self.stats.location = ''
        self.stats.channel = ''
        if obspy_trace is not None:
            self.from_obspy_trace(obspy_trace)

    def __str__(self):
        delta = self.stats.delta
        ndigits = _n_significant_digits(delta)
        delta_str = f'{delta:.{ndigits}f}'
        fmin_str = f'{self.freq[0]:.{ndigits}f}' if self.freq.size else '...'
        fmax_str = f'{self.freq[-1]:.{ndigits}f}' if self.freq.size else '...'
        delta_logspaced = self.stats.delta_logspaced
        ndigits = _n_significant_digits(delta_logspaced)
        delta_logspaced_str = f'{delta_logspaced:.{ndigits}f}'
        fmin_logspaced_str =\
            f'{self.freq_logspaced[0]:.{ndigits}f}'\
            if self.freq_logspaced.size else '...'
        fmax_logspaced_str =\
            f'{self.freq_logspaced[-1]:.{ndigits}f}'\
            if self.freq_logspaced.size else '...'
        return (
            f'{self.id} | '
            f'{self.stats.npts} samples, {fmin_str}-{fmax_str} Hz | '
            f'{delta_str} Hz sample interval | '
            f'{self.stats.npts_logspaced} samples logspaced, '
            f'{fmin_logspaced_str}-{fmax_logspaced_str} Hz | '
            f'{delta_logspaced_str} log10([Hz]) sample interval logspaced '
        )

    def __repr__(self):
        return f'Spectrum {self}'

    def __gt__(self, other):
        return self.id > other.id

    def __lt__(self, other):
        return self.id < other.id

    def __ge__(self, other):
        return self.id >= other.id

    def __le__(self, other):
        return self.id <= other.id

    @property
    def id(self):
        """Return the id of the spectrum."""
        return (
            f'{self.stats.network}.{self.stats.station}.'
            f'{self.stats.location}.{self.stats.channel}'
        )

    @id.setter
    def id(self, value):
        """Set the id of the spectrum."""
        try:
            net, sta, loc, cha = value.split('.')
        except ValueError as e:
            raise ValueError(f'Not a valid SEED id: {value}') from e
        self.stats.network = net
        self.stats.station = sta
        self.stats.location = loc
        self.stats.channel = cha

    def get_id(self):
        """Return the id of the spectrum."""
        return self.id

    @property
    def data(self):
        """Return the array containing the amplitude spectrum."""
        return self._data

    @data.setter
    def data(self, value):
        """Set the array containing the amplitude spectrum."""
        if not isinstance(value, np.ndarray):
            value = np.array(value)
        self._data = value
        self.stats.npts = len(value)

    @property
    def data_mag(self):
        """
        Return the array containing the amplitude spectrum in mangitude units.
        """
        return self._data_mag

    @data_mag.setter
    def data_mag(self, value):
        """
        Set the array containing the amplitude spectrum in magnitude units.
        """
        if not isinstance(value, np.ndarray):
            value = np.array(value)
        if len(value) > 0 and len(value) != len(self.data):
            raise ValueError('data_mag must have the same length as data')
        self._data_mag = value

    @property
    def data_logspaced(self):
        """
        Return the array containing the amplitude spectrum in logspaced
        frequencies.
        """
        return self._data_logspaced

    @data_logspaced.setter
    def data_logspaced(self, value):
        """
        Set the array containing the amplitude spectrum in logspaced
        frequencies.
        """
        if not isinstance(value, np.ndarray):
            value = np.array(value)
        self._data_logspaced = value
        self.stats.npts_logspaced = len(value)

    @property
    def data_mag_logspaced(self):
        """
        Return the array containing the amplitude spectrum in logspaced
        frequencies in magnitude units.
        """
        return self._data_mag_logspaced

    @data_mag_logspaced.setter
    def data_mag_logspaced(self, value):
        """
        Set the array containing the amplitude spectrum in logspaced
        frequencies in magnitude units.
        """
        if not isinstance(value, np.ndarray):
            value = np.array(value)
        if len(value) > 0 and len(value) != len(self.data_logspaced):
            raise ValueError(
                'data_mag_logspaced must have the same length as '
                'data_logspaced'
            )
        self._data_mag_logspaced = value

    @property
    def freq(self):
        """Return the frequency axis of the spectrum."""
        return self._freq

    @freq.setter
    def freq(self, value):
        """Set the frequency axis of the spectrum."""
        if not isinstance(value, np.ndarray):
            value = np.array(value)
        if len(value) > 0:
            delta = np.diff(value)
            if not np.isclose(delta, delta[0], rtol=0.01).all():
                raise ValueError('Frequency axis must be evenly spaced')
            self.stats.delta = float(delta[0])
        else:
            self.stats.delta = 1.
        self._freq = value

    @property
    def freq_logspaced(self):
        """Return the logspaced frequency axis of the spectrum."""
        return self._freq_logspaced

    @freq_logspaced.setter
    def freq_logspaced(self, value):
        """Set the logspaced frequency axis of the spectrum."""
        if not isinstance(value, np.ndarray):
            value = np.array(value)
        if len(value) > 0:
            delta_logspaced = np.diff(np.log10(value))
            if not np.isclose(
                delta_logspaced, delta_logspaced[0], rtol=0.01
            ).all():
                raise ValueError(
                    'Logspaced frequency axis must be evenly spaced')
            self.stats.delta_logspaced = float(delta_logspaced[0])
        else:
            self.stats.delta_logspaced = 1.
        self._freq_logspaced = value

    def copy(self):
        """Return a copy of the spectrum."""
        spec_copy = Spectrum()
        spec_copy.stats = AttributeDict(self.stats)
        spec_copy.data = self.data.copy()
        spec_copy.data_mag = self.data_mag.copy()
        spec_copy.data_logspaced = self.data_logspaced.copy()
        spec_copy.data_mag_logspaced = self.data_mag_logspaced.copy()
        spec_copy.freq = self.freq.copy()
        spec_copy.freq_logspaced = self.freq_logspaced.copy()
        return spec_copy

    def slice(self, fmin, fmax, nearest_sample=True, pad=False,
              fill_value=None):
        """
        Slice the spectrum between fmin and fmax.

        :param fmin: Minimum frequency.
        :param fmax: Maximum frequency.
        :param nearest_sample: If True, the slice will include the nearest
            frequency to fmin and fmax.
        :param pad: If True, the slice will be padded with the value of
            fill_value until fmin and fmax are included.
        :param fill_value: The value to use for padding.

        :note: Only the linear spaced frequencies, data and data_mag are
            sliced. If the original spectrum contains logspaced frequencies,
            data, and data_mag, those are not preserved in the sliced spectrum.
        """
        freq = self.freq
        slice_condition = (freq >= fmin) & (freq <= fmax)
        if nearest_sample:
            # find the closest frequency to fmin:
            idx = (np.abs(freq - fmin)).argmin()
            slice_condition[idx] = True
            # find the closest frequency to fmax:
            idx = (np.abs(freq - fmax)).argmin()
            slice_condition[idx] = True
        freqs_slice = freq[slice_condition]
        data_slice = self.data[slice_condition]
        data_mag_slice =\
            self.data_mag[slice_condition] if self.data_mag.size\
            else np.array([])
        if pad:
            # add values to the slice until fmin and fmax are included
            while freqs_slice[0] > fmin:
                freqs_slice = np.insert(
                    freqs_slice, 0, freqs_slice[0] - self.stats.delta)
                data_slice = np.insert(data_slice, 0, fill_value)
                if data_mag_slice.size:
                    data_mag_slice = np.insert(data_mag_slice, 0, fill_value)
            while freqs_slice[-1] < fmax:
                freqs_slice = np.append(
                    freqs_slice, freqs_slice[-1] + self.stats.delta)
                data_slice = np.append(data_slice, fill_value)
                if data_mag_slice.size:
                    data_mag_slice = np.append(data_mag_slice, fill_value)
        spec_slice = Spectrum()
        spec_slice.stats = AttributeDict(self.stats)
        spec_slice.data = data_slice
        spec_slice.data_mag = data_mag_slice
        spec_slice.freq = freqs_slice
        return spec_slice

    def plot(self, **kwargs):
        """Plot the amplitude spectrum."""
        # pylint: disable=import-outside-toplevel
        import matplotlib.pyplot as plt
        plt.loglog(self.freq, self.data, **kwargs)
        plt.grid(True)
        plt.xlabel('Frequency (Hz)')
        plt.ylabel('Amplitude')
        plt.show()

    def from_obspy_trace(self, trace):
        """Compute the spectrum from an ObsPy Trace object."""
        # pylint: disable=import-outside-toplevel
        from obspy.core.trace import Trace
        if not isinstance(trace, Trace):
            raise TypeError('Only ObsPy Trace objects are supported')
        signal = trace.data
        delta = trace.stats.delta
        amp, freq = signal_fft(signal, delta)
        # remove DC component (freq=0)
        self.data = np.abs(amp)[1:]
        self.freq = freq[1:]
        # copy the trace metadata
        self.stats.station = trace.stats.station
        self.stats.network = trace.stats.network
        self.stats.location = trace.stats.location
        self.stats.channel = trace.stats.channel

    def _write_hdf5(self, group):
        """
        Write the spectrum to an HDF5 group.

        :param group: The HDF5 group to write to.
        """
        for attr, value in self.stats.items():
            # check if value is dict-like
            if hasattr(value, 'items'):
                # basic support for dict-like attributes, no nested dicts
                _keys = list(value.keys())
                _values = list(value.values())
                try:
                    group.attrs[f'_dict_{attr}_keys'] = _keys
                    group.attrs[f'_dict_{attr}_values'] = _values
                except TypeError:
                    warnings.warn(
                        f'Attribute "{attr}" is not a supported type and will '
                        'be ignored'
                    )
                    continue
            # check if it is a list-like
            elif hasattr(value, '__iter__'):
                if len(value) > 0:
                    type0 = type(value[0])
                    if not all(isinstance(v, type0) for v in value):
                        raise ValueError(
                            f'All values of attribute "{attr}" must be of the '
                            'same type'
                        )
                group.attrs[attr] = value
            # check if it is a number
            elif isinstance(value, (int, float)):
                group.attrs[attr] = value
            # ignore other types
            else:
                warnings.warn(
                    f'Attribute "{attr}" is not a supported type and will be '
                    'ignored'
                )
                continue
        group.create_dataset('freq', data=self.freq)
        group.create_dataset('data', data=self.data)
        group.create_dataset('data_mag', data=self.data_mag)
        group.create_dataset('freq_logspaced', data=self.freq_logspaced)
        group.create_dataset('data_logspaced', data=self.data_logspaced)
        group.create_dataset(
            'data_mag_logspaced', data=self.data_mag_logspaced)

    # pylint: disable=redefined-builtin
    def write(self, filename, format='HDF5', append=False):
        """
        Write the spectrum to a file.

        :param filename: The name of the file to write to.
        :param format: The format to use. Currently, only 'HDF5' is supported.
        :param append: If True, append the spectrum to an existing file.
        """
        if format == 'HDF5':
            if append:
                fp = h5py.File(filename, 'a')
                try:
                    lastgroup = sorted(fp.keys())[-1]
                    newgroup = f'spectrum_{int(lastgroup[-4:])+1:04d}'
                except IndexError:
                    newgroup = 'spectrum_0000'
            else:
                fp = h5py.File(filename, 'w')
                newgroup = 'spectrum_0000'
            self._write_hdf5(fp.create_group(newgroup))
            fp.close()
        else:
            raise ValueError(f'Unsupported format: {format}')


class SpectrumStream(list):
    """
    A class to handle a collection of amplitude spectra.
    """
    def __str__(self):
        return (
            f'SpectrumStream with {len(self)} Spectrum objects:\n'
            + '\n'.join(f'{s}' for s in self)
        )

    def __repr__(self):
        return self.__str__()

    def sort(self, reverse=False):
        """Sort the SpectrumStream in place."""
        super().sort(reverse=reverse)

    def append(self, spectrum):
        """Append a spectrum to the collection."""
        if not isinstance(spectrum, Spectrum):
            raise TypeError('Only Spectrum objects can be appended')
        super().append(spectrum)

    def select(self, **kwargs):
        """Select a subset of the SpectrumStream."""
        selected = SpectrumStream()
        for spectrum in self:
            for key, value in kwargs.items():
                if key == 'id':
                    stored_value = spectrum.get_id()
                else:
                    stored_value = getattr(spectrum.stats, key)
                if not fnmatch.fnmatch(stored_value, value):
                    break
            else:
                selected.append(spectrum)
        return selected

    # pylint: disable=redefined-builtin
    def write(self, filename, format='HDF5'):
        """
        Write the SpectrumStream to a file.

        :param filename: The name of the file to write to.
        :param format: The format to use. Currently, only 'HDF5' is supported.
        """
        if format == 'HDF5':
            append = False
            for spectrum in self:
                spectrum.write(filename, format, append)
                append = True
        else:
            raise ValueError(f'Unsupported format: {format}')


def _read_spectrum_from_hdf5_group(group):
    """
    Read a Spectrum object from an HDF5 group.

    :param group: The HDF5 group to read from.
    :return: The Spectrum object.
    """
    spectrum = Spectrum()
    for attr, value in group.attrs.items():
        # basic support for dict-like attributes, no nested dicts
        if attr.startswith('_dict_'):
            dict_attr = attr.split('_')[2]
            if dict_attr in spectrum.stats:
                # already processed
                continue
            _keys_attr = f'_dict_{dict_attr}_keys'
            _values_attr = f'_dict_{dict_attr}_values'
            if _keys_attr not in group.attrs or\
                    _values_attr not in group.attrs:
                continue
            keys = group.attrs[_keys_attr]
            values = group.attrs[_values_attr]
            spectrum.stats[dict_attr] = dict(zip(keys, values))
            continue
        spectrum.stats[attr] = value
    spectrum.freq = group['freq']
    spectrum.data = group['data']
    spectrum.data_mag = group['data_mag']
    spectrum.freq_logspaced = group['freq_logspaced']
    spectrum.data_logspaced = group['data_logspaced']
    spectrum.data_mag_logspaced = group['data_mag_logspaced']
    return spectrum


# pylint: disable=redefined-builtin
def read_spectra(filename, format='HDF5'):
    """
    Read a SpectrumStream from a file.

    :param filename: The name of the file to read from.
    :param format: The format to use. Currently, only 'HDF5' is supported.
    :return: The SpectrumStream object.
    """
    if format == 'HDF5':
        with h5py.File(filename, 'r') as fp:
            spectra = SpectrumStream()
            for group in fp.values():
                spectra.append(_read_spectrum_from_hdf5_group(group))
        return spectra
    else:
        raise ValueError(f'Unsupported format: {format}')
