# -*- coding: utf-8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Define Spectrum and SpectrumStream classes,
similar to ObsPy's Trace and Stream.

Provides the high-level function read_spectra() to read
SpectrumStream objects from HDF5 or TEXT files.

:copyright:
    2012-2024 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
import os
import copy
import fnmatch
import warnings
import math
import h5py
import yaml
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
        stats = _normalize_metadata_object(self.stats)
        for attr, value in stats.items():
            # convert dictionaries to strings using YAML
            if hasattr(value, 'items'):
                value = yaml.dump(
                    value,
                    Dumper=_HDF5HeaderDumper,
                    default_flow_style=True
                ).replace('\n', '')
            # if value is a list-like,
            # check if all elements are of the same type
            elif hasattr(value, '__iter__') and len(value) > 0:
                type0 = type(value[0])
                if not all(isinstance(v, type0) for v in value):
                    raise ValueError(
                        f'All values of attribute "{attr}" must be of the '
                        'same type'
                    )
            # ignore unsupported types
            elif not isinstance(value, (int, float, bool, str)):
                warnings.warn(
                    f'Attribute "{attr}" is not a supported type '
                    f'({type(value)}) and will be ignored'
                )
                continue
            group.attrs[attr] = value
        group.create_dataset('freq', data=self.freq)
        group.create_dataset('data', data=self.data)
        group.create_dataset('data_mag', data=self.data_mag)
        group.create_dataset('freq_logspaced', data=self.freq_logspaced)
        group.create_dataset('data_logspaced', data=self.data_logspaced)
        group.create_dataset(
            'data_mag_logspaced', data=self.data_mag_logspaced)

    def _write_text(self, filename):
        """
        Write the spectrum to a TEXT file.

        :param filename: The name of the file to write to.
        """
        with open(filename, 'w', encoding='utf-8') as fp:
            fp.write('# %SOURCESPEC TEXT SPECTRUM FORMAT 1.0\n')
            fp.write('# %BEGIN STATS YAML\n')
            stats = _normalize_metadata_object(self.stats)
            stats_str = yaml.safe_dump(
                stats,
                sort_keys=False
            ).rstrip()
            for line in stats_str.split('\n'):
                fp.write(f'# {line}\n')
            fp.write(
                '# %END STATS YAML\n'
                '# %BEGIN LINSPACED DATA\n'
                '# frequency(Hz) data data_mag\n'
            )
            if self.data_mag.size:
                data_mag = self.data_mag
            else:
                data_mag = np.ones_like(self.data) * np.nan
            for freq, data, data_mag in zip(self.freq, self.data, data_mag):
                fp.write(f'{freq:.6f} {data:.6f} {data_mag:.6f}\n')
            fp.write(
                '# %END LINSPACED DATA\n'
                '# %BEGIN LOGSPACED DATA\n'
                '# frequency_logspaced(Hz) data_logspaced data_mag_logspaced\n'
            )
            if self.data_mag_logspaced.size:
                data_mag_logspaced = self.data_mag_logspaced
            else:
                data_mag_logspaced = np.ones_like(self.data_logspaced) * np.nan
            for freq_logspaced, data_logspaced, data_mag_logspaced in zip(
                    self.freq_logspaced, self.data_logspaced,
                    data_mag_logspaced):
                fp.write(
                    f'{freq_logspaced:.6f} {data_logspaced:.6f} '
                    f'{data_mag_logspaced:.6f}\n'
                )
            fp.write('# %END LOGSPACED DATA\n')

    # pylint: disable=redefined-builtin
    def write(self, filename, format='HDF5', append=False):
        """
        Write the spectrum to a file.

        :param filename: The name of the file to write to.
        :param format: The format to use. One of 'HDF5' or 'TEXT'.
            Default is 'HDF5'.
        :param append: If True, append the spectrum to an existing file.
            Only valid for HDF5 format.
        """
        if format == 'HDF5':
            if append:
                fp = h5py.File(filename, 'a')
                fp.attrs['format'] = 'SourceSpec SpectrumStream HDF5'
                fp.attrs['version'] = '1.0'
                fp.attrs['url'] = 'https://sourcespec.seismicsource.org'
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
        elif format == 'TEXT':
            if append:
                raise ValueError('Cannot append to a TEXT file')
            self._write_text(filename)
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
        :param format: The format to use. One of 'HDF5' or 'TEXT'.
        """
        if format == 'HDF5':
            append = False
            for spectrum in self:
                spectrum.write(filename, format, append)
                append = True
        elif format == 'TEXT':
            if len(self) == 1:
                self[0].write(filename, format)
            else:
                for n, spectrum in enumerate(self):
                    _root, _ext = os.path.splitext(filename)
                    _filename = f'{_root}_{n:04d}{_ext}'
                    spectrum.write(_filename, format)
        else:
            raise ValueError(f'Unsupported format: {format}')


# ---- Reading/writing functions and helper functions ----

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


class _HDF5HeaderDumper(yaml.SafeDumper):
    """
    A YAML dumper used for writing the HDF5 header.
    """


def _default_yaml_representer(dumper, data):
    """
    Default YAML representer for unsupported types.

    :param dumper: The YAML dumper.
    :param data: The data to represent.
    :return: The YAML representation of the data.
    """
    return dumper.represent_scalar('tag:yaml.org,2002:str', str(data))


def _quoted_representer(dumper, data):
    """
    YAML representer for strings, with quotes.

    :param dumper: The YAML dumper.
    :param data: The data to represent.
    :return: The YAML representation of the data.
    """
    return dumper.represent_scalar('tag:yaml.org,2002:str', data, style="'")


# register the representers
yaml.representer.SafeRepresenter.add_representer(
    None, _default_yaml_representer)
_HDF5HeaderDumper.add_representer(str, _quoted_representer)
_HDF5HeaderDumper.add_representer(None, _default_yaml_representer)


def _normalize_metadata_object(obj):
    """
    Normalize a metadata object to use:
    - dictionaries instead of custom objects;
    - standard floats instead of numpy floats;
    - standard booleans instead of numpy booleans.

    All the other types are left unchanged.

    :param obj: The object to normalize.
    :return: A dictionary or the original value.
    """
    if hasattr(obj, 'items'):
        return {
            key: _normalize_metadata_object(val) for key, val in obj.items()
        }
    if isinstance(obj, (np.float64, np.float32)):
        obj = float(obj)
    if isinstance(obj, np.bool_):
        obj = bool(obj)
    return obj


def _read_spectrum_from_hdf5_group(group):
    """
    Read a Spectrum object from an HDF5 group.

    :param group: The HDF5 group to read from.
    :return: The Spectrum object.
    """
    spectrum = Spectrum()
    for attr, value in group.attrs.items():
        # convert strings back to dictionaries, using YAML
        if (
            isinstance(value, str)
            and value.startswith('{') and value.endswith('}')
        ):
            try:
                value = yaml.safe_load(value)
            except yaml.YAMLError:
                warnings.warn(
                    f'Attribute "{attr}" is not a supported type and will be '
                    'ignored'
                )
        spectrum.stats[attr] = value
    spectrum.freq = group['freq']
    spectrum.data = group['data']
    spectrum.data_mag = group['data_mag']
    spectrum.freq_logspaced = group['freq_logspaced']
    spectrum.data_logspaced = group['data_logspaced']
    spectrum.data_mag_logspaced = group['data_mag_logspaced']
    return spectrum


def _read_stats_from_text_lines(lines):
    """
    Read the stats block from a TEXT file.

    :param lines: The lines to read from.
    :return: The stats block.
    """
    _stats_lines = []
    _stats_block = False
    for line in lines:
        if line.startswith('# %BEGIN STATS YAML'):
            _stats_block = True
        elif line.startswith('# %END STATS YAML'):
            _stats_block = False
            break
        elif _stats_block:
            _stats_lines.append(line[2:])
    return yaml.safe_load(''.join(_stats_lines))


def _read_data_block_from_text_lines(lines, start_string, end_string):
    """
    Read a data block from a TEXT file.

    :param lines: The lines to read from.
    :param start_string: The string marking the start of the data block.
    :param end_string: The string marking the end of the data block.
    :return: The data block.
    """
    _data_lines = []
    _data_block = False
    for line in lines:
        if line.startswith(start_string):
            _data_block = True
        elif line.startswith(end_string):
            _data_block = False
            break
        elif _data_block:
            if line.startswith('#'):
                continue
            _data_lines.append(line.split())
    freq, data, data_mag = np.array(_data_lines, dtype=float).T
    return freq, data, data_mag


def _read_spectrum_from_text_file(filename):
    """
    Read a Spectrum object from a TEXT file.

    :param filename: The name of the file to read from.
    :return: The Spectrum object.
    """
    with open(filename, 'r', encoding='utf-8') as fp:
        lines = fp.readlines()
        stats = _read_stats_from_text_lines(lines)
        freq, data, data_mag = _read_data_block_from_text_lines(
            lines, '# %BEGIN LINSPACED DATA', '# %END LINSPACED DATA'
        )
        freq_logspaced, data_logspaced, data_mag_logspaced =\
            _read_data_block_from_text_lines(
                lines, '# %BEGIN LOGSPACED DATA', '# %END LOGSPACED DATA'
            )
    spectrum = Spectrum()
    spectrum.stats = AttributeDict(stats)
    spectrum.freq = freq
    spectrum.data = data
    if not np.isnan(data_mag).all():
        spectrum.data_mag = data_mag
    spectrum.freq_logspaced = freq_logspaced
    spectrum.data_logspaced = data_logspaced
    if not np.isnan(data_mag_logspaced).all():
        spectrum.data_mag_logspaced = data_mag_logspaced
    return spectrum


# pylint: disable=redefined-builtin
def read_spectra(filename, format='HDF5'):
    """
    Read a SpectrumStream from a file.

    :param filename: The name of the file to read from.
    :param format: The format to use. One of 'HDF5' or 'TEXT'.
        Default is 'HDF5'.
    :return: The SpectrumStream object.
    """
    if format == 'HDF5':
        with h5py.File(filename, 'r') as fp:
            spectra = SpectrumStream()
            for group in fp.values():
                spectra.append(_read_spectrum_from_hdf5_group(group))
        return spectra
    elif format == 'TEXT':
        spectra = SpectrumStream()
        spectra.append(_read_spectrum_from_text_file(filename))
        return spectra
    else:
        raise ValueError(f'Unsupported format: {format}')
