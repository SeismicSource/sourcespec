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


class SpectrumStream(list):
    """
    A class to handle a collection of amplitude spectra.
    """
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
