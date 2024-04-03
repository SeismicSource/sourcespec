.. _spectral_file_formats:

#####################
Spectral File Formats
#####################

Spectra produced by SourceSpec can be optionally saved in a `HDF5`_ file
named ``EVID.spectra.hdf5`` (where ``EVID`` is the event ID; see the
``save_spectra`` option in :ref:`configuration_file:Configuration File`).

Additionally, spectral residuals are always saved to a file named
``EVID.residuals.hdf5``, and average residuals produced by ``source_residuals``
are saved to a file named ``residual_mean.hdf5``.

All these files use an `HDF5`_ based format, detailed in the following
section. They can be read using the function :func:`spectrum.read_spectra`,
which returns a :class:`spectrum.SpectrumStream` object, consisting of
:class:`spectrum.Spectrum` objects.

Example code:

.. code-block:: python

    >>> from sourcespec.spectrum import read_spectra
    >>> specst = read_spectra('38471103.spectra.hdf5')
    >>> print(specst)
    SpectrumStream with 439 Spectrum objects:
    CI.CCA..HHE | 200 samples, 0.2-39.9 Hz | 0.2 Hz sample interval | 59 samples logspaced, 0.20-41.70 Hz | 0.04 log10([Hz]) sample interval logspaced
    CI.CCA..HHH | 200 samples, 0.2-39.9 Hz | 0.2 Hz sample interval | 59 samples logspaced, 0.20-41.70 Hz | 0.04 log10([Hz]) sample interval logspaced
    CI.CCA..HHN | 200 samples, 0.2-39.9 Hz | 0.2 Hz sample interval | 59 samples logspaced, 0.20-41.70 Hz | 0.04 log10([Hz]) sample interval logspaced
    CI.CCA..HHS | 200 samples, 0.2-39.9 Hz | 0.2 Hz sample interval | 59 samples logspaced, 0.20-41.70 Hz | 0.04 log10([Hz]) sample interval logspaced
    CI.CCA..HHZ | 200 samples, 0.2-39.9 Hz | 0.2 Hz sample interval | 59 samples logspaced, 0.20-41.70 Hz | 0.04 log10([Hz]) sample interval logspaced
    CI.CCA..HHh | 200 samples, 0.2-39.9 Hz | 0.2 Hz sample interval | 59 samples logspaced, 0.20-41.70 Hz | 0.04 log10([Hz]) sample interval logspaced
    CI.CCA..HHs | 200 samples, 0.2-39.9 Hz | 0.2 Hz sample interval | 59 samples logspaced, 0.20-41.70 Hz | 0.04 log10([Hz]) sample interval logspaced
    ...
    >>> print(specst[0])
    Spectrum CI.CCA..HHE | 200 samples, 0.2-39.9 Hz | 0.2 Hz sample interval | 59 samples logspaced, 0.20-41.70 Hz | 0.04 log10([Hz]) sample interval logspaced
    >>> print(specst[0].stats)
    {'delta': 0.1996007984031936,
     'npts': 200,
     'delta_logspaced': 0.040000000000000036,
     'npts_logspaced': 59,
     'station': 'CCA',
     'network': 'CI',
     'location': '',
     'channel': 'HHE',
     'azimuth': 198.9047069007039,
     'coeff': 1282817000215832.2,
     'coords': {'elevation': 0.71,
     'latitude': 35.15251922607422,
     'longitude': -118.01648712158203},
     ...


Additionally, :class:`spectrum.SpectrumStream` and :class:`spectrum.Spectrum`
objects can be saved to a text file format, with one ``Spectrum`` per file.
The details of this format are explained below.

Here is an example of how to read an HDF5 file and save it to a text file:

.. code-block:: python

    >>> from sourcespec.spectrum import read_spectra
    >>> specst = read_spectra('38471103.spectra.hdf5')
    >>> specst.write('38471103.spectra.txt', format='TEXT')
    >>> from os import listdir
    >>> print('\n'.join(sorted(listdir('.'))))
    ...
    38471103.spectra_0000.txt
    38471103.spectra_0001.txt
    38471103.spectra_0002.txt
    38471103.spectra_0003.txt
    38471103.spectra_0004.txt
    38471103.spectra_0005.txt
    38471103.spectra_0006.txt
    38471103.spectra_0007.txt
    38471103.spectra_0008.txt
    38471103.spectra_0009.txt
    38471103.spectra_0010.txt
    ...
    >>> specst0 = read_spectra('38471103.spectra_0000.txt', format='TEXT')
    >>> print(specst0)
    SpectrumStream with 1 Spectrum objects:
    CI.CCA..HHE | 200 samples, 0.2-39.9 Hz | 0.2 Hz sample interval | 59 samples logspaced, 0.20-41.70 Hz | 0.04 log10([Hz]) sample interval logspaced


HDF5 File Format
----------------
In the HDF5 file format, all the spectra are stored in a group named
``spectra``. This will allow for storing additional data types in the future.
Within the ``spectra`` group, each :class:`spectrum.Spectrum` object is stored
in a `group`_ named ``spectrum_NNNNN_NET.STA.LOC.CHAN``, where ``NNNN`` is the
index of the spectrum in the original :class:`spectrum.SpectrumStream` object.
For each group, metadata is stored in the `attributes`_ section, and data is
stored into 6 `datasets`_, as illustrated below:

.. code-block:: none

    HDF5 ──> attributes
      └── spectra
          ├── spectrum_NNNNN_NET.STA.LOC.CHAN ──> attributes
          |   ├── data
          |   ├── data_logspaced
          |   ├── data_mag
          |   ├── data_mag_logspaced
          |   ├── freq
          |   └── freq_logspaced
          ├── spectrum_NNNNN_NET.STA.LOC.CHAN ──> attributes
          ...

The mandatory metadata fields are:

- ``network``: the network code;
- ``station``: the station code;
- ``location``: the location code;
- ``channel``: the channel code;
- ``delta``: the sample interval for linearly spaced frequencies;
- ``npts``: the number of samples for linearly spaced frequencies and data;
- ``delta_logspaced``: the sample interval for logspaced frequencies (set to 1
  if logspaced frequencies are not used);
- ``npts_logspaced``: the number of samples for logspaced frequencies and data
  (set to 0 if logspaced frequencies are not used).

Other metadata fields might be present. Dictionary-like metadata fields are
stored as `YAML`_ strings.

The 6 datasets are:

- ``data`` (mandatory): spectral amplitude;
- ``data_logspaced`` (optional): spectral amplitude for logspaced frequencies;
- ``data_mag`` (optional): spectral amplitude in magnitude units;
- ``data_mag_logspaced`` (optional): spectral amplitude in magnitude units for
  logspaced frequencies;
- ``freq`` (mandatory): linearly spaced frequencies;
- ``freq_logspaced`` (optional): logspaced frequencies.

Example code for reading a SourceSpec HDF5 file using ``h5py``:

.. code-block:: python

    >>> import h5py
    >>> fp = h5py.File('38471103.spectra.hdf5', 'r')
    >>> spectra = fp['spectra']
    >>> print('\n'.join(spectra.keys()))
    spectrum_00000_CI.CCA..HHE
    spectrum_00001_CI.CCA..HHH
    spectrum_00002_CI.CCA..HHN
    spectrum_00003_CI.CCA..HHS
    ...
    >>> spec = spectra['spectrum_00000_CI.CCA..HHE']
    >>> print(spec.attrs['network'])
    CI
    >>> print(spec.attrs['station'])
    CCA
    >>> print(spec.attrs['channel'])
    HHE
    >>> print(spec.attrs['coords'])
    {'elevation': 0.959, 'latitude': 35.34149169921875, 'longitude': -116.87464141845703}
    >>> print(spec['freq'][...])
    [ 0.1996008   0.3992016   0.5988024   0.79840319  0.99800399  1.19760479
      1.39720559  1.59680639  1.79640719  1.99600798  2.19560878  2.39520958
    ...
    >>> print(spec['data'][...])
    [2.49035173e+15 1.37636948e+15 1.44746675e+15 1.63566457e+15
     6.93049162e+14 1.01315194e+15 9.36761128e+14 8.00776096e+14
    ...

TEXT File Format
----------------
The TEXT file format is not used internally by SourceSpec, but it can be useful
to convert the HDF5 files to a more human-readable format (see above for an
example on reading an HDF5 file and converting it to TEXT format).

The format is structured as follows:

- A header section in `YAML`_ format. The header is identified by the two
  lines ``# %BEGIN STATS YAML`` and ``# %END STATS YAML``. Each line in the
  header starts with a ``#`` character, which should be removed when
  using a YAML parser.
- One or two data sections, each with three columns: ``frequency (Hz)``,
  ``data``, ``data_mag`` (if ``data_mag`` is not present, it is set to ``nan``):

  - linearly spaced data, between ``# %BEGIN LINSPACED DATA`` and
    ``# %END LINSPACED DATA``;
  - (optional) logspaced data, between ``# %BEGIN LOGSPACED DATA`` and
    ``# %END LOGSPACED DATA``.

Here's an example TEXT file (with ellipses for brevity):

.. code-block::

    # %SOURCESPEC TEXT SPECTRUM FORMAT 1.0
    # %BEGIN STATS YAML
    # delta: 0.1996007984031936
    # npts: 200
    # delta_logspaced: 0.040000000000000036
    # npts_logspaced: 59
    # station: CCA
    # network: CI
    # location: ''
    # channel: HHE
    # azimuth: 198.9047069007039
    # coeff: 1282817000215832.2
    # coords:
    #   elevation: 0.71
    #   latitude: 35.15251922607422
    #   longitude: -118.01648712158203
    ...
    # %END STATS YAML
    # %BEGIN LINSPACED DATA
    # frequency(Hz) data data_mag
    0.199601 60680538429002.429688 3.122033
    0.399202 111489590392894.000000 3.298156
    0.598802 109341550136192.765625 3.292523
    0.798403 46532921128263.000000 3.045174
    0.998004 69772021841544.039062 3.162454
    ...
    # %END LINSPACED DATA
    # %BEGIN LOGSPACED DATA
    # frequency_logspaced(Hz) data_logspaced data_mag_logspaced
    0.199601 60680538429002.437500 3.122033
    0.218858 65978522113754.125000 3.146268
    0.239973 71686845260565.812500 3.170293
    0.263125 77968569379397.437500 3.194613
    0.288511 84857989097347.140625 3.219128
    ...
    # %END LOGSPACED DATA


.. _HDF5: https://www.hdfgroup.org/
.. _group: https://docs.h5py.org/en/stable/high/group.html
.. _attributes: https://docs.h5py.org/en/stable/high/attr.html
.. _datasets: https://docs.h5py.org/en/stable/high/dataset.html
.. _YAML: http://yaml.org/
