.. _getting_started:

###############
Getting Started
###############

For the impatient
~~~~~~~~~~~~~~~~~

.. note::

   Note that the default config parameters are suited for a M<5 earthquake
   recorded within ~100 km. Adjust ``win_length``, ``noise_pre_time``, the
   frequency bands (``bp_freqmin_*``, ``bp_freqmax_*``, ``freq1_*``,
   ``freq2_*``) and the bounds on ``fc`` and ``t_star``, according to your
   problem.

Use case: miniSEED + StationXML + QuakeML
------------------------------------------

If you have seismic recordings in `miniSEED`_ format (e.g., ``traces.mseed``),
metadata in `StationXML`_ format (e.g., ``station.xml``) and event information
in `QuakeML`_ format (e.g., ``event.xml``), then:

1. Generate a config file via ``source_spec -S``;
2. Edit the config file variable ``station_metadata`` to point to
   ``station.xml`` file;
3. Run ``source_spec -t traces.mseed -q event.xml``.

Use case: SAC + PAZ + SourceSpec Event File
--------------------------------------------

If you have seismic recordings in `SAC`_ format (e.g., in a directory named
``sac_data``), metadata as `SAC polezero (PAZ)`_ (e.g., in a directory named
``paz``) and event information in any format, then:

1. Generate a config file via ``source_spec -S``;
2. Edit the config file variable ``station_metadata`` to point to the ``paz``
   directory;
3. Generate a sample :ref:`source_spec_event_file:SourceSpec Event File` using
   ``source_spec -y``; this will create a file named ``ssp_event.yaml``;
4. Edit the file ``ssp_event.yaml`` with your event information;
5. Run ``source_spec -t sac_data -H ssp_event.yaml``.


Command line arguments
~~~~~~~~~~~~~~~~~~~~~~

After successfully installed SourceSpec (see :ref:`installation:Installation`),
you can get help on the command line arguments used by each code by typing from
your terminal:

::

   source_spec -h

(or ``source_model -h``, or ``source_residuals -h``).

``source_spec`` and ``source_model`` require you to provide the path to
seismic traces via the ``--trace_path`` command line argument (see
:ref:`file_formats:File formats`).

Information on the seismic event can be stored in the trace header
(`SAC <https://ds.iris.edu/ds/support/faq/17/sac-file-format/>`__
format), or provided through a
`QuakeML <https://quake.ethz.ch/quakeml/>`__ file (``--qmlfile``) or,
alternatively (``--hypocenter``), through
a :ref:`source_spec_event_file:SourceSpec Event File`,
a `HYPO71 <https://pubs.er.usgs.gov/publication/ofr72224>`__ file, or
a `HYPOINVERSE-2000 <https://pubs.er.usgs.gov/publication/ofr02171>`__
file. See :ref:`file_formats:File Formats` for more
information on the supported file formats.

Configuration file
~~~~~~~~~~~~~~~~~~

``source_spec`` and ``source_model`` require a configuration file. The
default file name is ``source_spec.conf``, other file names can be
specified via the ``--configfile`` command line argument.

You can generate a sample configuration file through:

::

   source_spec -S

Take your time to go through the generated configuration file (named
``source_spec.conf``): the comments within the file will guide you on
how to set up the different parameters.

More details are in section :ref:`configuration_file:Configuration File`.


.. File format links:
.. _miniSEED: http://ds.iris.edu/ds/nodes/dmc/data/formats/miniseed/
.. _SAC: https://ds.iris.edu/ds/support/faq/17/sac-file-format/
.. _SAC file header: https://ds.iris.edu/files/sac-manual/manual/file_format.html
.. _QuakeML: https://quake.ethz.ch/quakeml/
.. _HYPO71: https://pubs.er.usgs.gov/publication/ofr72224
.. _HYPOINVERSE-2000: https://pubs.er.usgs.gov/publication/ofr02171
.. _StationXML: http://docs.fdsn.org/projects/stationxml/en/latest/
.. _Dataless SEED: https://ds.iris.edu/ds/nodes/dmc/data/formats/dataless-seed/
.. _SEED resp: https://ds.iris.edu/ds/nodes/dmc/data/formats/resp/
.. _SAC polezero (PAZ): https://www.jakewalter.net/sacresponse.html
.. _pickle: https://docs.python.org/3/library/pickle.html
.. _Cartopy: https://scitools.org.uk/cartopy/docs/latest
.. _SQLite: https://www.sqlite.org
.. _YAML: https://yaml.org