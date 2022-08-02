.. _getting_started:

###############
Getting started
###############

For the impatient
~~~~~~~~~~~~~~~~~

If you have seismic recordings in
`miniSEED <http://ds.iris.edu/ds/nodes/dmc/data/formats/miniseed/>`__
format (e.g., ``traces.mseed``), metadata in
`StationXML <http://docs.fdsn.org/projects/stationxml/en/latest/>`__
format (e.g., ``station.xml``) and event information in
`QuakeML <https://quake.ethz.ch/quakeml/>`__ format (e.g.,
``event.xml``), then:

1. Generate a config file via ``source_spec -S``;
2. Edit the config file variable ``station_metadata`` to point to
   ``station.xml`` file;
3. Run ``source_spec -t traces.mseed -q event.xml``.

Command line arguments
~~~~~~~~~~~~~~~~~~~~~~

After successfully installed SourceSpec (see
`Installation <#installation>`__ below), you can get help on the command
line arguments used by each code by typing from your terminal:

::

   source_spec -h

(or ``source_model -h``, or ``source_residuals -h``).

``source_spec`` and ``source_model`` require you to provide the path to
seismic traces via the ``--trace_path`` command line argument (see
`Supported file formats <#supported-file-formats>`__ below).

Information on the seismic event can be stored in the trace header
(`SAC <https://ds.iris.edu/ds/support/faq/17/sac-file-format/>`__
format), or provided through a
`QuakeML <https://quake.ethz.ch/quakeml/>`__ file (``--qmlfile``) or a
`HYPO71 <https://pubs.er.usgs.gov/publication/ofr72224>`__ or
`HYPOINVERSE-2000 <https://pubs.er.usgs.gov/publication/ofr02171>`__
file (``--hypocenter``). See :ref:`supported_file_formats:Supported File
Formats` for more information on the supported file formats.

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