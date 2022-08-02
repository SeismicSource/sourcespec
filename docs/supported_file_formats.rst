.. _supported_file_formats:

######################
Supported File Formats
######################

Trace formats
~~~~~~~~~~~~~

SourceSpec can read all the `trace formats supported by
ObsPy <https://docs.obspy.org/packages/autogen/obspy.core.stream.read.html>`__.

Two very common choices are:

-  `miniSEED <http://ds.iris.edu/ds/nodes/dmc/data/formats/miniseed/>`__
-  `SAC <https://ds.iris.edu/ds/support/faq/17/sac-file-format/>`__

The SAC format can carry additional information in its header, like
event location and origin time, phase picks, instrument sensitivity.

Event formats
~~~~~~~~~~~~~

SourceSpec can read event information (event ID, location, origin time)
in the following formats:

-  `QuakeML <https://quake.ethz.ch/quakeml/>`__:
   SourceSpec will also read phase picks and focal mechanism, if available
-  `HYPO71 <https://pubs.er.usgs.gov/publication/ofr72224>`__
-  `HYPOINVERSE-2000 <https://pubs.er.usgs.gov/publication/ofr02171>`__:
   SourceSpec will also read phase picks, if available

Event information can also be stored in the SAC file headers (header
fields: ``EVLA``, ``EVLO``, ``EVDP``, ``O``, ``KEVNM``).

Phase pick formats
~~~~~~~~~~~~~~~~~~

Phase picks for P and S waves can be read from one of the following
formats:

-  `QuakeML <https://quake.ethz.ch/quakeml/>`__
-  `HYPO71 <https://pubs.er.usgs.gov/publication/ofr72224>`__
-  `HYPOINVERSE-2000 <https://pubs.er.usgs.gov/publication/ofr02171>`__

Phase picks can also be stored in the SAC file headers (header fields:
``A`` and ``T0``).

Station metadata formats
~~~~~~~~~~~~~~~~~~~~~~~~

Station metadata (coordinates, instrumental response) can be provided in
one of the following formats:

-  `StationXML <http://docs.fdsn.org/projects/stationxml/en/latest/>`__
-  `Dataless
   SEED <https://ds.iris.edu/ds/nodes/dmc/data/formats/dataless-seed/>`__
-  `SEED RESP <https://ds.iris.edu/ds/nodes/dmc/data/formats/resp/>`__
-  `SAC polezero (PAZ) <https://www.jakewalter.net/sacresponse.html>`__

Note that SEED RESP and PAZ formats do not contain station coordinates,
which should therefore be in the trace header (traces in SAC format).

The station metadata file name or file directory is provided in the
configuration file through the parameter ``station_metadata``.

Alternatively, instrument sensitivity can be provided in the SAC header
or as a constant in the configuration file. In both cases, use the
configuration parameter ``sensitivity``.

Output files
~~~~~~~~~~~~

The SourceSpec main code, ``source_spec`` will produce the following
output files (``EVID`` is replaced by the actual event ID):

-  ``EVID.ssp.out``: text file containing the estimated earthquake
   source parameters (per station and average)
-  ``EVID.ssp.log``: log file in text format (including the command line
   arguments, for
   `reproducibility <https://en.wikipedia.org/wiki/Reproducibility>`__)
-  ``EVID.ssp.conf``: the input config file (for
   `reproducibility <https://en.wikipedia.org/wiki/Reproducibility>`__)
-  ``EVID-residuals.pickle``: station residuals in `Python pickle
   format <https://docs.python.org/3/library/pickle.html>`__
-  ``EVID.xml``: updated StationXML file with the results of the
   SourceSpec inversion (only if an input StationXML file is provided)

The following plots will be created, in png or pdf format:

-  ``EVID.traces.png[.pdf]``: trace plots
-  ``EVID.ssp.png[.pdf]``: spectral plots
-  ``EVID.sspweight.png[.pdf]``: spectral weight plots
-  ``EVID.boxplot.png[.pdf]``: `box
   plots <https://en.wikipedia.org/wiki/Box_plot>`__ for the earthquake
   source parameters retrieved at each station
-  Misfit plots, when using “grid search” or “importance sampling” for
   the spectral inversion

As an option, station maps can be created (requires
`Cartopy <https://scitools.org.uk/cartopy/docs/latest>`__):

-  ``EVID.map_mag.png[.pdf]``: station map with symbols colored by
   estimated moment magnitude
-  ``EVID.map_fc.png[.pdf]``: station map with symbols colored by
   estimated corner frequency

As an option, the retrieved source parameters (per station and average)
can be appended to a `SQLite <https://www.sqlite.org>`__ database, whose
path is defined in the configuration file.

Finally, always as an option, ``source_spec`` can generate a report in
HTML format.