.. _file_formats:

############
File Formats
############

Trace formats
~~~~~~~~~~~~~

SourceSpec can read all the `trace formats supported by
ObsPy <https://docs.obspy.org/packages/autogen/obspy.core.stream.read.html>`__.

Two very common choices are:

-  `miniSEED`_
-  `SAC`_

The SAC format can carry additional information in its header, like
event location and origin time, phase picks, instrument sensitivity.

Event formats
~~~~~~~~~~~~~

SourceSpec can read event information (event ID, location, origin time)
in the following formats:

-  `QuakeML`_: SourceSpec will also read phase picks and focal mechanism,
   if available
-  `HYPO71`_
-  `HYPOINVERSE-2000`_: SourceSpec will also read phase picks, if available

Event information can also be stored in the `SAC file header`_ (header
fields: ``EVLA``, ``EVLO``, ``EVDP``, ``O``, ``KEVNM``).

Phase pick formats
~~~~~~~~~~~~~~~~~~

Phase picks for P and S waves can be read from one of the following
formats:

-  `QuakeML`_
-  `HYPO71`_
-  `HYPOINVERSE-2000`_

Phase picks can also be stored in the `SAC file header`_, using the header
fields ``A`` and ``T0`` through ``T9``. A pick label can be specified (header
fields ``KA`` and ``KT0`` through ``KT9``) to identify the pick; the pick label
can be a standard 4-characters SAC label (e.g., ``"IPU0"``, ``" S 1"``) or a
label starting with ``"P"`` or ``"S"`` (lowercase or uppercase, e.g., ``"P"``,
``"pP"``, ``"Pg"``, ``"S"``, ``"Sn"``).
Picks with labels that cannot be parsed by SourceSpec will be ignored.
If no label is specified, then SourceSpec will assume that ``A`` is the P-pick
and ``T0`` is the S-pick.

Station metadata formats
~~~~~~~~~~~~~~~~~~~~~~~~

Station metadata (coordinates, instrumental response) can be provided in
one of the following formats:

-  `StationXML`_
-  `Dataless SEED`_
-  `SEED RESP`_
-  `SAC polezero (PAZ)`_

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

-  ``EVID.ssp.yaml``: `YAML`_ file containing the estimated spectral parameters
   (summary values and per station values)
-  ``EVID.ssp.out`` (*deprecated*): text file containing the estimated spectral
   parameters (summary values and per station values)
-  ``EVID.ssp.log``: log file in text format (including the command line
   arguments, for `reproducibility`_)
-  ``EVID.ssp.conf``: the input config file (for `reproducibility`_)
-  ``EVID-residuals.pickle``: station residuals in Python `pickle`_ format
-  ``EVID.ssp.h``: hypocenter file in `HYPO71`_ format with the estimated
   moment magnitude (only if an input HYPO71 file is provided)
-  ``EVID.xml``: updated `QuakeML`_ file with the results of the SourceSpec
   inversion (only if an input QuakeML file is provided)

The following plots will be created, in png, pdf or svg format:

-  ``EVID.traces.png[.pdf,.svg]``: trace plots
-  ``EVID.ssp.png[.pdf,.svg]``: spectral plots
-  ``EVID.sspweight.png[.pdf,.svg]``: spectral weight plots
-  ``EVID.boxplot.png[.pdf,.svg]``: `box plots`_ for the earthquake source
   parameters retrieved at each station
-  Misfit plots, when using “grid search” or “importance sampling” for
   the spectral inversion

As an option, station maps can be created (requires `Cartopy`_):

-  ``EVID.map_mag.png[.pdf,.svg]``: station map with symbols colored by
   estimated moment magnitude
-  ``EVID.map_fc.png[.pdf,.svg]``: station map with symbols colored by
   estimated corner frequency

As an option, the retrieved source parameters (per station and average)
can be appended to a `SQLite`_ database, whose path is defined in the
configuration file.

Finally, always as an option, ``source_spec`` can generate a report in
HTML format.

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

.. Method links:
.. _reproducibility: https://en.wikipedia.org/wiki/Reproducibility
.. _box plots: https://en.wikipedia.org/wiki/Box_plot