.. _output_file_formats:

###################
Output File Formats
###################

The SourceSpec main code, ``source_spec`` will produce the following
output files (``EVID`` is replaced by the actual event ID):

-  ``EVID.ssp.yaml``: `YAML`_ file containing the estimated spectral parameters
   (summary values and per station values)
-  ``EVID.ssp.out`` (*deprecated*): text file containing the estimated spectral
   parameters (summary values and per station values)
-  ``EVID.ssp.log``: log file in text format (including the command line
   arguments, for `reproducibility`_)
-  ``EVID.ssp.conf``: the input config file (for `reproducibility`_)
-  ``EVID.residuals.hdf5``: station residuals in
   :ref:`spectral_file_formats:HDF5 File Format`
-  ``EVID.spectra.hdf5``: (optional) spectra in
   :ref:`spectral_file_formats:HDF5 File Format`
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
.. _Cartopy: https://scitools.org.uk/cartopy/docs/latest
.. _SQLite: https://www.sqlite.org
.. _YAML: https://yaml.org

.. Method links:
.. _reproducibility: https://en.wikipedia.org/wiki/Reproducibility
.. _box plots: https://en.wikipedia.org/wiki/Box_plot