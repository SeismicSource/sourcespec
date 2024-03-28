.. SourceSpec documentation master file, created by
   sphinx-quickstart on Fri Oct 25 17:47:32 2013.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

SourceSpec documentation
========================

Earthquake source parameters from P- or S-wave displacement spectra

:Copyright: 2011-2023 Claudio Satriano satriano@ipgp.fr
:Release: |release|
:Date:    |today|

SourceSpec is a collection of command line tools to compute earthquake source
parameters (seismic moment, corner frequency, radiated energy, source size,
static stress drop, apparent stress) from the inversion of P-wave and S-wave
displacement spectra recorded at one or more seismic stations.
SourceSpec also computes attenuation parameters (t-star, quality factor) and,
as a bonus, local magnitude.

See :cite:t:`Madariaga2011` for a primer on earthquake source parameters and
scaling laws.

Go to section :ref:`theoretical_background:Theoretical Background`
to get more information on how the code works.

SourceSpec is written in Python and requires a working Python
environment to run (see :ref:`installation:Installation`).
However, since SourceSpec is based on command line, you don’t have to
know how to code in Python to use it.

The SourceSpec package is made of several command line tools:

-  ``source_spec``: Compute earthquake source parameters from the
   inversion of P- or S-wave spectra.
-  ``source_model``: Direct modelling of P- or S-wave spectra, based on
   user-defined earthquake source parameters.
-  ``source_residuals``: Compute station residuals from ``source_spec``
   output.
- ``clipping_detection``: Test the clipping detection algorithm.
- ``plot_sourcepars``: 1D or 2D plot of source parameters from a sqlite
  parameter file.


Contents:

.. toctree::
   :maxdepth: 2

   theoretical_background
   signal_processing
   clipping_detection
   getting_started
   configuration_file
   file_formats
   source_spec_event_file
   spectral_file_formats
   installation
   sample_runs
   getting_help
   contributing
   citing
   api
   changelog
   bibliography
   citing_literature


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
