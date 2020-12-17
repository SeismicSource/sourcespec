.. SourceSpec documentation master file, created by
   sphinx-quickstart on Fri Oct 25 17:47:32 2013.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

SourceSpec documentation
========================

SourceSpec is a collection of command line programs (written in Python) to
determine earthquake source parameters (seismic moment :math:`M_0`, corner
frequency :math:`f_c`) and the inelastic attenuation term (:math:`t^*`), from
the modeling of waveform spectra.

Other parameters (source radius :math:`r_0`, stress drop :math:`\Delta \sigma`)
are computed from the inverted ones. The quality factor :math:`Q` is determined
from :math:`t^*`.

As a bonus, local magnitude :math:`M_l` is computed as well.

SourceSpec is composed of the following programs:

* ``source_spec``: inverts the S-wave displacement spectra from station
  recordings of a single event.
* ``ssp_residuals``: computes station residuals from ``source_spec`` output.
* ``source_model``: direct spectral modelling.


Contents:

.. toctree::
   :maxdepth: 2

   source_spec
   configuration_file
   api


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
