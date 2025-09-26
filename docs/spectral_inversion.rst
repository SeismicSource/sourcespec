.. _spectral_inversion_quality_and_uncertainty:

############################################
Spectral Inversion, Quality, and Uncertainty
############################################


Inversion algorithms
====================

The spectral model in SourceSpec is defined by three main parameters:
:math:`M_w` (moment magnitude),
:math:`f^{p|s}_c` (corner frequency for P or S waves), and
:math:`t^*` (attenuation).

See the :ref:`theoretical_background:Theoretical Background` section for
details on the spectral model.

Different algorithms can be used to invert for these parameters.
The algorithm is chosen with the ``inv_algorithm`` option in the
configuration file (see :ref:`configuration_file:Configuration File`).
The inversion is performed in moment magnitude :math:`M_w` units (logarithmic
amplitude).

Available inversion algorithms:

- **TNC**: `Truncated Newton algorithm <https://en.wikipedia.org/wiki/Truncated_Newton_method>`__ (with bounds).
- **LM**: `Levenberg–Marquardt algorithm <https://en.wikipedia.org/wiki/Levenberg–Marquardt_algorithm>`__.
  If bounds are provided, the `Trust Region Reflective algorithm <https://en.wikipedia.org/wiki/Trust_region>`__ is used instead.
- **BH**: `Basin-hopping algorithm <https://en.wikipedia.org/wiki/Basin-hopping>`__.
- **GS**: `Grid search <https://en.wikipedia.org/wiki/Hyperparameter_optimization#Grid_search>`__.
- **IS**: `Importance sampling <http://alomax.free.fr/nlloc/octtree/OctTree.html>`__ of the misfit grid, using a `k-d tree <https://en.wikipedia.org/wiki/K-d_tree>`__.


Inversion weighting
===================

The inversion can be weighted using different schemes, selected with the
``weighting`` option in the configuration file
(see :ref:`configuration_file:Configuration File`).

Weights are normalized between 0 and 1 and are shown in the output figures
(if plotting is enabled).

Available weighting schemes:

- **'noise'**: Applies weights based on spectral signal-to-noise ratio.
- **'frequency'**: Applies a constant weight for :math:`f \leq f_\text{weight}`;
  for :math:`f > f_\text{weight}`, a weight of 1 is applied.
  Controlled by the parameters ``f_weight`` and ``weight``.
- **'inv_frequency'**: Computes weights as :math:`1/(f-f_0+0.25)^{0.25}`
  for :math:`f \leq f_1`.  
  Weights are set to 0 for :math:`f < f_0` or :math:`f > f_1`.
  Here, ``f_0`` and ``f_1`` are the first and last frequencies where the
  signal-to-noise ratio exceeds 3 (or the spectrum bounds if no noise window
  is available).
- **'no_weight'**: No weighting (all weights = 1).


Spectrum quality
================

The quality of each spectrum is evaluated using the spectral
signal-to-noise ratio (SSNR):

.. math::

   SSNR(f) = \frac{S(f)}{N(f)}

where :math:`S(f)` is the signal amplitude spectrum and :math:`N(f)` is the
noise amplitude spectrum.

From the SSNR, two quality parameters are derived:

1. **SSNR_mean**: Mean SSNR within the full frequency band (or within the
   user-defined ``spectral_sn_freq_range``).
2. **SSNR_max**: Maximum SSNR within the same frequency band.

These values are reported in the output YAML file and in the SQLite database
(if enabled).  
See :ref:`output_file_formats:Output File Formats` for details.

A minimum threshold for **SSNR_mean** can be set with
``spectral_sn_min`` (see :ref:`configuration_file:Configuration File`).
Spectra below this threshold are excluded from the inversion.


Quality of spectral fit
=======================

Normalized root mean square error (RMSN)
----------------------------------------

The fit between observed and synthetic spectra is evaluated with the
normalized root mean square error (RMSN):

.. math::

   RMSN = \frac{RMS}{\sigma_w}

where :math:`RMS` is the root mean square error between observed and
synthetic spectra, and :math:`\sigma_w` is the weighted standard deviation of
the observed spectrum.

The RMS is given by:

.. math::
   \sqrt{\frac{1}{n} \sum_{i=1}^{n}
   \left(S_{obs}(f_i) - S_{syn}(f_i)\right)^2}

with :math:`S_{obs}(f_i)` and :math:`S_{syn}(f_i)` the observed and synthetic
amplitudes at frequency :math:`f_i`, and :math:`n` the number of frequency
samples.

The weighted standard deviation :math:`\sigma_w` is:

.. math::
   \sqrt{\frac{1}{\sum_{i=1}^{n} w_i}
   \sum_{i=1}^{n} w_i
   \left(S_{obs}(f_i) - \overline{S_{obs}}\right)^2}

where :math:`w_i` is the weight for each frequency sample.

RMSN values are reported in the output YAML file and SQLite database
(if enabled).


Quality of fit
--------------

From the RMSN, a fit quality parameter :math:`Q_{fit}` is defined:

.. math::

    Q_{fit} = 100 \times e^{-RMSN}

:math:`Q_{fit}` ranges from 0 to 100, with higher values indicating better fits.
It is reported in the output YAML file and SQLite database (if enabled).

The minimum acceptable fit quality can be set with
``pi_quality_of_fit_min`` (see :ref:`configuration_file:Configuration File`).
Fits below this threshold are excluded from the final results.


Spectral parameters: single spectra and event summaries
=======================================================

For each inverted spectrum, SourceSpec computes the following parameters
(see :ref:`theoretical_background:Theoretical Background` for details):

- **Moment magnitude** (:math:`M_w`): obtained directly from inversion.
- **Corner frequency** (:math:`f_c`): obtained directly from inversion.
- **Attenuation parameter** (:math:`t^*`): obtained directly from inversion.
- **Radiated energy** (:math:`E_r`): measured directly from the spectrum.
- **Seismic moment** (:math:`M_0`): derived from :math:`M_w`.
- **Source radius** (:math:`a`): derived from :math:`M_0` and :math:`f_c`.
- **Static stress drop** (:math:`\Delta\sigma`): derived from :math:`M_0` and :math:`a`.
- **Apparent stress** (:math:`\sigma_a`): derived from :math:`E_r` and :math:`M_0`.
- **Quality factor** (:math:`Q`): derived from :math:`t^*`.

For each spectrum, both the parameter estimates and their uncertainties (from
the chosen inversion algorithm) are stored in the output YAML file and in the
SQLite database (if enabled).
See :ref:`output_file_formats:Output File Formats` for details.

In addition, **event-level summary values** and associated uncertainties are
computed across spectra using the mean, weighted mean, or user-specified
percentiles. Outliers are rejected using the
`interquartile range <https://en.wikipedia.org/wiki/Interquartile_range>`__
(IQR) method (see the parameter ``nIQR`` in the configuration file).

A boxplot of the distribution of spectral parameters across stations is
produced if plotting is enabled.


General quality of the inversion
================================

The overall reliability of the inversion is characterized by several
summary parameters:

- **Number of spectra inverted**: The number of spectra that passed quality
  checks and contributed to the inversion.
- **Primary azimuthal gap**: the largest gap in station azimuthal coverage;
  measures the overall quality of the distribution.
- **Secondary azimuthal gap**: the largest gap that would remain if any one
  station were removed; measures the robustness of the coverage.
- **Mean RMSN**: The average normalized root mean square error (RMSN) across
  all inverted spectra (see the RMSN definition above).
- **Mean fit quality**: The average spectral fit quality (:math:`Q_{fit}`)
  across all inverted spectra (see the definition of :math:`Q_{fit}` above).
- **Spectral dispersion (RMSN)**: The normalized RMSN quantifying how much
  the observed spectra deviate from the event summary spectrum (see details
  below).
- **Spectral dispersion score**: A quality score derived from the spectral
  dispersion RMSN (see details below).

All these parameters are written to the output YAML file and the SQLite
database (if enabled).

Spectral dispersion (RMSN)
--------------------------

The spectral dispersion RMSN quantifies how much individual spectra deviate
from the event summary spectrum. It is defined as:

.. math::
    RMSN_{disp} = \frac{RMS_{disp}}{IPR_{80}}

Here, :math:`RMS_{disp}` is the root mean square error between all observed
spectra and the event summary spectrum:

.. math::
    RMS_{disp} = \sqrt{\frac{1}{mn} \sum_{j=1}^{n} \sum_{i=1}^{m}
    \left(S_{obs,j}(f_i) - S_{event}(f_i)\right)^2}

where:

- :math:`S_{obs,j}(f_i)` = observed amplitude of spectrum :math:`j` at frequency :math:`f_i`,
- :math:`S_{event}(f_i)` = amplitude of the event summary spectrum at frequency :math:`f_i`,
- :math:`m` = number of frequency samples per spectrum,
- :math:`n` = number of spectra.

The denominator, :math:`IPR_{80}`, is the 80% inter-percentile range of all
data points in the observed spectra, i.e., the range between the 10th and 90th
percentiles.

Spectral dispersion score
-------------------------

This score is derived from the spectral dispersion RMSN
(:math:`RMSN_{disp}`) and ranges from 0 to 100, with higher values indicating
better overall fit among spectra.

It is defined as:

.. math::
    Score_{disp} = 100 \times e^{-RMSN_{disp}}
