.. _theoretical_background:

######################
Theoretical Background
######################

Overview
========

``source_spec`` inverts the P- or S-wave displacement amplitude spectra from
station recordings of a single event.

For each station, the code computes P- or S-wave displacement amplitude spectra
for each component (e.g., Z, N, E), within a predefined time window.

.. math::

   S(f) = \left| \int_{t0}^{t1} s(t) e^{-i 2 \pi f t} dt \right|

Note that the Fourier amplitude spectrum of ground displacement :math:`S(f)`
has the dimensions of displacement (:math:`s(t)`) multiplied by
time (:math:`dt`).

The same thing is done for a noise time window: noise spectrum is used to
compute spectral signal-to-noise ratio (and possibly reject low SNR spectra)
and, optionally, to weight the spectral inversion.

.. figure:: imgs/example_trace.svg
  :alt: Example trace plot
  :width: 600

  Example three-component trace plot (in velocity), showing noise and S-wave
  windows.

The component spectra are combined through the root sum of squares:

.. math::

    S(f) = \sqrt{S^2_z(f) + S^2_n(f) + S^2_e(f)}

where :math:`f` is the frequency and :math:`S_x(f)` is the P- or S-wave
spectrum for component :math:`x`.

.. figure:: imgs/example_spectrum.svg
  :alt: Example spectrum plot
  :width: 600

  Example displacement spectrum for noise and S-wave, including inversion
  results.


Spectral model
==============

The Fourier amplitude spectrum of the P- or S-wave displacement in far field
can be modelled as the product of a source term :cite:p:`Brune1970` and a
propagation term (geometric and anelastic attenuation of body waves):

.. math::

   S(f) =
          \frac{1}{\mathcal{G}(r)}
          \times
          \frac{2 R_{\Theta\Phi}}
               {4 \pi \rho_h^{1/2} \rho_r^{1/2} c_h^{5/2} c_r^{1/2}}
          \times
          M_O
          \times
          \frac{1}{1+\left(\frac{f}{f_c}\right)^2}
          \times
          e^{- \pi f t^*}

where:

- :math:`\mathcal{G}(r)` is the geometrical spreading coefficient (see below)
  and :math:`r` is the hypocentral distance;
- the coefficient :math:`2` is the free surface amplification factor;
- :math:`R_{\Theta\Phi}` is the radiation pattern coefficient for P- or S-waves
  (average or computed from focal mechanism, if available);
- :math:`\rho_h` and :math:`\rho_r` are the medium densities at the hypocenter
  and at the receiver, respectively;
- :math:`c_h` and :math:`c_r` are the P- or S-wave velocities at the hypocenter
  and at the receiver, respectively;
- :math:`M_O` is the seismic moment;
- :math:`f` is the frequency;
- :math:`f_c` is the corner frequency;
- :math:`t^*` is an attenuation parameter which includes anelastic path
  attenuation (quality factor) and station-specific effects.

Geometrical spreading
---------------------
The geometrical spreading coefficient :math:`\mathcal{G}(r)` can be defined in
one of the following ways (see the ``geom_spred_model`` option in
:ref:`configuration_file:Configuration File`):

- :math:`\mathcal{G}(r) = r^n`: :math:`n` can be any positive number.
  :math:`n=1` (default value) is the theoretical value for a body wave in a
  homogeneous full-space;
  :math:`n=0.5` is the theoretical value for a surface wave in a homogeneous
  half-space.

- Following :cite:t:`Boatwright2002` (eq. 8), to account for the mixture of
  body waves, Lg waves and surface waves at regional distances
  (:math:`r < 200 km`), a two-part geometrical spreading coefficient:

  - body wave spreading (:math:`\mathcal{G}(r) = r`) for hypocentral distances
    below a cutoff distance :math:`r_0`;
  - frequency dependent spreading for hypocentral distances above the
    cutoff distance :math:`r_0`.

More precisely, the expression derived from :cite:t:`Boatwright2002` is:

.. math::

  \mathcal{G}(r) =
  \begin{cases}
    r  &  r \le r_0\\
    r_0 (r/r_0)^{\gamma (f)}  &  r > r_0
  \end{cases}

with

.. math::

  \gamma (f) =
  \begin{cases}
    0.5  &  f \le 0.20 Hz\\
    0.5 + 2 \log_{10} (5f)  &  0.20 < f < 0.25 Hz\\
    0.7  &  f \ge 0.25 Hz\\
  \end{cases}

Note that here we use the square root of eq. 8 in :cite:t:`Boatwright2002`,
since we correct the spectral amplitude and not the energy.


Building spectra
================

In ``source_spec``, the observed spectrum of component :math:`x` (vertical or
horizontal), :math:`S_x(f)` is converted into moment magnitude units
:math:`M_w`.

The first step is to multiply the spectrum for the geometrical spreading
coefficient and convert it to seismic moment units:

.. math::

   M_x(f) \equiv
   \mathcal{G}(r) \times
   \frac{4 \pi \rho_h^{1/2} \rho_r^{1/2} c_h^{5/2} c_r^{1/2}}
        {2 R_{\Theta\Phi}}
   \times S_x(f) =
          M_O \times
          \frac{1}{1+\left(\frac{f}{f_c}\right)^2}
          \times
          e^{- \pi f t^*}


Then the spectrum is converted in units of magnitude
(the :math:`Y_x (f)` vector used in the inversion):

.. math::

   Y_x(f) \equiv
          \frac{2}{3} \times
          \left( \log_{10} M_x(f) - 9.1 \right)

The data vector is compared to the theoretical model:

.. math::

   Y_x(f) =
          \frac{2}{3}
          \left[ \log_{10} \left(
                    M_O \times
                    \frac{1}{1+\left(\frac{f}{f_c}\right)^2}
                    \times
                    e^{- \pi f t^*}
                    \right) - 9.1 \right] =

          =
          \frac{2}{3} (\log_{10} M_0 - 9.1) +
          \frac{2}{3} \left[ \log_{10} \left(
                    \frac{1}{1+\left(\frac{f}{f_c}\right)^2} \right) +
                    \log_{10} \left( e^{- \pi f t^*} \right)
                    \right]


Finally coming to the following model used for the inversion:

.. math::

   Y_x(f) =
          M_w +
          \frac{2}{3} \left[ - \log_{10} \left(
                    1+\left(\frac{f}{f_c}\right)^2 \right) -
                    \pi \, f t^* \log_{10} e
                    \right]

where :math:`M_w \equiv \frac{2}{3} (\log_{10} M_0 - 9.1)`.


Inverted parameters
===================

The parameters determined from the spectral inversion are :math:`M_w`,
:math:`f_c` and :math:`t^*`.

The inversion is performed in moment magnitude :math:`M_w` units (logarithmic
amplitude). Different inversion algorithms can be used:

-  TNC: `truncated Newton
   algorithm <https://en.wikipedia.org/wiki/Truncated_Newton_method>`__
   (with bounds)
-  LM: `Levenberg-Marquardt
   algorithm <https://en.wikipedia.org/wiki/Levenberg–Marquardt_algorithm>`__
   (warning: `Trust Region Reflective
   algorithm <https://en.wikipedia.org/wiki/Trust_region>`__ will be
   used instead if bounds are provided)
-  BH: `basin-hopping
   algorithm <https://en.wikipedia.org/wiki/Basin-hopping>`__
-  GS: `grid
   search <https://en.wikipedia.org/wiki/Hyperparameter_optimization#Grid_search>`__
-  IS: `importance
   sampling <http://alomax.free.fr/nlloc/octtree/OctTree.html>`__ of
   misfit grid, using `k-d
   tree <https://en.wikipedia.org/wiki/K-d_tree>`__


Other computed parameters
=========================

Starting from the inverted parameters :math:`M_0` ( :math:`M_w` ),
:math:`fc`, :math:`t^*` and following the equations in :cite:t:`Madariaga2011`
and :cite:t:`Lancieri2012`, other quantities are computed for each station:

-  the Brune static stress drop :math:`\Delta \sigma`
-  the source radius :math:`a`
-  the radiated energy :math:`E_r`
-  the apparent stress :math:`\sigma_a`
-  the quality factor :math:`Q_0` of P- or S-waves

As a bonus, local magnitude :math:`M_l` can be computed as well.

Event summaries (mean, weighted mean, percentiles) are computed from single
station estimates. For mean and weighted mean estimation, outliers are rejected
based on the `interquartile
range <https://en.wikipedia.org/wiki/Interquartile_range>`__ rule.


Source radius and Brune static stress drop
------------------------------------------
The Brune static stress drop :math:`\Delta \sigma` is computed under the
assumption of a circular rupture of radius :math:`a`. The model of
:cite:t:`Brune1970` provides an expression for the source radius (equation 31
in :cite:t:`Madariaga2011`):

.. math::

   a = 0.3724 \frac{\beta_h}{f_c}

where :math:`\beta_h` is the S-wave velocity at the hypocenter (in :math:`m/s`)
and :math:`f_c` is the corner frequency (in :math:`Hz`) estimated from the
spectral inversion.

The Brune static stress drop is then computed using the circular crack model,
as discussed in :cite:t:`Madariaga2011` (equation 27):

.. math::

   \Delta \sigma =
   \frac{7}{16}
   \frac{M_0}{a^3}

where :math:`M_0` is the seismic moment (in :math:`N \cdot m`) and
:math:`a` is the source radius (in :math:`m`).


Radiated energy
---------------
The radiated energy :math:`E_r` is computed from the integral of the squared
velocity amplitude spectrum: :math:`\dot{S}^2(f) = [ 2 \pi f S(f) ]^2`.

Following :cite:t:`Boatwright2002` (equation 1) and :cite:t:`Lancieri2012`
(equation 3), the radiated energy is computed as:

.. math::

   E_r = 8 \pi \mathcal{G}^2(r) C^2 \rho_r c_r
            \int_{0}^{f_{max}} e^{2 \pi f t^*} \dot{S}^2(f) df

where :math:`\mathcal{G}^2(r)` is the squared geometrical spreading coefficient
(see above), :math:`C` is a constant discussed below, :math:`\rho_r` and
:math:`c_r` are, respectively, the density and P- or S-wave velocity [#f1]_
at the receiver (their product is the seismic impedance), :math:`f_{max}` is
the maximum frequency used to compute the energy (see
:ref:`configuration_file:Configuration File` for details on the ``max_freq_Er``
parameter), and the exponential term in the integrand is the squared correction
for anelastic attenuation.

The constant :math:`C` is defined in :cite:t:`Boatwright2002` (equation 2) as:

.. math::

   C = \frac{\left<R_{\Theta\Phi}\right>}{R_{\Theta\Phi} F}

where :math:`\left<R_{\Theta\Phi}\right>` is the root mean square radiation
pattern computed on the focal sphere, :math:`R_{\Theta\Phi}` is the
radiation pattern coefficient for the given station, and :math:`F` is the
free surface amplification factor.
Here we assume :math:`F = 2` and :math:`\left<R_{\Theta\Phi}\right> = 1`
(hence, :math:`C = 1/2`).
The latter assumption means that we rely on the averaging between measurements
of radiated energy at different stations, instead of precise measurements at a
single station.

Noise correction
++++++++++++++++
To account for low frequency noise, below the corner frequency, under the
hypothesis that energy is additive and that noise is stationary, we compute
a corrected energy as:

.. math::

   \tilde{E}_r = E_r- E_{r,noise}

where :math:`E_r` is the observed radiated energy and :math:`E_{r,noise}` is
the radiated energy computed from the noise spectrum.

Finite bandwidth correction
+++++++++++++++++++++++++++
The final step is to correct the radiated energy for the finite bandwidth
of the observed spectrum. Following :cite:t:`Lancieri2012` (equation 4), and
:cite:t:`DiBona1988`, the noise-corrected radiated energy is divided by
the following factor:

.. math::

  R = \frac{2}{\pi}
    \left[
      \frac{-f_{max}/f_c}{1+(f_{max}/f_c)^2} + \arctan(f_{max}/f_c)
    \right]

where :math:`f_c` is the corner frequency and :math:`f_{max}` is the maximum
frequency used to compute the energy.


Apparent stress
---------------

The apparent stress :math:`\sigma_a` is computed as (:cite:t:`Madariaga2011`,
eq. 18):

.. math::

   \sigma_a = \mu_h \frac{E_r}{M_0}

where :math:`\mu_h` is the shear modulus (or rigidity, in :math:`Pa`) near the
hypocenter, :math:`E_r` is the radiated energy (in :math:`N \cdot m`), and
:math:`M_0` is the seismic moment (in :math:`N \cdot m`).

The value of :math:`\mu_h` is computed from the shear wave velocity
(:math:`\beta_h`) and the density (:math:`\rho_h`) at the hypocenter,
using the following expression:

.. math::

   \mu_h = \rho_h \beta_h^2


Quality factor
--------------
The retrieved attenuation parameter :math:`t^*` is converted to the P- or
S-wave quality factor :math:`Q_0^{[p|s]}` using the following expression:

.. math::

   Q_0^{[p|s]} = \frac{tt_{[p|s]}(r)}{t^*}

where :math:`tt_{[p|s]}(r)` is the P- or S-wave travel time from source to
station and :math:`r` is the hypocentral distance.


Station Residuals
-----------------
Station-specific effects can be determined by running ``source_spec`` on several
events and computing the average of station residuals between observed and
inverted spectra. These averages are obtained through the command
``source_residuals``; the resulting residuals file can be used for a second run
of ``source_spec`` (see the ``residuals_filepath`` option in
:ref:`configuration_file:Configuration File`).


.. rubric:: Footnotes

.. [#f1] SourceSpec can compute radiated energy from either the P- or S-wave
   displacement spectra, depending on the value chosen for the configuration
   parameter ``wave_type`` (see :ref:`configuration_file:Configuration File`).
   However, when using P waves, the code will warn that radiated energy
   computed from P waves might be underestimated.