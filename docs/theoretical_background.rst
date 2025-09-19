.. _theoretical_background:

######################
Theoretical Background
######################

Overview
========

``source_spec`` inverts the P- or S-wave displacement amplitude spectra from
station recordings of a single event.

For each station, the code computes P- or S-wave displacement amplitude spectra
for each component :math:`x` (e.g., Z, N, E), within a predefined time window.

.. math::

   S_x^{p|s}(f) =
      \left| \int_{t^{p|s}_0}^{t^{p|s}_1} s_x(t) e^{-i 2 \pi f t} dt \right|

where the exponent :math:`p|s` means that we are considering either P- or
S-waves, :math:`t^{p|s}_0` and :math:`t^{p|s}_1` are the start and end times of
the P- or S-wave time window, :math:`s_x(t)` is the displacement time series
for component :math:`x`, and :math:`f` is the frequency.

Note that the Fourier amplitude spectrum of ground displacement has the
dimensions of displacement (:math:`s_x(t)`) multiplied by time (:math:`dt`).

The same thing is done for a noise time window: noise spectrum is used to
compute spectral signal-to-noise ratio (and possibly reject low SNR spectra)
and, optionally, to weight the spectral inversion.

.. figure:: imgs/example_trace.svg
  :alt: Example trace plot
  :width: 600

  Example three-component trace plot (in velocity), showing noise and S-wave
  windows.

The component spectra are combined through the root sum of squares
(e.g., Z, N, E):

.. math::

    S^{p|s}(f) =
      \sqrt{
         \left( S^{p|s}_z(f) \right)^2 +
         \left( S^{p|s}_n(f) \right)^2 +
         \left( S^{p|s}_e(f) \right)^2
      }

(This is actually done later in the code, after converting the spectra to
magnitude units, see below.)

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

   S^{p|s}(f) =
          \frac{1}{\mathcal{G}(r)}
          \times
          \frac{F R_{\Theta\Phi}}
               {4 \pi \rho_h^{1/2} \rho_r^{1/2} c_h^{5/2} c_r^{1/2}}
          \times
          M_O
          \times
          \frac{1}{1+\left(\frac{f}{f^{p|s}_c}\right)^2}
          \times
          e^{- \pi f t^*}

where:

- :math:`\mathcal{G}(r)` is the geometrical spreading coefficient (see below)
  and :math:`r` is the hypocentral distance;
- :math:`F` is the free surface amplification factor (generally assumed to be
  :math:`2`);
- :math:`R_{\Theta\Phi}` is the radiation pattern coefficient for P- or S-waves
  (average or computed from focal mechanism, if available);
- :math:`\rho_h` and :math:`\rho_r` are the medium densities at the hypocenter
  and at the receiver, respectively;
- :math:`c_h` and :math:`c_r` are the P- or S-wave velocities at the hypocenter
  and at the receiver, respectively;
- :math:`M_O` is the seismic moment;
- :math:`f` is the frequency;
- :math:`f^{p|s}_c` is the corner frequency for P- or S-waves;
- :math:`t^*` is an attenuation parameter which includes anelastic path
  attenuation (quality factor) and station-specific effects.

Geometrical spreading
---------------------
The geometrical spreading coefficient :math:`\mathcal{G}(r)` can be defined,
for local and regional distances, in one of the following ways (see the
``geom_spred_model`` option in :ref:`configuration_file:Configuration File`):

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

For teleseismic distances (see the option
``geom_spread_min_teleseismic_distance``
in :ref:`configuration_file:Configuration File`), the geometrical spreading
coefficient is defined as in :cite:t:`Okal1992` (eq. 4):

.. math::

   \mathcal{G}(\Delta) = \frac{a}{g(\Delta)}

where :math:`\Delta` is the great circle distance between the source and the
receiver, :math:`a` is the Earth radius and :math:`g(\Delta)` is defined as:

.. math::

   g(\Delta) = \left(
      \frac{\rho_h c_h}{\rho_r c_r}
      \frac{\sin i_h}{\sin \Delta}
      \frac{1}{\cos i_r}
      \left| \frac{d i_h}{d \Delta} \right|
   \right)^{1/2}

where :math:`\rho_h` and :math:`\rho_r` are the medium densities at the
hypocenter and at the receiver, respectively, :math:`c_h` and :math:`c_r` are
the P- or S-wave velocities at the hypocenter and at the receiver,
respectively, :math:`i_h` and :math:`i_r` are the takeoff angle (hypocenter) and
the incidence angle (receiver), respectively,
and :math:`\frac{d i_h}{d \Delta}` is the variation of the takeoff angle within
a ray tube of width :math:`\Delta` (see :cite:t:`Okal1992` for details).

Building spectra
================

In ``source_spec``, the observed spectrum of component :math:`x` (vertical or
horizontal), :math:`S^{p|s}_x(f)` is converted into moment magnitude units
:math:`M_w`.

The first step is to multiply the spectrum for the geometrical spreading
coefficient and convert it to seismic moment units:

.. math::

   M^{p|s}_x(f) \equiv
   \mathcal{G}(r) \times
   \frac{4 \pi \rho_h^{1/2} \rho_r^{1/2} c_h^{5/2} c_r^{1/2}}
        {F R_{\Theta\Phi}}
   \times S^{p|s}_x(f)

Then the spectrum is converted in units of magnitude:

.. math::

   Y^{p|s}_x(f) \equiv
          \frac{2}{3} \times
          \left( \log_{10} M^{p|s}_x(f) - 9.1 \right)

And the final data vector :math:`Y^{p|s}(f)` is obtained by combining the
three components (e.g., Z, N, E) through the root sum of squares:

.. math::

   Y^{p|s}(f) =
          \sqrt{
             \left( Y^{p|s}_z(f) \right)^2 +
             \left( Y^{p|s}_n(f) \right)^2 +
             \left( Y^{p|s}_e(f) \right)^2
          }

The data vector is compared to the theoretical model :math:`M^{p|s}_{theo}(f)`
which incorporates the Brune's spectrum of the seismic moment and the inelastic
attenuation. This model is also converted in magnitude units:

.. math::

   Y^{p|s}(f) =
          \frac{2}{3}
          \left[ \log_{10} \left( M^{p|s}_{theo}(f) \right) - 9.1 \right]
          =
          \frac{2}{3}
          \left[ \log_{10} \left(
                    M_O \times
                    \frac{1}{1+\left(\frac{f}{f^{p|s}_c}\right)^2}
                    \times
                    e^{- \pi f t^*}
                    \right) - 9.1 \right] =

          =
          \frac{2}{3} \left( \log_{10} M_0 - 9.1 \right) +
          \frac{2}{3} \left[ \log_{10} \left(
                    \frac{1}{1+\left(\frac{f}{f^{p|s}_c}\right)^2} \right) +
                    \log_{10} \left( e^{- \pi f t^*} \right)
                    \right]


Finally coming to the following model used for the inversion:

.. math::

   Y^{p|s}(f) =
          M_w +
          \frac{2}{3} \left[ - \log_{10} \left(
                    1+\left(\frac{f}{f^{p|s}_c}\right)^2 \right) -
                    \pi \, f t^* \log_{10} e
                    \right]

where :math:`M_w \equiv \frac{2}{3} (\log_{10} M_0 - 9.1)`.


Inverted parameters
===================

The parameters determined from the spectral inversion are :math:`M_w`,
:math:`f^{p|s}_c` and :math:`t^*`.

The inversion is performed in moment magnitude :math:`M_w` units (logarithmic
amplitude). More details on the inversion method can be found in
:ref:`spectral_inversion:Spectral Inversion, Quality and Uncertainty`.


Other computed parameters
=========================

Starting from the inverted parameters :math:`M_0` ( :math:`M_w` ),
:math:`f^{p|s}_c`, :math:`t^*` and following the equations in :cite:t:`Madariaga2011`
and :cite:t:`Lancieri2012`, other quantities are computed for each station:

-  the static stress drop :math:`\Delta \sigma`
-  the source radius :math:`a`
-  the radiated energy :math:`E_r`
-  the apparent stress :math:`\sigma_a`
-  the quality factor :math:`Q_0` of P- or S-waves

As a bonus, local magnitude :math:`M_l` can be computed as well.

Event summaries (mean, weighted mean, percentiles) are computed from single
station estimates. For mean and weighted mean estimation, outliers are rejected
based on the `interquartile
range <https://en.wikipedia.org/wiki/Interquartile_range>`__ rule.


Source radius
-------------
The source radius is computed, assuming a circular rupture model, from the
corner frequency :math:`f^{p|s}_c` (:cite:t:`Kaneko2014`, equation 2):

.. math::

   a = k^{p|s} \frac{\beta_h}{f^{p|s}_c}

where :math:`\beta_h` is the shear wave speed at the hypocenter (in :math:`m/s`),
:math:`f^{p|s}_c` is the corner frequency (in :math:`Hz`) estimated from the
spectral inversion of P or S waves and :math:`k^{p|s}` is a constant which
depends on the source model.

:cite:t:`Brune1970` provides an expression for :math:`k^s` in the case of a
static circular crack (equation 31 in :cite:t:`Madariaga2011`):

.. math::

   k^s_{Brune} = 0.3724

:cite:t:`Kaneko2014` compiled a table including their own values for
:math:`k^p` and :math:`k^s` as well as values obtained from other authors.
The values are given as a function of the rupture velocity :math:`V_r` of a
dynamic circular crack (or a static crack, when :math:`V_r` is infinite):

.. list-table:: Table 1 of :cite:t:`Kaneko2014`
   :header-rows: 1

   * - :math:`V_r/\beta_h`
     - :math:`k^p_{K\&S}`
     - :math:`k^s_{K\&S}`
     - :math:`k^p_{Mada}`
     - :math:`k^s_{Mada}`
     - :math:`k^s_{Brune}`
     - :math:`k^p_{S\&H}`
     - :math:`k^s_{S\&H}`
   * - Infinite
     -
     -
     -
     -
     - 0.3724
     -
     -
   * - 0.9
     - 0.38
     - 0.26
     - 0.32
     - 0.21
     -
     - 0.42
     - 0.29
   * - 0.8
     - 0.35
     - 0.26
     -
     -
     -
     - 0.39
     - 0.28
   * - 0.7
     - 0.32
     - 0.26
     -
     -
     -
     - 0.36
     - 0.27
   * - 0.6
     - 0.30
     - 0.25
     -
     -
     -
     - 0.34
     - 0.27
   * - 0.5
     - 0.28
     - 0.22
     -
     -
     -
     - 0.31
     - 0.24

Where "K\&S" stands for :cite:t:`Kaneko2014`,
"Mada" for :cite:t:`Madariaga1976`, and "S\&H" for :cite:t:`Sato1973`.


Static stress drop
------------------
The static stress drop :math:`\Delta \sigma` is computed under the
assumption of a circular rupture of radius :math:`a`, as discussed in
:cite:t:`Madariaga2011` (equation 27):

.. math::

   \Delta \sigma =
   \frac{7}{16}
   \frac{M_0}{a^3}

where :math:`M_0` is the seismic moment (in :math:`N \cdot m`) and
:math:`a` is the source radius (in :math:`m`).


Radiated energy
---------------
The computation of the radiated energy :math:`E_r` starts with the integral of
the squared velocity amplitude spectrum:
:math:`[\dot{S}^{p|s}(f)]^2 = [ 2 \pi f S^{p|s}(f) ]^2`.

Following :cite:t:`Boatwright2002` (equation 1) and :cite:t:`Lancieri2012`
(equation 3), the P- or S-wave radiated energy is computed as:

.. math::

   \tilde{\tilde{E}}_r^{p|s} = 8 \pi \mathcal{G}^2(r) C^2 \rho_r c_r
            \int_{f_{min}}^{f_{max}} e^{2 \pi f t^*} [\dot{S}^{p|s}(f)]^2 df

where :math:`\mathcal{G}^2(r)` is the squared geometrical spreading coefficient
(see above), :math:`C` is a constant discussed below, :math:`\rho_r` and
:math:`c_r` are, respectively, the density and P- or S-wave velocity
at the receiver (their product is the seismic impedance), :math:`f_{min}` and
:math:`f_{max}` are the minimum and maximum frequency used to compute the
energy (see :ref:`configuration_file:Configuration File` for details on the
``Er_freq_range`` parameter), and the exponential term in the integrand is the
squared correction for anelastic attenuation.
The double tilde on top of :math:`\tilde{\tilde{E}}_r^{p|s}` means that the
radiated energy needs to be further corrected for noise and finite bandwidth
(see below).

The constant :math:`C` is defined in :cite:t:`Boatwright2002` (equation 2) as:

.. math::

   C = \frac{\left<R_{\Theta\Phi}\right>}{R_{\Theta\Phi} F}

where :math:`\left<R_{\Theta\Phi}\right>` is the root mean square P- or S-wave
radiation pattern computed on the focal sphere, :math:`R_{\Theta\Phi}` is the
radiation pattern coefficient for the given station, and :math:`F` is the
free surface amplification factor.
If a focal mechanism is not available, then it is assumed
:math:`R_{\Theta\Phi} = \left<R_{\Theta\Phi}\right>` and, hence,
:math:`C = 1/F`.
This assumption means that we rely on the averaging between measurements
of radiated energy at different stations, instead of precise measurements at a
single station.

Noise correction
++++++++++++++++
To account for low frequency noise, below the corner frequency, under the
hypothesis that energy is additive and that noise is stationary, we compute
a noise-corrected energy as:

.. math::

   \tilde{E}^{p|s}_r = \tilde{\tilde{E}}^{p|s}_r - \tilde{\tilde{E}}^{noise}_r

where the first term is the radiated energy computed from the P- or S-wave
spectrum and the second term is the radiated energy computed from the noise
spectrum. If the above difference is negative, then the measure is rejected,
since the noise is too large compared to the signal.

Finite bandwidth correction
+++++++++++++++++++++++++++
Next step is to correct the radiated energy for the missing energy above
:math:`f_{max}`, not taken into account in the integral of the squared velocity
amplitude spectrum (finite bandwidth correction).
Following :cite:t:`Lancieri2012` (equation 4), and :cite:t:`DiBona1988`, the
noise-corrected radiated energy is divided by the theoretical ratio :math:`R`
between the estimated radiated energy and the true radiated energy, defined as:

.. math::

  R = \frac{2}{\pi}
    \left[
      \arctan(f_{max}/f^{p|s}_c) -
      \frac{f_{max}/f^{p|s}_c}{1+(f_{max}/f^{p|s}_c)^2}
    \right]

where :math:`f_{max}` is the maximum frequency used to compute the energy
integral and :math:`f^{p|s}_c` is the P- or S-wave corner frequency.

The values of R range between 0 (for :math:`f_{max}/f^{p|s}_c \to 0`) and 1
(for :math:`f_{max}/f^{p|s}_c \to \infty`).

The corrected radiated energy for P- or S-waves is then:

.. math::

   E^{p|s}_r = \frac{\tilde{E}^{p|s}_r}{R}


Energy partition
++++++++++++++++
The final step is to account for the partition of energy between P and S waves.
Following :cite:t:`Boatwright1986` (equations 8 and 15) the ratio between the
radiated energy measured from S-waves and the radiated energy measured from
P-waves is:

.. math::

   \frac{E_r^s}{E_r^p} = 15.6

The final estimate of the radiated energy is then:

.. math::

   E_r = \left( 1 + 15.6 \right) E_r^p

or

.. math::

   E_r = \left( 1 + \frac{1}{15.6} \right) E_r^s

depending on whether the radiated energy is computed from P or S waves.

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
S-wave quality factor :math:`Q_0^{p|s}` using the following expression:

.. math::

   Q_0^{p|s} = \frac{tt^{p|s}(r)}{t^*}

where :math:`tt^{p|s}(r)` is the P- or S-wave travel time from source to
station and :math:`r` is the hypocentral distance.


Station Residuals
=================
Station-specific effects can be determined by running ``source_spec`` on several
events and computing the average of station residuals between observed and
inverted spectra. These averages are obtained through the command
``source_residuals``; the resulting residuals file can be used for a second run
of ``source_spec`` (see the ``residuals_filepath`` option in
:ref:`configuration_file:Configuration File`).