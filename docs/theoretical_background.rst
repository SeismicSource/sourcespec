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

The same thing is done for a noise time window: noise spectrum is used to
compute spectral signal-to-noise ratio (and possibily reject low SNR spectra)
and, optionnaly, to weight the spectral inversion.

.. figure:: imgs/example_trace.svg
  :alt: Example trace plot
  :width: 600

  Example three-component trace plot (in velocity), showing noise and S-wave
  windows.

The component spectra are combined through the root sum of squares:

.. math::

    S(f) = \sqrt{S^2_z(f) + S^2_n(f) + S^2_e(f)}

where :math:`f` is the frequency and :math:`S_x(f)``` is the P- or S-wave
spectrum for component :math:`x`.

.. figure:: imgs/example_spectrum.svg
  :alt: Example spectrum plot
  :width: 600

  Example displacement spectrum for noise and S-wave, including inversion
  results.


Spectral model
==============

The Fourier amplitude spectrum of the S-wave displacement in far field can be
modelled as the product of a source term :cite:p:`Brune1970` and a
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
one of the following ways:

- :math:`\mathcal{G}(r) = r^n`: :math:`n` can be any positive number.
  :math:`n=1` (default value) is the theoretical value for a body wave in a
  homogeneous full-space;
  :math:`n=0.5` is the theoretical value for a surface wave in a homogeneous
  half-space.

- Follwing :cite:t:`Boatwright2002`, eq. 8:

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
    0.5 + 2 \log (5f)  &  0.20 < f < 0.25 Hz\\
    0.7  &  f \ge 0.25 Hz\\
  \end{cases}

Note that here we use the square root of eq. 8 in Boatwright et al. (2002),
since we correct the spectral amplitude and not the energy.


Building spectra
================

In ``source_spec``, the observed spectrum of component :math:`x`,
:math:`S_x(f)` is converted into moment magnitude units :math:`M_w`.

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

The data vector is compared to the teoretical model:

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


Inversion procedure
===================

The parameters to determine are :math:`M_w`, :math:`f_c` and :math:`t^*`.

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

Starting from the inverted parameters :math:`M_0` ( :math:`M_w` ),
:math:`fc`, :math:`t^*` and following the equations in :cite:t:`Madariaga2011`,
other quantities are computed for each station:

-  the Brune stress drop
-  the source radius
-  the quality factor :math:`Q_0` of P- or S-waves

Finally, the radiated energy :math:`E_r` can be mesured from the
displacement spectra, following the approach described in
:cite:t:`Lancieri2012`.

As a bonus, local magnitude :math:`M_l` can be computed as well.

Event averages are computed from single station estimates. Outliers are
rejected based on the `interquartile
range <https://en.wikipedia.org/wiki/Interquartile_range>`__ rule.


Attenuation
-----------
The retrieved attenuation parameter :math:`t^*` is converted to the P- or
S-wave quality factor :math:`Q_0^{[P|S]}` using the following expression:

.. math::

   Q_0^{[P|S]} = \frac{tt_{[P|S]}(r)}{t^*}

where :math:`tt_{[P|S]}(r)` is the P- or S-wave travel time from source to
station and :math:`r` is the hypocentral distance.

Station-specific effects can be determined by running ``source_spec`` on several
events and computing the average of station residuals between observed and
inverted spectra. These averages are obtained through the command
``source_residuals``; the resulting residuals file can be used for a second run
of ``source_spec`` (see the ``residuals_filepath`` option in
:ref:`configuration_file:Configuration File`).
