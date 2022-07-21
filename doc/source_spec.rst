.. _source_spec:

###########
SourceSpec
###########

Earthquake source parameters from inversion of P- or S-wave spectra.

:copyright:
    2011-2022 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement, Version 2.1
    (http://www.cecill.info/index.en.html)

Overview
========

``source_spec`` inverts the P- or S-wave displacement spectra from
station recordings of a single event.

Spectral model
==============

The Fourier spectrum of the S-wave displacement in far field can be
modelled as the product of a source term (Brune model) and a
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

- Follwing Boatwright et al. (2002), eq. 8:

  - body wave spreading (:math:`\mathcal{G}(r) = r`) for hypocentral distances
    below a cutoff distance :math:`r_0`;
  - frequency dependent spreading for hypocentral distances above the
    cutoff distance :math:`r_0`.

More precisely, the expression derived from Boatwright et al. (2002) is:

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

In ``source_spec``, the observed spectrum :math:`S(f)` is converted into
moment magnitude units :math:`M_w`.

The first step is to multiply the spectrum for the geometrical spreading
coefficient and convert it to seismic moment units:

.. math::

   M(f) \equiv
   \mathcal{G}(r) \times
   \frac{4 \pi \rho_h^{1/2} \rho_r^{1/2} c_h^{5/2} c_r^{1/2}}
        {2 R_{\Theta\Phi}}
   \times S(f) =
          M_O \times
          \frac{1}{1+\left(\frac{f}{f_c}\right)^2}
          \times
          e^{- \pi f t^*}


Then the spectrum is converted in units of magnitude
(the :math:`Y_{data} (f)` vector used in the inversion):

.. math::

   Y_{data}(f) \equiv
            \frac{2}{3} \times
            \left( \log_{10} M(f) - 9.1 \right)

The data vector is compared to the teoretical model:

.. math::

   Y_{data}(f) =
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

   Y_{data}(f) =
            M_w +
            \frac{2}{3} \left[ - \log_{10} \left(
                      1+\left(\frac{f}{f_c}\right)^2 \right) -
                      \pi \, f t^* \log_{10} e
                      \right]

where :math:`M_w \equiv \frac{2}{3} (\log_{10} M_0 - 9.1)`.


Inversion procedure
===================

The parameters to determine are :math:`M_w`, :math:`f_c` and :math:`t^*`.

The retrieved attenuation parameter :math:`t^*` is then converted to the P- or
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
