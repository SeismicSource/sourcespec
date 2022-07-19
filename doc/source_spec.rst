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
          \frac{1}{r}
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

- :math:`r` is the hypocentral distance (:math:`\frac{1}{r}` is geometric
  attenuation);
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


In ``source_spec``, the observed spectra :math:`S(f)` are converted into
moment magnitude :math:`M_w`.

The first step is to multiply the spectrum for the hypocentral distance
and convert them to seismic moment units:

.. math::

   M(f) \equiv
   r \times
   \frac{4 \pi \rho_h^{1/2} \rho_r^{1/2} c_h^{5/2} c_r^{1/2}}
        {2 R_{\Theta\Phi}}
   \times S(f) =
          M_O \times
          \frac{1}{1+\left(\frac{f}{f_c}\right)^2}
          \times
          e^{- \pi f t^*}


Then the spectrum is converted in unities of magnitude
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

The parameters to determine are :math:`M_w`, :math:`f_c` and :math:`t^*`.

The retrieved attenuation parameter :math:`t^*` is then converted to the P- or
S-wave quality factor :math:`Q_0^{[P|S]}` using the following expression:

.. math::

   Q_0^{[P|S]} = \frac{tt_{[P|S]}(r)}{t^*}

where :math:`tt_{[P|S]}(r)` is the P- or S-wave travel time from source to
station.

Station-specific effects can be determined by running ``source_spec`` on several
events and computing the average of station residuals between observed and
inverted spectra. These averages are obtained through the command
``source_residuals``; the resulting residuals file can be used for a second run
of ``source_spec`` (see the ``residuals_filepath`` option in
:ref:`configuration_file:Configuration File`).
