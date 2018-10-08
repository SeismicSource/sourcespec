.. _source_spec:

###########
source_spec
###########

Earthquake source parameters from inversion of S-wave spectra.

:copyright:
    2012 Claudio Satriano <satriano@ipgp.fr>

    2013-2014 Claudio Satriano <satriano@ipgp.fr>,
              Emanuela Matrullo <matrullo@geologie.ens.fr>,
              Agnes Chounet <chounet@ipgp.fr>

    2015-2018 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement, Version 2.1
    (http://www.cecill.info/index.en.html)

Overview
========

``source_spec`` inverts the S-wave displacement spectra from
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
               {4 \pi \rho_h^{1/2} \rho_r^{1/2} \beta_h^{5/2} \beta_r^{1/2}}
          \times
          M_O
          \times
          \frac{1}{1+\left(\frac{f}{f_c}\right)^2}
          \times
          \exp \left( \frac{-\pi r f}{Q_O V_S} \right)

where
:math:`r` is the hypocentral distance;
:math:`R_{\Theta\Phi}` is the radiation pattern coefficient for S-waves;
:math:`\rho_h` and :math:`\rho_r` are the medium densities at the hypocenter
and at the receiver, respectively;
:math:`\beta_h` and :math:`\beta_r` are the S-wave velocities at the hypocenter
and at the receiver, respectively;
:math:`M_O` is the seismic moment;
:math:`f` is the frequency;
:math:`f_c` is the corner frequency;
:math:`V_S` is the average S-wave velocity along the wave propagation path;
:math:`Q_O` is the quality factor.



In ``source_spec``, the observed spectra :math:`S(f)` are converted in
moment magnitude :math:`M_w`.

The first step is to multiply the spectrum for the hypocentral distance
and convert them to seismic moment units:

.. math::

   M(f) \equiv
   r \times
   \frac{4 \pi \rho_h^{1/2} \rho_r^{1/2} \beta_h^{5/2} \beta_r^{1/2}}
        {2 R_{\Theta\Phi}}
   \times S(f) =
          M_O \times
          \frac{1}{1+\left(\frac{f}{f_c}\right)^2}
          \times
          \exp \left( \frac{-\pi r f}{Q_O V_S} \right)


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
                      \exp \left( \frac{-\pi r f}{Q_O V_S} \right)
                      \right) - 9.1 \right] =

            =
            \frac{2}{3} (\log_{10} M_0 - 9.1) +
            \frac{2}{3} \left[ \log_{10} \left(
                      \frac{1}{1+\left(\frac{f}{f_c}\right)^2} \right) +
                      \log_{10} \left(
                      \exp \left( \frac{-\pi r f}{Q_O V_S} \right) \right)
                      \right]


Finally coming to the following model used for the inversion:

.. math::

   Y_{data}(f) =
            M_w +
            \frac{2}{3} \left[ - \log_{10} \left(
                      1+\left(\frac{f}{f_c}\right)^2 \right) -
                      \pi \, f t^* \log_{10} e
                      \right]

Where :math:`M_w \equiv \frac{2}{3} (\log_{10} M_0 - 9.1)`
and :math:`t^* \equiv \frac{r}{Q_O V_S}`.

The parameters to determine are :math:`M_w`, :math:`f_c` and :math:`t^*`.
