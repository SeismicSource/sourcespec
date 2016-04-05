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

    2015-2016 Claudio Satriano <satriano@ipgp.fr>
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

   S(f) = M_O \times \frac{2 R_{\Theta\Phi}}{4 \pi \rho \beta^3}
          \times
          \frac{1}{1+\left(\frac{f}{f_c}\right)^2}
          \times
          \left[ \exp \left( \frac{-\pi r f}{Q_O V_S} \right)
                 \frac{1}{r} \right]

where :math:`f` is the freqeuncy, :math:`r` is the hypocentral
distance, :math:`M_O` is the seismic moment, :math:`f_c` is the
corner frequency; :math:`R_{\Theta\Phi}` is the radiation pattern
coefficient for S-waves, :math:`\rho` is the average density of the
medium, :math:`\beta` and :math:`V_S` are the S-wave speed at the
source and the average S-wave speed along the wave propagation path,
respectively; finally, :math:`Q_O` is the quality factor.



In source_spec, the observed spectra :math:`S(f)` are converted in
moment magnitude :math:`Mw`.

The first step is to multiply the spectrum for the hypocentral distance
and convert them to seismic moment units:

.. math::

   r \times
   \frac{4 \pi \rho \beta^3}{2 R_{\Theta\Phi}} \times
   S(f) =
          M_O \times
          \frac{1}{1+\left(\frac{f}{f_c}\right)^2}
          \times
          \exp \left( \frac{-\pi r f}{Q_O V_S} \right)


Then the spectrum is converted in unities of magnitude
(the :math:`Y_{data}` vector used in the inversion):

.. math::

   Y_{data} =
            \frac{2}{3} \times
            \left[ \log_{10} \left(
                      r \times
                      \frac{4 \pi \rho \beta^3}{2 R_{\Theta\Phi}} \times
                      S(f)
                      \right) - 9.1 \right]


   Y_{data} =
            \frac{2}{3}
            \left[ \log_{10} \left(
                      M_O \times
                      \frac{1}{1+\left(\frac{f}{f_c}\right)^2}
                      \times
                      \exp \left( \frac{-\pi r f}{Q_O V_S} \right)
                      \right) - 9.1 \right]


   Y_{data} =
            \frac{2}{3} (\log_{10} M_0 - 9.1) +
            \frac{2}{3} \left[ \log_{10} \left(
                      \frac{1}{1+\left(\frac{f}{f_c}\right)^2} \right) +
                      \log_{10} \left(
                      \exp \left( \frac{-\pi r f}{Q_O V_S} \right) \right)
                      \right]


Finally coming to the following model used for the inversion:

.. math::

   Y_{data} =
            M_w +
            \frac{2}{3} \left[ - \log_{10} \left(
                      1+\left(\frac{f}{f_c}\right)^2 \right) -
                      \pi \, f t^* \log_{10} e
                      \right]

Where :math:`Mw \equiv \frac{2}{3} (\log_{10} M_0 - 9.1)`
and :math:`t^* \equiv \frac{r}{Q_O V_S}`
