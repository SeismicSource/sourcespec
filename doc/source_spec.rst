##############
source_spec.py
##############

Overview
========

``source_spec.py`` inverts the S-wave displacement spectra from
station recordings of a single event.

Spectral model
==============

The Fourier spectrum of the S-wave displacement in far field can be
modelled as the product of a source term (Brune model) and a
propagation term (geometric and anelastic attenuation of body waves):

.. math::

   S(f) = M_O \times \frac{2 R_{\Theta\Phi}}{4 \pi \rho \beta^3}
          \times
          \left( \frac{1}{1+\left(\frac{f}{f_c}\right)^2} \right)
          \times
          \left( \exp \left( \frac{-\pi r f}{Q_O V_S} \right)
                 \frac{1}{r} \right)

where :math:`f` is the freqeuncy, :math:`r` is the hypocentral
distance, :math:`M_O` is the seismic moment, :math:`f_c` is the
corner frequency; :math:`R_{\Theta\Phi}` is the radiation pattern
coefficient for S-waves, :math:`\rho` is the average density of the
medium, :math:`\beta` and :math:`V_S` are the S-wave speed at the
source and the average S-wave speed along the wave propagation path,
respectively; finally, :math:`Q_O` is the quality factor.

.. automodule:: source_spec
   :members:
