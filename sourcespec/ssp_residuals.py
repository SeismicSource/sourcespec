# -*- coding: utf8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Spectral residual routine for sourcespec.

:copyright:
    2013-2014 Claudio Satriano <satriano@ipgp.fr>,
              Agnes Chounet <chounet@ipgp.fr>

    2015-2026 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
import os
import contextlib
import logging
from sourcespec.spectrum import SpectrumStream
from sourcespec.ssp_spectral_model import spectral_model
from sourcespec.ssp_util import (
    mag_to_moment, select_trace, set_spectrum_processing_info
)
logger = logging.getLogger(__name__.rsplit('.', maxsplit=1)[-1])


def spectral_residuals(config, spec_st, weight_st, sspec_output):
    """
    Compute spectral residuals with respect to an average spectral model.
    Saves a stream of residuals and weights to disk in HDF5 format.

    :param config: Configuration object
    :type config: :class:`~sourcespec.config.Config`
    :param spec_st: Stream of spectra
    :type spec_st: :class:`~sourcespec.spectrum.SpectrumStream`
    :param weight_st: Stream of spectral weights
    :type weight_st: :class:`~sourcespec.spectrum.SpectrumStream`
    :param sspec_output: Output of the source spectral parameter estimation
    :type sspec_output: :class:`~sourcespec.ssp_data_types.SourceSpecOutput`
    """
    # get reference summary values
    summary_values = sspec_output.reference_values()
    params_name = ('Mw', 'fc', 't_star')
    sourcepar_summary = {p: summary_values[p] for p in params_name}
    residuals = SpectrumStream()
    for station in {x.stats.station for x in spec_st}:
        spec_st_sel = spec_st.select(station=station)
        for spec in spec_st_sel:
            if spec.stats.channel[-1] != 'H':
                continue
            if spec.stats.ignore:
                continue
            xdata = spec.freq
            synth_mean_mag = spectral_model(xdata, **sourcepar_summary,
                                            n=config.source_falloff_power)
            res = spec.copy()
            # remove existing data
            del res.data
            del res.data_logspaced
            del res.data_mag
            del res.data_mag_logspaced
            # Residuals are computed in magnitude units,
            # but res.data (moment units) has to be set before res.data_mag
            _res_mag = spec.data_mag - synth_mean_mag
            res.data = mag_to_moment(_res_mag)
            res.data_mag = _res_mag
            set_spectrum_processing_info(res, config)
            res.stats.data_type = 'residual'
            # Store spectral weight as a separate spectrum object
            weight_spec = select_trace(
                weight_st, spec.id, spec.stats.instrtype)
            if weight_spec is not None:
                weight_copy = weight_spec.copy()
                set_spectrum_processing_info(weight_copy, config)
                residuals.append(weight_copy)
            # remove t_star_model from stats before saving,
            # as it might be a function, which is unsupported by HDF5
            with contextlib.suppress(KeyError):
                del res.stats['t_star_model']
            # remove Q_model, as it might be None
            with contextlib.suppress(KeyError):
                del res.stats['Q_model']
            residuals.append(res)
    # Save residuals to HDF5 file
    evid = config.event.event_id
    res_file = os.path.join(config.options.outdir, f'{evid}.residuals.hdf5')
    residuals.write(res_file, format='HDF5')
    logger.info(f'Spectral residuals saved to: {res_file}')
