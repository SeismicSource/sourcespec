# -*- coding: utf8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Spectral residual routine for sourcespec.

:copyright:
    2013-2014 Claudio Satriano <satriano@ipgp.fr>,
              Agnes Chounet <chounet@ipgp.fr>

    2015-2024 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
import os
import logging
from sourcespec._version import get_versions
from sourcespec.spectrum import SpectrumStream
from sourcespec.ssp_spectral_model import spectral_model
from sourcespec.ssp_util import mag_to_moment
logger = logging.getLogger(__name__.rsplit('.', maxsplit=1)[-1])


def spectral_residuals(config, spec_st, sspec_output):
    """
    Compute spectral residuals with respect to an average spectral model.
    Saves a stream of residuals to disk in HDF5 format.

    :param config: Configuration object
    :type config: :class:`~sourcespec.config.Config`
    :param spec_st: Stream of spectra
    :type spec_st: :class:`~sourcespec.spectrum.SpectrumStream`
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
            xdata = spec.freq
            synth_mean_mag = spectral_model(xdata, **sourcepar_summary)
            res = spec.copy()
            res.data_mag = spec.data_mag - synth_mean_mag
            res.data = mag_to_moment(res.data_mag)
            res.stats.software = 'SourceSpec'
            res.stats.software_version = get_versions()['version']
            residuals.append(res)
    # Save residuals to HDF5 file
    evid = config.event.event_id
    res_file = os.path.join(config.options.outdir, f'{evid}.residuals.hdf5')
    residuals.write(res_file, format='HDF5')
    logger.info(f'Spectral residuals saved to: {res_file}')
