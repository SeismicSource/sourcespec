# -*- coding: utf8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Direct spectral modelling.

:copyright:
    2013-2014 Claudio Satriano <satriano@ipgp.fr>,
              Agnes Chounet <chounet@ipgp.fr>

    2015-2026 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""


def make_synth(config, spec_st, trace_spec=None):
    """
    Make synthetic spectra.
    """
    # pylint: disable=import-outside-toplevel
    # Lazy-import modules for speed
    import math
    import numpy as np
    from copy import deepcopy
    from .spectrum import Spectrum
    from .ssp_spectral_model import spectral_model, objective_func
    from .ssp_util import mag_to_moment, moment_to_mag
    fdelta = 0.01
    fmin = config.options.fmin
    fmax = config.options.fmax + fdelta

    residuals = []
    for n, (fc, mag, Mo, t_star, alpha) in enumerate(zip(
            config.options.fc, config.options.mag, config.options.Mo,
            config.options.t_star, config.options.alpha)):
        spec = Spectrum()
        if trace_spec:
            spec.stats = deepcopy(trace_spec.stats)
            spec.freq = trace_spec.freq.copy()
        else:
            spec.freq = np.arange(fmin, fmax, fdelta)

        if math.isnan(Mo):
            Mo = mag_to_moment(mag)
        else:
            mag = moment_to_mag(Mo)

        spec.stats.station = f'S{n:02d}'
        spec.stats.instrtype = 'Synth'
        spec.stats.channel = f'{spec.stats.channel[:-1]}S'
        spec.stats.par = {
            'Mw': mag, 'fc': fc, 't_star': t_star, 'alpha': alpha}

        freq = spec.freq
        freq_logspaced = trace_spec.freq_logspaced
        spec.freq_logspaced = freq_logspaced
        # spec.data and spec.data_logspaced have to be set before
        # spec.data_mag and spec.data_mag_logspaced, so that the Spectrum
        # object can check for consistency
        data_mag = spectral_model(freq, mag, fc, t_star, alpha)
        spec.data = mag_to_moment(data_mag)
        spec.data_mag = data_mag
        data_mag_logspaced = spectral_model(
            freq_logspaced, mag, fc, t_star, alpha)
        spec.data_logspaced = mag_to_moment(data_mag_logspaced)
        spec.data_mag_logspaced = data_mag_logspaced
        spec_st.append(spec)

        if trace_spec:
            objective_func2 = objective_func(freq, trace_spec.data_mag,
                                             np.ones_like(trace_spec.data_mag))
            print(Mo, mag, fc, t_star, alpha,
                  objective_func2((mag, fc, t_star, alpha)))
            residuals.append([Mo, mag, fc, t_star,
                             objective_func2((mag, fc, t_star, alpha))])


def main():
    """
    Main function.
    """
    # pylint: disable=import-outside-toplevel
    # Lazy-import modules for speed
    from .ssp_parse_arguments import parse_args
    options = parse_args(progname='source_model')
    from .config import config, configure_cli
    from .ssp_setup import ssp_exit
    plot_show = bool(options.plot)
    conf_overrides = {
        'plot_show': plot_show,
        'plot_save': False,
        'html_report': False
    }
    configure_cli(
        options, progname='source_model', config_overrides=conf_overrides)
    from .spectrum import SpectrumStream
    from .ssp_read_traces import read_traces
    from .ssp_process_traces import process_traces
    from .ssp_build_spectra import build_spectra

    # We don't use weighting in source_model
    config.weighting = 'no_weight'
    if len(config.options.trace_path) > 0:
        st = read_traces()
        # Deconvolve, filter, cut traces:
        proc_st = process_traces(config, st)
        # Build spectra (amplitude in magnitude units)
        spec_st, _specnoise_st, _weight_st = build_spectra(config, proc_st)
        if len(spec_st) == 0:
            ssp_exit()
        # We keep just horizontal component:
        spec_st = SpectrumStream(
            [x for x in spec_st if x.stats.channel[-1] == 'H'])

        for trace_spec in spec_st:
            orientation = trace_spec.stats.channel[-1]
            if orientation != 'H':
                continue
            make_synth(config, spec_st, trace_spec)
    else:
        spec_st = SpectrumStream()
        make_synth(config, spec_st)

    from .ssp_plot_spectra import plot_spectra
    from .ssp_plot_traces import plot_traces
    plot_traces(config, proc_st, ncols=2, block=False)
    plot_spectra(config, spec_st, ncols=1, stack_plots=True)


if __name__ == '__main__':
    main()
