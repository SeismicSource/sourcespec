# -*- coding: utf-8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Build spectral objects.

:copyright:
    2012 Claudio Satriano <satriano@ipgp.fr>

    2013-2014 Claudio Satriano <satriano@ipgp.fr>,
              Emanuela Matrullo <matrullo@geologie.ens.fr>,
              Agnes Chounet <chounet@ipgp.fr>

    2015-2025 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
import logging
import math
import numpy as np
# cumtrapz was deprecated in scipy 1.12.0 and removed in 1.14.0
try:
    from scipy.integrate import cumtrapz
except ImportError:
    from scipy.integrate import cumulative_trapezoid as cumtrapz
from scipy.interpolate import interp1d
from obspy.core import Stream
from .config import config
from .spectrum import Spectrum, SpectrumStream
from .ssp_setup import ssp_exit
from .ssp_util import (
    smooth, cosine_taper, moment_to_mag, MediumProperties)
from .ssp_geom_spreading import (
    geom_spread_r_power_n, geom_spread_r_power_n_segmented,
    geom_spread_boatwright, geom_spread_teleseismic)
from .ssp_process_traces import filter_trace
from .ssp_correction import station_correction
from .ssp_radiation_pattern import get_radiation_pattern_coefficient
logger = logging.getLogger(__name__.rsplit('.', maxsplit=1)[-1])


class SpectrumIgnored(Exception):
    """Spectrum ignored exception"""
    def __init__(self, message, reason):
        # Call the base class constructor with the parameters it needs
        super().__init__(message)
        # Now for your custom code...
        self.reason = reason


def _get_nint(trace):
    """
    Return the number of integrations to be performed on the trace.

    :param trace: ObsPy trace object
    :type trace: :class:`obspy.core.trace.Trace`

    :return: Number of integrations
    :rtype: int

    :raises ValueError: If the instrument type is unknown
    """
    if config.trace_units == 'auto':
        instrtype = trace.stats.instrtype
    else:
        instrtype = config.trace_units
    if instrtype == 'acc':
        nint = 2
    elif instrtype in ['broadb', 'shortp', 'vel']:
        nint = 1
    elif instrtype == 'disp':
        nint = 0
    else:
        raise ValueError
    return nint


def _time_integrate(trace):
    """
    Integrate the trace in time domain.

    Trace is filtered after each integration.

    :param trace: ObsPy trace object
    :type trace: :class:`obspy.core.trace.Trace
    """
    nint = _get_nint(trace)
    trace.detrend(type='constant')
    trace.detrend(type='linear')
    for _ in range(nint):
        trace.data = cumtrapz(trace.data) * trace.stats.delta
        trace.stats.npts -= 1
        filter_trace(trace)


def _frequency_integrate(spec):
    """
    Integrate the spectrum in frequency domain.

    :param spec: Spectrum object
    :type spec: :class:`~sourcespec.spectrum.Spectrum`
    """
    nint = _get_nint(spec)
    for _ in range(nint):
        spec.data /= (2 * math.pi * spec.freq)


def _cut_spectrum(spec):
    """
    Cut the spectrum to the frequency range specified in the configuration.

    :param spec: Spectrum object
    :type spec: :class:`~sourcespec.spectrum.Spectrum`

    :return: Cut spectrum
    :rtype: :class:`~sourcespec.spectrum.Spectrum`

    :raises RuntimeError: If the instrument type is unknown
    """
    # see if there is a station-specific frequency range
    station = spec.stats.station
    try:
        freq1 = float(config[f'freq1_{station}'])
        freq2 = float(config[f'freq2_{station}'])
    except KeyError:
        try:
            instrtype = spec.stats.instrtype
            freq1 = float(config[f'freq1_{instrtype}'])
            freq2 = float(config[f'freq2_{instrtype}'])
        except KeyError as e:
            raise RuntimeError(
                f'{spec.id}: Unknown instrument type: {instrtype}: '
                'skipping spectrum'
            ) from e
    return spec.slice(freq1, freq2)


def _compute_h(spec_st, code, vertical_channel_codes=None, wave_type='S'):
    """
    Compute the component 'H' from geometric mean of the stream components.

    (which can also be all three components)

    :param spec_st: Stream of spectra
    :type spec_st: :class:`~sourcespec.spectrum.SpectrumStream`
    :param code: Band+instrument code
    :type code: str
    :param vertical_channel_codes: List of vertical channel codes
        (default: ['Z'] if None)
    :type vertical_channel_codes: list
    :param wave_type: Wave type ('S' or 'P', default: 'S')
    :type wave_type: str

    :return: Spectrum object with the 'H' component
    :rtype: :class:`~sourcespec.spectrum.Spectrum`
    """
    if vertical_channel_codes is None:
        vertical_channel_codes = ['Z']
    spec_h = None
    for spec in spec_st:
        # this avoids taking a component from co-located station:
        # ('code' is band+instrument code)
        channel = spec.stats.channel
        if channel[:-1] != code:
            continue
        # avoid reusing a previous 'H' channel
        if channel[-1] == 'H':
            continue
        # only use transverse component for SH
        if wave_type == 'SH' and channel[-1] != 'T':
            continue
        # do not use transverse component for SV
        if wave_type == 'SV' and channel[-1] == 'T':
            continue
        if spec_h is None:
            spec_h = spec.copy()
            spec_h.data = np.power(spec_h.data, 2)
            spec_h.data_logspaced = np.power(spec_h.data_logspaced, 2)
            spec_h.stats.channel = f'{code}H'
        else:
            spec_h.data += np.power(spec.data, 2)
            spec_h.data_logspaced += np.power(spec.data_logspaced, 2)
    if spec_h is not None:
        spec_h.data = np.sqrt(spec_h.data)
        spec_h.data_logspaced = np.sqrt(spec_h.data_logspaced)
    return spec_h


def _check_data_len(trace):
    """
    Check if data length is sufficient.

    :param trace: ObsPy trace object
    :type trace: :class:`obspy.core.trace.Trace`

    :raises RuntimeError: If data length is insufficient
    """
    traceId = trace.get_id()
    trace_cut = trace.copy()
    if config.wave_type[0] == 'S':
        t1 = trace.stats.arrivals['S1'][1]
        t2 = trace.stats.arrivals['S2'][1]
    elif config.wave_type[0] == 'P':
        t1 = trace.stats.arrivals['P1'][1]
        t2 = trace.stats.arrivals['P2'][1]
    trace_cut.trim(starttime=t1, endtime=t2, pad=True, fill_value=0)
    npts = len(trace_cut.data)
    if npts == 0:
        raise RuntimeError(
            f'{traceId}: No data for the selected cut interval: '
            'skipping trace')
    nzeros = len(np.where(trace_cut.data == 0)[0])
    if nzeros > npts / 4:
        raise RuntimeError(
            f'{traceId}: Signal window is zero for more than 25%: '
            'skipping trace')


def _cut_signal_noise(trace):
    """
    Cut signal and noise windows from trace.

    :param trace: ObsPy trace object
    :type trace: :class:`obspy.core.trace.Trace`

    :return: Signal and noise traces
    :rtype: tuple of :class:`obspy.core.trace.Trace
    """
    trace_signal = trace.copy()
    trace_noise = trace.copy()

    # Integrate in time domain, if required.
    # (otherwise frequency-domain integration is performed later)
    if config.time_domain_int:
        _time_integrate(trace_signal)
        _time_integrate(trace_noise)

    # trim...
    if config.wave_type[0] == 'S':
        t1 = trace.stats.arrivals['S1'][1]
        t2 = trace.stats.arrivals['S2'][1]
    elif config.wave_type[0] == 'P':
        t1 = trace.stats.arrivals['P1'][1]
        t2 = trace.stats.arrivals['P2'][1]
    trace_signal.trim(starttime=t1, endtime=t2, pad=True, fill_value=0)
    trace_signal.stats.type = 'signal'
    # Noise time window for weighting function:
    noise_t1 = trace.stats.arrivals['N1'][1]
    noise_t2 = trace.stats.arrivals['N2'][1]
    trace_noise.trim(
        starttime=noise_t1, endtime=noise_t2, pad=True, fill_value=0)
    trace_noise.stats.type = 'noise'
    # ...taper...
    cosine_taper(trace_signal.data, width=config.taper_halfwidth)
    cosine_taper(trace_noise.data, width=config.taper_halfwidth)

    # Be sure that both traces have same length:
    npts = min(len(trace_signal), len(trace_noise))
    if len(trace_signal) <= len(trace_noise):
        # truncate noise to signal length
        trace_noise.data = trace_noise.data[:npts]
    else:
        tr_noise_id = trace_noise.get_id()[:-1]
        if config.weighting == 'noise':
            msg = f'{tr_noise_id}: truncating signal window to noise length!'
            trace_signal.data = trace_signal.data[:npts]
            # recompute signal window end time, so that it's plotted correctly
            _recompute_time_window(
                trace, config.wave_type[0], npts, keep='start')
        else:
            msg = f'{tr_noise_id}: zero-padding noise window to signal length'
            # Notes:
            # 1. no risk of ringing, as noise has been tapered
            # 2. we use np.pad instead of obspy trim method
            #    to avoid potential rounding errors with start/end times
            # 3. we just pad at the end for simplicity (no changes
            #    in trace headers required) and it is identical to
            #    symmetric padding (tested)
            pad_len = len(trace_signal) - len(trace_noise)
            trace_noise.data = np.pad(trace_noise.data, (0, pad_len))
            # recompute noise window start time, so that it's plotted correctly
            _recompute_time_window(trace, 'N', len(trace_signal), keep='end')
        logger.warning(msg)
    return trace_signal, trace_noise


def _recompute_time_window(trace, wave_type, npts, keep='start'):
    """
    Recompute start or end time of signal or noise window,
    based on new number of points.

    :param trace: ObsPy trace object
    :type trace: :class:`obspy.core.trace.Trace`
    :param wave_type: Wave type ('P', 'S', 'N')
    :type wave_type: str
    :param npts: Number of points
    :type npts: int
    :param keep: Keep 'start' or 'end' of window (default: 'start')
    :type keep: str

    :raises ValueError: If keep is not 'start' or 'end'
    """
    length = npts * trace.stats.delta
    if keep == 'end':
        label, _ = trace.stats.arrivals[f'{wave_type}1']
        t1 = trace.stats.arrivals[f'{wave_type}2'][1] - length
        trace.stats.arrivals[f'{wave_type}1'] = (label, t1)
    elif keep == 'start':
        label, _ = trace.stats.arrivals[f'{wave_type}2']
        t2 = trace.stats.arrivals[f'{wave_type}1'][1] + length
        trace.stats.arrivals[f'{wave_type}2'] = (label, t2)
    else:
        raise ValueError('keep must be "start" or "end"')


def _check_noise_level(trace_signal, trace_noise):
    """
    Check noise level.

    :param trace_signal: Signal trace
    :type trace_signal: :class:`obspy.core.trace.Trace`
    :param trace_noise: Noise trace
    :type trace_noise: :class:`obspy.core.trace.Trace`

    :raises RuntimeError: If noise level is too low
    """
    traceId = trace_signal.get_id()
    trace_signal_rms = ((trace_signal.data**2).sum())**0.5
    # Scale trace_noise_rms to length of signal window,
    # based on length of non-zero noise window
    try:
        scale_factor = float(len(trace_signal)) / len(trace_noise.data != 0)
    except ZeroDivisionError:
        scale_factor = 1
    trace_noise_rms = ((trace_noise.data**2 * scale_factor).sum())**0.5
    if (
        trace_noise_rms / trace_signal_rms < 1e-6 and
        config.weighting == 'noise'
    ):
        # Skip trace if noise level is too low and if noise weighting is used
        raise RuntimeError(
            f'{traceId}: Noise level is too low or zero: '
            'station will be skipped')


def _geometrical_spreading_coefficient(spec):
    """
    Return the geometrical spreading coefficient for the given spectrum.

    :param spec: Spectrum object
    :type spec: :class:`~sourcespec.spectrum.Spectrum`

    :return: Geometrical spreading coefficient
    :rtype: float

    :raises ValueError: If the geometrical spreading model is unknown
    """
    hypo_dist_in_km = spec.stats.hypo_dist
    epi_dist_in_km = spec.stats.epi_dist
    # set geometrical spreading distance to a very large value if it is None
    geom_spread_min_dist = config.geom_spread_min_teleseismic_distance or 1e99
    geom_spread_model =\
        'teleseismic' if epi_dist_in_km >= geom_spread_min_dist\
        else config.geom_spread_model
    logger.info(f'{spec.id}: geometrical spreading model: {geom_spread_model}')
    if geom_spread_model == 'teleseismic':
        angular_distance = spec.stats.gcarc
        source_depth_in_km = spec.stats.event.hypocenter.depth.value_in_km
        station_depth_in_km = -spec.stats.coords.elevation
        phase = config.wave_type[0]
        return geom_spread_teleseismic(
            angular_distance, source_depth_in_km, station_depth_in_km, phase)
    if config.geom_spread_model == 'r_power_n':
        exponent = config.geom_spread_n_exponent
        return geom_spread_r_power_n(hypo_dist_in_km, exponent)
    if config.geom_spread_model == 'r_power_n_segmented':
        exponents = config.geom_spread_n_exponents
        hinge_distances = config.geom_spread_n_distances
        if len(hinge_distances) != len(exponents):
            raise ValueError(
                f'The number of exponents must be equal to the number of '
                f'hinge distances. You provided {len(exponents)} exponents '
                f'and {len(hinge_distances)} hinge distances'
            )
        return geom_spread_r_power_n_segmented(hypo_dist_in_km, exponents,
                                               hinge_distances)
    if config.geom_spread_model == 'boatwright':
        cutoff_dist_in_km = config.geom_spread_cutoff_distance
        return geom_spread_boatwright(
            hypo_dist_in_km, cutoff_dist_in_km, spec.freq)
    raise ValueError(
        f'Unknown geometrical spreading model: {config.geom_spread_model}')


# store log messages to avoid duplicates
PROPERTY_LOG_MESSAGES = []


def _displacement_to_moment(stats):
    """
    Return the coefficient for converting displacement to seismic moment.

    From Aki&Richards,1980

    :param stats: Stats object
    :type stats: :class:`~sourcespec.spectrum.AttributeDict`
    """
    phase = config.wave_type[0]
    lon = stats.coords.longitude
    lat = stats.coords.latitude
    depth = -stats.coords.elevation
    medium_properties = MediumProperties(lon, lat, depth, config)
    depth_string = medium_properties.to_string('station depth', depth)
    v_name = f'v{phase.lower()}'
    v_source = config.event.hypocenter[v_name]
    v_source_string = medium_properties.to_string(f'{v_name}_source', v_source)
    v_station = medium_properties.get(mproperty=v_name, where='stations')
    stats.v_station = v_station
    stats.v_station_type = phase
    v_station_string = medium_properties.to_string(
        f'{v_name}_station', v_station)
    rho_source = config.event.hypocenter.rho
    rho_source_string = medium_properties.to_string('rho_source', rho_source)
    rho_station = medium_properties.get(mproperty='rho', where='stations')
    stats.rho_station = rho_station
    rho_station_string = medium_properties.to_string(
        'rho_station', rho_station)
    specid = '.'.join((
        stats.network, stats.station, stats.location, stats.channel))
    msg = (
        f'{specid}: {depth_string}, '
        f'{v_source_string}, {v_station_string}'
    )
    if msg not in PROPERTY_LOG_MESSAGES:
        logger.info(msg)
        PROPERTY_LOG_MESSAGES.append(msg)
    msg = (
        f'{specid}: {depth_string}, '
        f'{rho_source_string}, {rho_station_string}'
    )
    if msg not in PROPERTY_LOG_MESSAGES:
        logger.info(msg)
        PROPERTY_LOG_MESSAGES.append(msg)
    v_source *= 1000.
    v_station *= 1000.
    v3 = v_source**(5. / 2) * v_station**(1. / 2)
    rho = rho_source**0.5 * rho_station**0.5
    fsa = config.free_surface_amplification
    return 4 * math.pi * v3 * rho / (fsa * stats.radiation_pattern)


def _smooth_spectrum(spec, smooth_width_decades=0.2):
    """
    Smooth spectrum in a log10-freq space.

    :param spec: Spectrum object
    :type spec: :class:`~sourcespec.spectrum.Spectrum`
    :param smooth_width_decades: Width of the smoothing window in decades
        (default: 0.2)
    :type smooth_width_decades: float
    """
    # 1. Generate log10-spaced frequencies
    freq = spec.freq
    _log_freq = np.log10(freq)
    # frequencies in logarithmic spacing
    log_df = _log_freq[-1] - _log_freq[-2]
    freq_logspaced =\
        10**(np.arange(_log_freq[0], _log_freq[-1] + log_df, log_df))
    # 2. Reinterpolate data using log10 frequencies
    # make sure that extrapolation does not create negative values
    f = interp1d(freq, spec.data, fill_value='extrapolate')
    data_logspaced = f(freq_logspaced)
    data_logspaced[data_logspaced <= 0] = np.min(spec.data)
    # 3. Smooth log10-spaced data points
    npts = max(1, int(round(smooth_width_decades / log_df)))
    data_logspaced = smooth(data_logspaced, window_len=npts)
    # 4. Reinterpolate to linear frequencies
    # make sure that extrapolation does not create negative values
    f = interp1d(freq_logspaced, data_logspaced, fill_value='extrapolate')
    data = f(freq)
    data[data <= 0] = np.min(spec.data)
    spec.data = data
    # 5. Optimize the sampling rate of log spectrum,
    #    based on the width of the smoothing window
    # make sure that extrapolation does not create negative values
    log_df = smooth_width_decades / 5
    freq_logspaced =\
        10**(np.arange(_log_freq[0], _log_freq[-1] + log_df, log_df))
    spec.freq_logspaced = freq_logspaced
    data_logspaced = f(freq_logspaced)
    data_logspaced[data_logspaced <= 0] = np.min(spec.data)
    spec.data_logspaced = data_logspaced


def _build_spectrum(trace):
    """
    Build a spectrum from a trace.

    :param trace: ObsPy trace object
    :type trace: :class:`obspy.core.trace.Trace`

    :return: Spectrum object
    :rtype: :class:`~sourcespec.spectrum.Spectrum`

    :raises RuntimeError: If an error occurs while building the spectrum
    """
    spec = Spectrum(obspy_trace=trace)
    spec.stats.instrtype = trace.stats.instrtype
    spec.stats.coords = trace.stats.coords
    spec.stats.event = trace.stats.event
    spec.stats.hypo_dist = trace.stats.hypo_dist
    spec.stats.epi_dist = trace.stats.epi_dist
    spec.stats.gcarc = trace.stats.gcarc
    spec.stats.azimuth = trace.stats.azimuth
    spec.stats.travel_times = trace.stats.travel_times
    spec.stats.takeoff_angles = trace.stats.takeoff_angles
    spec.stats.ignore = trace.stats.ignore
    # Integrate in frequency domain, if no time-domain
    # integration has been performed
    if not config.time_domain_int:
        _frequency_integrate(spec)
    # cut the spectrum
    spec = _cut_spectrum(spec)
    # correct geometrical spreading
    try:
        geom_spread = _geometrical_spreading_coefficient(spec)
    except Exception as e:
        raise RuntimeError(
            f'{spec.id}: Error computing geometrical spreading: '
            f'skipping spectrum\n{str(e)}'
        ) from e
    spec.data *= geom_spread
    # store the radiation pattern coefficient in the spectrum stats
    spec.stats.radiation_pattern =\
        get_radiation_pattern_coefficient(spec.stats)
    # convert to seismic moment
    coeff = _displacement_to_moment(spec.stats)
    spec.data *= coeff
    # store coeff to correct back data in displacement units
    # for radiated_energy()
    spec.stats.coeff = coeff
    # smooth
    try:
        _smooth_spectrum(spec, config.spectral_smooth_width_decades)
    except ValueError as e:
        raise RuntimeError(
            f'{spec.id}: Error smoothing spectrum: '
            f'skipping spectrum\n{str(e)}'
        ) from e
    return spec


def _build_uniform_weight(spec):
    """
    Build a uniform spectral weight.

    :param spec: Spectrum object
    :type spec: :class:`~sourcespec.spectrum.Spectrum`

    :return: Uniform weight
    :rtype: :class:`~sourcespec.spectrum.Spectrum`
    """
    weight = spec.copy()
    weight.snratio = None
    weight.data = np.ones_like(weight.data)
    weight.data_logspaced = np.ones_like(weight.data_logspaced)
    return weight


def _build_weight_from_frequency(spec):
    """
    Build spectral weights from frequency.

    :param spec: Spectrum object
    :type spec: :class:`~sourcespec.spectrum.Spectrum`

    :return: Spectral weights
    :rtype: :class:`~sourcespec.spectrum.Spectrum`
    """
    weight = spec.copy()
    freq = weight.freq
    weight.data = np.ones_like(weight.data)
    weight.data[freq <= config.f_weight] = config.weight
    weight.data /= np.max(weight.data)
    freq_logspaced = weight.freq_logspaced
    weight.data_logspaced = np.ones_like(weight.data_logspaced)
    weight.data_logspaced[freq_logspaced <= config.f_weight] = config.weight
    weight.data_logspaced /= np.max(weight.data_logspaced)
    return weight


def _build_weight_from_inv_frequency(spec, power=0.25):
    """
    Build spectral weights from inverse frequency (raised to a power < 1)

    :param spec: Spectrum object
    :type spec: :class:`~sourcespec.spectrum.Spectrum`
    :param power: Power to which the inverse frequency is raised
        (default: 0.25)
    :type power: float

    :return: Spectral weights
    :rtype: :class:`~sourcespec.spectrum.Spectrum`

    :raises ValueError: If power is >= 1
    """
    if power >= 1:
        raise ValueError('pow must be < 1')
    # Note: weight.data is used for plotting,
    #       weight.data_logspaced for actual weighting
    weight = spec.copy()
    freq = weight.freq
    weight.data *= 0
    # Limit non-zero weights to fmin/fmax from spectral_snratio if available
    snr_fmin = getattr(spec.stats, 'spectral_snratio_fmin', None)
    snr_fmax = getattr(spec.stats, 'spectral_snratio_fmax', None)
    i0 = np.where(freq >= snr_fmin)[0][0] if snr_fmin else 0
    i1 = np.where(freq <= snr_fmax)[0][-1] if snr_fmax else len(freq) - 1
    # Build weights as if frequencies always start from 0.25 Hz
    # to obtain similar curves regardless of fmin
    # and to avoid too much weight for very low frequencies
    weight.data[i0: i1 + 1] = 1. / (freq[i0: i1 + 1] - freq[i0] + 0.25)**power
    weight.data /= np.max(weight.data)
    freq_logspaced = weight.freq_logspaced
    weight.data_logspaced *= 0
    i0 = np.where(freq_logspaced >= snr_fmin)[0][0] if snr_fmin else 0
    i1 = np.where(freq_logspaced <= snr_fmax)[0][-1] if snr_fmax\
        else len(freq_logspaced) - 1
    weight.data_logspaced[i0: i1 + 1] =\
        1. / (freq_logspaced[i0: i1 + 1] - freq_logspaced[i0] + 0.25)**power
    weight.data_logspaced /= np.max(weight.data_logspaced)
    return weight


def _build_weight_from_ratio(spec, specnoise, smooth_width_decades):
    """
    Build spectral weights from the ratio of signal to noise.

    :param spec: signal spectrum
    :type spec: :class:`~sourcespec.spectrum.Spectrum`
    :param specnoise: noise spectrum
    :type specnoise: :class:`~sourcespec.spectrum.Spectrum`
    :param smooth_width_decades: Width of the smoothing window in decades
    :type smooth_width_decades: float

    :return: Spectral weights
    :rtype: :class:`~sourcespec.spectrum.Spectrum`
    """
    weight = spec.copy()
    weight.data /= specnoise.data
    # save signal-to-noise ratio before log10, smoothing, and normalization
    weight.snratio = weight.data.copy()
    # The inversion is done in magnitude units,
    # so let's take log10 of weight
    weight.data = np.log10(weight.data)
    # Weight spectrum is smoothed once more
    _smooth_spectrum(weight, smooth_width_decades)
    weight.data /= np.max(weight.data)
    # slightly taper weight at low frequencies, to avoid overestimating
    # weight at low frequencies, in cases where noise is underestimated
    cosine_taper(
        weight.data,
        min(0.25, weight.stats.delta / 4),
        left_taper=True)
    # Zero out weight below 0.2 Hz, so that it does not affect the fit
    weight.data[weight.data <= 0.2] = 1e-9
    return weight


def _build_weight_from_noise(spec, specnoise):
    """
    Build spectral weights from signal to noise ratio, if available.
    Otherwise, build uniform weights.

    :param spec: Signal spectrum
    :type spec: :class:`~sourcespec.spectrum.Spectrum`
    :param specnoise: Noise spectrum
    :type specnoise: :class:`~sourcespec.spectrum.Spectrum`

    :return: Spectral weights
    :rtype: :class:`~sourcespec.spectrum.Spectrum`
    """
    if specnoise is None or np.all(specnoise.data == 0):
        spec_id = spec.get_id()[:-1]
        logger.warning(
            f'{spec_id}: No available noise window: '
            'a uniform weight will be applied')
        weight = _build_uniform_weight(spec)
    else:
        weight = _build_weight_from_ratio(
            spec, specnoise, config.spectral_smooth_width_decades)
    # interpolate to log-frequencies
    f = interp1d(weight.freq, weight.data, fill_value='extrapolate')
    weight.data_logspaced = f(weight.freq_logspaced)
    weight.data_logspaced /= np.max(weight.data_logspaced)
    # Make sure weight is positive
    weight.data_logspaced[weight.data_logspaced <= 0] = 0.001
    return weight


def _build_weight_spectral_stream(spec_st, specnoise_st):
    """
    Build a stream of weights from a stream of spectra and a stream of
    noise spectra.

    :param spec_st: Stream of signal spectra
    :type spec_st: :class:`~sourcespec.spectrum.SpectrumStream`
    :param specnoise_st: Stream of noise spectra
    :type specnoise_st: :class:`~sourcespec.spectrum.SpectrumStream`

    :return: Stream of weights
    :rtype: :class:`~sourcespec.spectrum.SpectrumStream`
    """
    weight_st = SpectrumStream()
    spec_ids = {sp.id[:-1] for sp in spec_st if not sp.stats.ignore}
    for specid in spec_ids:
        try:
            spec_h = _select_spectra(spec_st, f'{specid}H')[0]
            specnoise_h = _select_spectra(specnoise_st, f'{specid}H')[0]
        except Exception:
            continue
        if config.weighting == 'noise':
            weight = _build_weight_from_noise(spec_h, specnoise_h)
        elif config.weighting == 'frequency':
            weight = _build_weight_from_frequency(spec_h)
        elif config.weighting == 'inv_frequency':
            weight = _build_weight_from_inv_frequency(spec_h)
        elif config.weighting == 'no_weight':
            weight = _build_uniform_weight(spec_h)
        weight_st.append(weight)
    return weight_st


def _select_spectra(spec_st, specid):
    """
    Select spectra from stream, based on specid.

    :param spec_st: Stream of spectra
    :type spec_st: :class:`~sourcespec.spectrum.SpectrumStream`
    :param specid: Spectrum ID
    :type specid: str

    :return: Stream of selected spectra
    :rtype: :class:`~sourcespec.spectrum.SpectrumStream`
    """
    network, station, location, channel = specid.split('.')
    channel = channel + '?' * (3 - len(channel))
    spec_st_sel = spec_st.select(
        network=network, station=station, location=location, channel=channel)
    spec_st_sel = SpectrumStream(
        sp for sp in spec_st_sel if not sp.stats.ignore)
    return spec_st_sel


def _build_H(spec_st, specnoise_st=None, vertical_channel_codes=None,
             wave_type='S'):
    """
    Add to spec_st and specnoise_st the "H" component.

    H component is obtained from the modulus of all the available components.

    :param spec_st: Stream of spectra
    :type spec_st: :class:`~sourcespec.spectrum.SpectrumStream`
    :param specnoise_st: Stream of noise spectra
    :type specnoise_st: :class:`~sourcespec.spectrum.SpectrumStream`

    :param vertical_channel_codes: List of vertical channel codes
        (default: ['Z'] if None)
    :type vertical_channel_codes: list
    :param wave_type: Wave type ('S' or 'P', default: 'S')
    :type wave_type: str
    """
    if vertical_channel_codes is None:
        vertical_channel_codes = ['Z']
    spec_ids = {sp.id[:-1] for sp in spec_st if not sp.stats.ignore}
    for specid in spec_ids:
        spec_st_sel = _select_spectra(spec_st, specid)
        specnoise_st_sel = _select_spectra(specnoise_st, specid)
        # 'code' is band+instrument code
        for code in {x.stats.channel[:-1] for x in spec_st_sel}:
            spec_h = _compute_h(
                spec_st_sel, code, vertical_channel_codes, wave_type)
            if spec_h is not None:
                spec_st.append(spec_h)
            specnoise_h = _compute_h(
                specnoise_st_sel, code, vertical_channel_codes, wave_type)
            if specnoise_h is not None:
                specnoise_st.append(specnoise_h)


def _check_spectral_sn_ratio(spec, specnoise):
    """
    Check spectral signal-to-noise ratio.

    :param spec: Signal spectrum
    :type spec: :class:`~sourcespec.spectrum.Spectrum`
    :param specnoise: Noise spectrum
    :type specnoise: :class:`~sourcespec.spectrum.Spectrum`

    :raises SpectrumIgnored: If the spectrum is to be ignored
    """
    spec_id = spec.get_id()
    weight = _build_weight_from_noise(spec, specnoise)
    freqs = weight.freq
    # if no noise window is available, snratio is not computed
    if weight.snratio is None:
        spec.stats.spectral_snratio = None
        return
    if config.spectral_sn_freq_range is not None:
        sn_fmin, sn_fmax = config.spectral_sn_freq_range[:2]
        valid_freqs_idx = np.where((sn_fmin <= freqs) * (freqs <= sn_fmax))
        valid_snratio = weight.snratio[valid_freqs_idx]
    else:
        valid_snratio = weight.snratio
    if len(valid_snratio) == 0:
        spec.stats.spectral_snratio = np.nan
        msg = (
            f'{spec_id}: no valid frequency to compute spectral S/N: '
            'ignoring spectrum'
        )
        reason = 'no valid frequency'
        raise SpectrumIgnored(msg, reason)
    spectral_snratio = valid_snratio.mean()
    spec.stats.spectral_snratio = spectral_snratio
    # Save frequency range where SNR > 3 so it can be used for building weights
    # Note: not sure if we could use config.spectral_sn_min here instead of 3
    snr_valid_freqs = freqs[weight.snratio >= 3]
    try:
        spec.stats.spectral_snratio_fmin = snr_valid_freqs[0]
        spec.stats.spectral_snratio_fmax = snr_valid_freqs[-1]
    except IndexError:
        spec.stats.spectral_snratio_fmin = None
        spec.stats.spectral_snratio_fmax = None
    logger.info(f'{spec_id}: average spectral S/N: {spectral_snratio:.2f}')
    ssnmin = config.spectral_sn_min or -np.inf
    if spectral_snratio < ssnmin:
        msg = (
            f'{spec_id}: spectral S/N smaller than {ssnmin:.2f}: '
            'ignoring spectrum')
        reason = 'low spectral S/N'
        raise SpectrumIgnored(msg, reason)


def _ignore_spectrum(msg, spec, specnoise):
    """
    Ignore spectrum. Set ignore flag and reason.

    :param msg: Ignore message
    :type msg: str
    :param spec: Signal spectrum to ignore
    :type spec: :class:`~sourcespec.spectrum.Spectrum`
    :param specnoise: Noise spectrum object to ignore
    """
    logger.warning(msg)
    spec.stats.ignore = True
    spec.stats.ignore_reason = msg.reason
    specnoise.stats.ignore = True
    specnoise.stats.ignore_reason = msg.reason


def _ignore_trace(msg, trace):
    """
    Ignore trace. Set ignore flag and reason.

    :param msg: Ignore message
    :type msg: str
    :param trace: Trace to ignore
    :type trace: :class:`obspy.core.trace.Trace`
    """
    # NOTE: no logger.warning here, because it is already done in
    # _ignore_spectrum()
    trace.stats.ignore = True
    trace.stats.ignore_reason = msg.reason


def _build_signal_and_noise_streams(st):
    """
    Build signal and noise streams.

    :param st: ObsPy Stream object
    :type st: :class:`obspy.core.stream.Stream`

    :return: Signal and noise streams
    :rtype: tuple of :class:`obspy.core.stream.Stream`
    """
    # remove traces with ignore flag
    traces = Stream([tr for tr in st if not tr.stats.ignore])
    signal_st = Stream()
    noise_st = Stream()
    for trace in sorted(traces, key=lambda tr: tr.id):
        try:
            _check_data_len(trace)
            trace_signal, trace_noise = _cut_signal_noise(trace)
            _check_noise_level(trace_signal, trace_noise)
            signal_st.append(trace_signal)
            noise_st.append(trace_noise)
        except RuntimeError as msg:
            # RuntimeError is for skipped traces
            logger.warning(msg)
            continue
    return signal_st, noise_st


def _trim_components(signal_st, noise_st, original_st):
    """
    Trim components of the same instrument to the same number of samples.

    Recompute time window of the signal and noise traces for correct plotting.

    :param signal_st: Stream of signal traces
    :type signal_st: :class:`obspy.core
    :param noise_st: Stream of noise traces
    :type noise_st: :class:`obspy.core
    :param original_st: Original stream
    :type original_st: :class:`obspy.core.stream.Stream`
    """
    for traceid in sorted({tr.id[:-1] for tr in signal_st}):
        st_sel = signal_st.select(id=f'{traceid}*') +\
            noise_st.select(id=f'{traceid}*')
        all_npts = {tr.stats.npts for tr in st_sel}
        if len(all_npts) == 1:
            continue
        logger.warning(
            f'{traceid}: components have different window lengths. '
            'Trimming signal and noise windows to the shortest one')
        npts = min(all_npts)
        for tr in st_sel:
            if tr.stats.type == 'signal':
                tr.data = tr.data[:npts]
            elif tr.stats.type == 'noise':
                tr.data = tr.data[-npts:]
        for tr in original_st.select(id=f'{traceid}*'):
            _recompute_time_window(tr, config.wave_type[0], npts, keep='start')
            _recompute_time_window(tr, 'N', npts, keep='end')


def _build_signal_and_noise_spectral_streams(signal_st, noise_st, original_st):
    """
    Build signal and noise spectral streams.

    Note: original_st is only used to keep track of ignored traces.

    :param signal_st: Stream of signal traces
    :type signal_st: :class:`obspy.core
    :param noise_st: Stream of noise traces
    :type noise_st: :class:`obspy.core
    :param original_st: Original stream
    :type original_st: :class:`obspy.core

    :return: Signal and noise spectral streams
    :rtype: tuple of :class:`~sourcespec.spectrum.SpectrumStream`
    """
    spec_st = SpectrumStream()
    specnoise_st = SpectrumStream()
    for trace_signal in sorted(signal_st, key=lambda tr: tr.id):
        trace_noise = noise_st.select(id=trace_signal.id)[0]
        try:
            spec = _build_spectrum(trace_signal)
            specnoise = _build_spectrum(trace_noise)
            _check_spectral_sn_ratio(spec, specnoise)
        except RuntimeError as msg:
            # RuntimeError is for skipped spectra
            logger.warning(msg)
            continue
        except SpectrumIgnored as msg:
            _ignore_spectrum(msg, spec, specnoise)
            trace_original = original_st.select(id=trace_signal.id)[0]
            _ignore_trace(msg, trace_original)
        spec_st.append(spec)
        specnoise_st.append(specnoise)
    if not spec_st:
        logger.error('No spectra left! Exiting.')
        ssp_exit()
    # build H component
    _build_H(
        spec_st, specnoise_st, config.vertical_channel_codes, config.wave_type)
    # convert the spectral amplitudes to moment magnitude
    for spec in spec_st:
        spec.data_mag = moment_to_mag(spec.data)
        spec.data_mag_logspaced = moment_to_mag(spec.data_logspaced)
    for specnoise in specnoise_st:
        specnoise.data_mag = moment_to_mag(specnoise.data)
    # apply station correction if a residual file is specified in config
    spec_st = station_correction(spec_st)
    return spec_st, specnoise_st


def _zero_pad(trace):
    """
    Zero-pad trace to spectral_win_length

    :param trace: ObsPy trace object
    :type trace: :class:`obspy.core.trace.Trace`
    """
    spec_win_len = config.spectral_win_length
    if spec_win_len is None:
        return
    if trace.stats.type == 'signal':
        wtype = config.wave_type[0]
        t1 = trace.stats.arrivals[f'{wtype}1'][1]
    elif trace.stats.type == 'noise':
        t1 = trace.stats.arrivals['N1'][1]
    trace.trim(starttime=t1, endtime=t1 + spec_win_len, pad=True, fill_value=0)


def build_spectra(st):
    """
    Build spectra and the ``spec_st`` object.

    Computes P- or S-wave (displacement) spectra from
    accelerometers and velocimeters, uncorrected for anelastic attenuation,
    corrected for instrumental constants, normalized by geometrical spreading.

    :param st: ObsPy Stream object
    :type st: :class:`obspy.core.stream.Stream`

    :return: spectra, noise spectra, and weights
    :rtype: tuple of :class:`~sourcespec.spectrum.Spectrum`
    """
    wave_type = config.wave_type
    logger.info(f'Building {wave_type}-wave spectra...')
    signal_st, noise_st = _build_signal_and_noise_streams(st)
    _trim_components(signal_st, noise_st, st)
    for trace in signal_st + noise_st:
        _zero_pad(trace)
    spec_st, specnoise_st = _build_signal_and_noise_spectral_streams(
        signal_st, noise_st, st)
    weight_st = _build_weight_spectral_stream(spec_st, specnoise_st)
    logger.info(f'Building {wave_type}-wave spectra: done')
    logger.info('---------------------------------------------------')
    return spec_st, specnoise_st, weight_st
