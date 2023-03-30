# -*- coding: utf-8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Build spectral objects.

:copyright:
    2012 Claudio Satriano <satriano@ipgp.fr>

    2013-2014 Claudio Satriano <satriano@ipgp.fr>,
              Emanuela Matrullo <matrullo@geologie.ens.fr>,
              Agnes Chounet <chounet@ipgp.fr>

    2015-2023 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
import logging
import numpy as np
import math
from scipy.integrate import cumtrapz
from scipy.interpolate import interp1d
from obspy.core import Stream
from sourcespec import spectrum
from sourcespec.ssp_setup import ssp_exit
from sourcespec.ssp_util import smooth, cosine_taper, moment_to_mag, get_vel
from sourcespec.ssp_process_traces import filter_trace
from sourcespec.ssp_correction import station_correction
from sourcespec.ssp_radiation_pattern import get_radiation_pattern_coefficient
logger = logging.getLogger(__name__.split('.')[-1])


class SpectrumIgnored(Exception):
    def __init__(self, message, reason):
        # Call the base class constructor with the parameters it needs
        super().__init__(message)
        # Now for your custom code...
        self.reason = reason


def _get_nint(config, trace):
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


def _time_integrate(config, trace):
    nint = _get_nint(config, trace)
    trace.detrend(type='constant')
    trace.detrend(type='linear')
    for _ in range(nint):
        trace.data = cumtrapz(trace.data) * trace.stats.delta
        trace.stats.npts -= 1
        filter_trace(config, trace)


def _frequency_integrate(config, spec):
    nint = _get_nint(config, spec)
    for _ in range(nint):
        spec.data /= (2 * math.pi * spec.get_freq())


def _cut_spectrum(config, spec):
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
    """
    if vertical_channel_codes is None:
        vertical_channel_codes = ['Z']
    spec_h = None
    for spec in spec_st.traces:
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
        # (only raidal and, optionally, vertical)
        if wave_type == 'SV' and channel[-1] == 'T':
            continue
        # only use vertical component for P
        if wave_type == 'P' and channel[-1] not in vertical_channel_codes:
            continue
        if spec_h is None:
            spec_h = spec.copy()
            spec_h.data = np.power(spec_h.data, 2)
            spec_h.data_log = np.power(spec_h.data_log, 2)
            spec_h.stats.channel = f'{code}H'
        else:
            spec_h.data += np.power(spec.data, 2)
            spec_h.data_log += np.power(spec.data_log, 2)
    if spec_h is not None:
        spec_h.data = np.sqrt(spec_h.data)
        spec_h.data_log = np.sqrt(spec_h.data_log)
    return spec_h


def _check_data_len(config, trace):
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
    if nzeros > npts/4:
        raise RuntimeError(
            f'{traceId}: Signal window is zero for more than 25%: '
            'skipping trace')


def _cut_signal_noise(config, trace):
    trace_signal = trace.copy()
    trace_noise = trace.copy()

    # Integrate in time domain, if required.
    # (otherwise frequency-domain integration is performed later)
    if config.time_domain_int:
        _time_integrate(config, trace_signal)
        _time_integrate(config, trace_noise)

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
            _recompute_time_window(trace, config.wave_type[0], npts)
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
            _recompute_time_window(
                trace, 'N', len(trace_signal), from_start=True)
        logger.warning(msg)

    # ...and zero pad to spectral_win_length
    if config.spectral_win_length is not None:
        trace_signal.trim(starttime=t1,
                          endtime=t1+config.spectral_win_length,
                          pad=True,
                          fill_value=0)
        trace_noise.trim(starttime=noise_t1,
                         endtime=noise_t1+config.spectral_win_length,
                         pad=True,
                         fill_value=0)

    return trace_signal, trace_noise


def _recompute_time_window(trace, wave_type, npts, from_start=False):
    """Recompute start or end time of signal or noise window,
    based on new number of points"""
    length = npts * trace.stats.delta
    if from_start:
        label, _ = trace.stats.arrivals[f'{wave_type}1']
        t1 = trace.stats.arrivals[f'{wave_type}2'][1] - length
        trace.stats.arrivals[f'{wave_type}1'] = (label, t1)
    else:
        label, _ = trace.stats.arrivals[f'{wave_type}2']
        t2 = trace.stats.arrivals[f'{wave_type}1'][1] + length
        trace.stats.arrivals[f'{wave_type}2'] = (label, t2)


def _check_noise_level(trace_signal, trace_noise, config):
    traceId = trace_signal.get_id()
    trace_signal_rms = ((trace_signal.data**2).sum())**0.5
    # Scale trace_noise_rms to length of signal window,
    # based on length of non-zero noise window
    try:
        scale_factor = float(len(trace_signal)) / len(trace_noise.data != 0)
    except ZeroDivisionError:
        scale_factor = 1
    trace_noise_rms = ((trace_noise.data**2 * scale_factor).sum())**0.5
    if trace_noise_rms/trace_signal_rms < 1e-6 and config.weighting == 'noise':
        # Skip trace if noise level is too low and if noise weighting is used
        raise RuntimeError(
            f'{traceId}: Noise level is too low or zero: '
            'station will be skipped')


def _geom_spread_r_power_n(hypo_dist_in_km, exponent):
    """rⁿ geometrical spreading coefficient."""
    dist = hypo_dist_in_km * 1e3
    return dist**exponent


def _geom_spread_boatwright(hypo_dist_in_km, cutoff_dist_in_km, freqs):
    """"
    Geometrical spreading coefficient from Boatwright et al. (2002), eq. 8.

    Except that we take the square root of eq. 8, since we correct amplitude
    and not energy.
    """
    dist = hypo_dist_in_km * 1e3
    cutoff_dist = cutoff_dist_in_km * 1e3
    if dist <= cutoff_dist:
        return dist
    else:
        return _boatwright_above_cutoff_dist(freqs, cutoff_dist, dist)


def _boatwright_above_cutoff_dist(freqs, cutoff_dist, dist):
    exponent = np.ones_like(freqs)
    low_freq = freqs <= 0.2
    mid_freq = np.logical_and(freqs > 0.2, freqs <= 0.25)
    high_freq = freqs >= 0.25
    exponent[low_freq] = 0.5
    exponent[mid_freq] = 0.5 + 2*np.log(5*freqs[mid_freq])
    exponent[high_freq] = 0.7
    return cutoff_dist * (dist/cutoff_dist)**exponent


def _geometrical_spreading_coefficient(config, spec):
    hypo_dist_in_km = spec.stats.hypo_dist
    if config.geom_spread_model == 'r_power_n':
        exponent = config.geom_spread_n_exponent
        return _geom_spread_r_power_n(hypo_dist_in_km, exponent)
    elif config.geom_spread_model == 'boatwright':
        cutoff_dist_in_km = config.geom_spread_cutoff_distance
        return _geom_spread_boatwright(
            hypo_dist_in_km, cutoff_dist_in_km, spec.get_freq())


# store log messages to avoid duplicates
velocity_log_messages = []


def _displacement_to_moment(stats, config):
    """
    Return the coefficient for converting displacement to seismic moment.

    From Aki&Richards,1980
    """
    phase = config.wave_type[0]
    if phase == 'P':
        v_hypo = config.hypo.vp
    elif phase == 'S':
        v_hypo = config.hypo.vs
    v_station = get_vel(
        stats.coords.longitude, stats.coords.latitude, -stats.coords.elevation,
        phase, config)
    specid = '.'.join((
        stats.network, stats.station, stats.location, stats.channel))
    msg = (
        f'{specid}: V{phase.lower()}_hypo: {v_hypo:.2f} km/s, '
        f'V{phase.lower()}_station: {v_station:.2f} km/s')
    global velocity_log_messages
    if msg not in velocity_log_messages:
        logger.info(msg)
        velocity_log_messages.append(msg)
    v_hypo *= 1000.
    v_station *= 1000.
    v3 = v_hypo**(5./2) * v_station**(1./2)
    rp = get_radiation_pattern_coefficient(stats, config)
    return 4 * math.pi * v3 * config.rho / (2 * rp)


def _smooth_spectrum(spec, smooth_width_decades=0.2):
    """Smooth spectrum in a log10-freq space."""
    # 1. Generate log10-spaced frequencies
    freq = spec.get_freq()
    _log_freq = np.log10(freq)
    # frequencies in logarithmic spacing
    log_df = _log_freq[-1] - _log_freq[-2]
    freq_logspace =\
        10**(np.arange(_log_freq[0], _log_freq[-1]+log_df, log_df))
    # 2. Reinterpolate data using log10 frequencies
    # make sure that extrapolation does not create negative values
    f = interp1d(freq, spec.data, fill_value='extrapolate')
    data_logspace = f(freq_logspace)
    data_logspace[data_logspace <= 0] = np.min(spec.data)
    # 3. Smooth log10-spaced data points
    npts = max(1, int(round(smooth_width_decades/log_df)))
    data_logspace = smooth(data_logspace, window_len=npts)
    # 4. Reinterpolate to linear frequencies
    # make sure that extrapolation does not create negative values
    f = interp1d(freq_logspace, data_logspace, fill_value='extrapolate')
    data = f(freq)
    data[data <= 0] = np.min(spec.data)
    spec.data = data
    # 5. Optimize the sampling rate of log spectrum,
    #    based on the width of the smoothing window
    # make sure that extrapolation does not create negative values
    log_df = smooth_width_decades/5
    freq_logspace =\
        10**(np.arange(_log_freq[0], _log_freq[-1]+log_df, log_df))
    spec.freq_log = freq_logspace
    data_logspace = f(freq_logspace)
    data_logspace[data_logspace <= 0] = np.min(spec.data)
    spec.data_log = data_logspace


def _build_spectrum(config, trace):
    spec = spectrum.do_spectrum(trace)
    stats = trace.stats
    spec.stats.instrtype = stats.instrtype
    spec.stats.coords = stats.coords
    spec.stats.hypo = stats.hypo
    spec.stats.hypo_dist = stats.hypo_dist
    spec.stats.epi_dist = stats.epi_dist
    spec.stats.ignore = stats.ignore
    spec.stats.travel_times = stats.travel_times
    # Integrate in frequency domain, if no time-domain
    # integration has been performed
    if not config.time_domain_int:
        _frequency_integrate(config, spec)
    # cut the spectrum
    spec = _cut_spectrum(config, spec)
    # correct geometrical spreading
    geom_spread = _geometrical_spreading_coefficient(config, spec)
    spec.data *= geom_spread
    # convert to seismic moment
    coeff = _displacement_to_moment(stats, config)
    spec.data *= coeff
    # store coeff to correct back data in displacement units
    # for radiated_energy()
    spec.coeff = coeff
    # smooth
    _smooth_spectrum(spec, config.spectral_smooth_width_decades)
    return spec


def _build_uniform_weight(spec):
    weight = spec.copy()
    weight.snratio = None
    weight.data = np.ones_like(weight.data)
    weight.data_log = np.ones_like(weight.data_log)
    return weight


def _build_weight_from_frequency(config, spec):
    weight = spec.copy()
    freq = weight.get_freq()
    weight.data = np.ones_like(weight.data)
    weight.data[freq <= config.f_weight] = config.weight
    weight.data /= np.max(weight.data)
    freq_log = weight.freq_log
    weight.data_log = np.ones_like(weight.data_log)
    weight.data_log[freq_log <= config.f_weight] = config.weight
    weight.data_log /= np.max(weight.data_log)
    return weight


def _build_weight_from_ratio(spec, specnoise, smooth_width_decades):
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
    cosine_taper(weight.data, weight.stats.delta/4, left_taper=True)
    # Make sure weight is positive
    weight.data[weight.data <= 0] = 0.001
    return weight


def _build_weight_from_noise(config, spec, specnoise):
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
    f = interp1d(weight.get_freq(), weight.data, fill_value='extrapolate')
    weight.data_log = f(weight.freq_log)
    weight.data_log /= np.max(weight.data_log)
    # Make sure weight is positive
    weight.data_log[weight.data_log <= 0] = 0.001
    return weight


def _build_weight_spectral_stream(config, spec_st, specnoise_st):
    """Build a stream of weights from a stream of spectra and a stream of
    noise spectra."""
    weight_st = Stream()
    spec_ids = {sp.id[:-1] for sp in spec_st if not sp.stats.ignore}
    for specid in spec_ids:
        try:
            spec_h = _select_spectra(spec_st, f'{specid}H')[0]
            specnoise_h = _select_spectra(specnoise_st, f'{specid}H')[0]
        except Exception:
            continue
        if config.weighting == 'noise':
            weight = _build_weight_from_noise(config, spec_h, specnoise_h)
        elif config.weighting == 'frequency':
            weight = _build_weight_from_frequency(config, spec_h)
        elif config.weighting == 'no_weight':
            weight = _build_uniform_weight(spec_h)
        weight_st.append(weight)
    return weight_st


def _select_spectra(spec_st, specid):
    """Select spectra from stream, based on specid."""
    network, station, location, channel = specid.split('.')
    channel = channel + '?'*(3-len(channel))
    spec_st_sel = spec_st.select(
        network=network, station=station, location=location, channel=channel)
    spec_st_sel = Stream(sp for sp in spec_st_sel if not sp.stats.ignore)
    return spec_st_sel


def _build_H(spec_st, specnoise_st=None, vertical_channel_codes=None,
             wave_type='S'):
    """
    Add to spec_st and specnoise_st the "H" component.

    H component is obtained from the modulus of all the available components.
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


def _check_spectral_sn_ratio(config, spec, specnoise):
    weight = _build_weight_from_noise(config, spec, specnoise)
    # if no noise window is available, snratio is not computed
    if weight.snratio is None:
        return
    if config.spectral_sn_freq_range is not None:
        sn_fmin, sn_fmax = config.spectral_sn_freq_range
        freqs = weight.get_freq()
        idx = np.where((sn_fmin <= freqs)*(freqs <= sn_fmax))
    else:
        idx = range(len(weight.snratio))
    spectral_snratio =\
        weight.snratio[idx].sum()/len(weight.snratio[idx])
    spec.stats.spectral_snratio = spectral_snratio
    spec_id = spec.get_id()
    logger.info(
        f'{spec_id}: spectral S/N: {spectral_snratio:.2f}')
    if config.spectral_sn_min:
        ssnmin = config.spectral_sn_min
        if spectral_snratio < ssnmin:
            msg = (
                f'{spec_id}: spectral S/N smaller than {ssnmin:.2f}: '
                'ignoring spectrum')
            reason = 'low spectral S/N'
            raise SpectrumIgnored(msg, reason)


def _ignore_spectrum(msg, trace, spec, specnoise):
    """Ignore spectrum. Set ignore flag and reason."""
    logger.warning(msg)
    trace.stats.ignore = True
    trace.stats.ignore_reason = msg.reason
    spec.stats.ignore = True
    spec.stats.ignore_reason = msg.reason
    specnoise.stats.ignore = True
    specnoise.stats.ignore_reason = msg.reason


def _build_signal_and_noise_streams(config, st):
    """Build signal and noise streams."""
    # remove traces with ignore flag
    traces = Stream([tr for tr in st if not tr.stats.ignore])
    signal_st = Stream()
    noise_st = Stream()
    for trace in sorted(traces, key=lambda tr: tr.id):
        try:
            _check_data_len(config, trace)
            trace_signal, trace_noise = _cut_signal_noise(config, trace)
            _check_noise_level(trace_signal, trace_noise, config)
            signal_st.append(trace_signal)
            noise_st.append(trace_noise)
        except RuntimeError as msg:
            # RuntimeError is for skipped traces
            logger.warning(msg)
            continue
    # trim components of the same instrument to the same number of samples
    for id in sorted({tr.id[:-1] for tr in signal_st}):
        st_sel = signal_st.select(id=f'{id}*') + noise_st.select(id=f'{id}*')
        all_npts = {tr.stats.npts for tr in st_sel}
        if len(all_npts) == 1:
            continue
        logger.warning(
            f'{id}: components have different window lengths. '
            'Trimming signal and noise windows to the shortest one')
        npts = min(all_npts)
        for tr in st_sel:
            if tr.stats.type == 'signal':
                tr.data = tr.data[:npts]
            elif tr.stats.type == 'noise':
                tr.data = tr.data[-npts:]
        for tr in st.select(id=f'{id}*'):
            _recompute_time_window(tr, config.wave_type[0], npts)
            _recompute_time_window(tr, 'N', npts, from_start=True)
    return signal_st, noise_st


def _build_signal_and_noise_spectral_streams(config, signal_st, noise_st):
    """Build signal and noise spectral streams."""
    spec_st = Stream()
    specnoise_st = Stream()
    for trace_signal in sorted(signal_st, key=lambda tr: tr.id):
        trace_noise = noise_st.select(id=trace_signal.id)[0]
        try:
            spec = _build_spectrum(config, trace_signal)
            specnoise = _build_spectrum(config, trace_noise)
            _check_spectral_sn_ratio(config, spec, specnoise)
        except RuntimeError as msg:
            # RuntimeError is for skipped spectra
            logger.warning(msg)
            continue
        except SpectrumIgnored as msg:
            _ignore_spectrum(msg, trace_signal, spec, specnoise)
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
        spec.data_log_mag = moment_to_mag(spec.data_log)
    for specnoise in specnoise_st:
        specnoise.data_mag = moment_to_mag(specnoise.data)
    # apply station correction if a residual file is specified in config
    spec_st = station_correction(spec_st, config)
    return spec_st, specnoise_st


def build_spectra(config, st):
    """
    Build spectra and the ``spec_st`` object.

    Computes P- or S-wave (displacement) spectra from
    accelerometers and velocimeters, uncorrected for anelastic attenuation,
    corrected for instrumental constants, normalized by geometrical spreading.
    """
    logger.info('Building spectra...')
    signal_st, noise_st = _build_signal_and_noise_streams(config, st)
    spec_st, specnoise_st = \
        _build_signal_and_noise_spectral_streams(config, signal_st, noise_st)
    weight_st = _build_weight_spectral_stream(config, spec_st, specnoise_st)
    logger.info('Building spectra: done')
    return spec_st, specnoise_st, weight_st
