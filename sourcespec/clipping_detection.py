# -*- coding: utf-8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Trace clipping score based on kernel density estimation.

:copyright:
    2023 Claudio Satriano <satriano@ipgp.fr>,
         Kris Vanneste <kris.vanneste@oma.be>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
import numpy as np
from scipy.stats import gaussian_kde
from scipy.signal import find_peaks


def _get_baseline(signal):
    """Get the signal baseline using a Savitzky-Golay filter."""
    from scipy.signal import savgol_filter
    # search for baselines with characteristic length of at least 1/3
    # of the signal
    wlen = len(signal) // 3
    if wlen % 2 == 0:
        wlen += 1
    return savgol_filter(signal, wlen, 3)


def _preprocess(trace, remove_linear_trend=False, remove_baseline=False):
    """Preprocess a trace for clipping detection."""
    trace = trace.copy()
    if remove_linear_trend:
        trace.detrend('linear')
    trace.detrend('demean')
    max_abs_data = np.max(np.abs(trace.data)) * 1.10  # 10% margin
    # Remove trace baseline if requested
    trace_baseline = None
    if remove_baseline:
        trace_baseline = _get_baseline(trace.data)
        trace.data -= trace_baseline
    return trace, trace_baseline, max_abs_data


def _get_kernel_density(trace, min_data, max_data, num_kde_bins=101):
    """Compute the kernel density of a trace."""
    # Compute gaussian kernel density using default bandwidth (Scott’s Rule)
    kde = gaussian_kde(trace.data)
    density_points = np.linspace(min_data, max_data, num_kde_bins)
    density = kde.pdf(density_points)
    density /= np.max(density)
    return density, density_points


def _get_histogram(trace):
    """Compute the histogram of a trace."""
    # Compute data histogram with a number of bins equal to 0.5% of data points
    # or 31, whichever is larger
    nbins = max(31, int(len(trace.data)*0.005))
    if nbins % 2 == 0:
        nbins += 1
    counts, bins = np.histogram(trace.data, bins=nbins)
    counts = counts/np.max(counts)
    bin_width = bins[1] - bins[0]
    return counts, bins, bin_width


def _get_distance_weight(density_points, order=2, min_weight=1, max_weight=5):
    """Compute the distance weight for the kernel density."""
    dist_weight = np.abs(density_points)**order
    dist_weight *= (max_weight-min_weight)/dist_weight.max()
    dist_weight += min_weight
    return dist_weight


def is_clipped(trace, sensitivity=3, clipping_percentile=10, debug=False):
    """
    Check if a trace is clipped, based on kernel density estimation.

    Kernel density estimation is used to find the peaks of the histogram of
    the trace data points. The peaks are then weighted by their distance from
    the trace average (which should be the most common value).
    The peaks with the highest weight are then checked for prominence,
    which is a measure of how much higher the peak is than the surrounding
    data. The prominence threshold is determined by the sensitivity parameter.
    If a peak is found in the range defined by lower_clip_bound, the trace
    is considered clipped or distorted.

    Parameters
    ----------
    trace : obspy.core.trace.Trace
        Trace to check.
    sensitivity : int
        Sensitivity level, from 1 (least sensitive) to 5 (most sensitive).
        (default: 3)
    clipping_percentile : float, between 0 and 100
        Percentile of trace amplitude range (expressed as percentage) to check
        for clipping. Default is 10, which means that the 10% highest and
        lowest values of the trace amplitude will be checked for clipping.
        A value of 0 means that no clipping check will be performed.
    debug : bool
        If True, plot trace, samples histogram and kernel density.

    Returns
    -------
    bool
        True if trace is clipped, False otherwise.
    """
    sensitivity = int(sensitivity)
    if sensitivity < 1 or sensitivity > 5:
        raise ValueError('sensitivity must be between 1 and 5')
    if clipping_percentile < 0 or clipping_percentile > 100:
        raise ValueError('clipping_percentile must be between 0 and 100')
    if clipping_percentile == 0:
        return False
    lower_clip_bound = 100 - clipping_percentile
    num_kde_bins = 101
    num_edge_bins = int(np.ceil(
        (num_kde_bins/2.) * (100 - lower_clip_bound) / 100.))
    trace, _, _ = _preprocess(
        trace, remove_linear_trend=False, remove_baseline=False)
    min_data = np.min(trace.data)
    max_data = np.max(trace.data)
    density, density_points = _get_kernel_density(
        trace, min_data, max_data, num_kde_bins)
    # Distance weight, parabolic, between 1 and 5
    dist_weight = _get_distance_weight(
        density_points, order=2, min_weight=1, max_weight=5)
    density_weight = density*dist_weight
    # Add 1 bin at start/end with value equal to min. of 5 first/last
    # bins to ensure local maxima at start/end are recognized as peaks
    density_weight = np.hstack([[density_weight[:num_edge_bins].min()],
                                density_weight,
                                [density_weight[-num_edge_bins:].min()]])
    # find peaks with minimum prominence based on clipping sensitivity
    min_prominence = [0.5, 0.3, 0.2, 0.15, 0.125]
    peaks, props = find_peaks(
        density_weight,
        prominence=min_prominence[sensitivity-1]
    )
    # Remove start/end bins again
    peaks = peaks[(peaks > 0) & (peaks < (num_kde_bins + 1))]
    peaks -= 1
    density_weight = density_weight[1:-1]
    npeaks = len(peaks)
    # Clipped peaks are peaks in the edge bins
    npeaks_clipped = np.sum(
        peak < num_edge_bins or peak > (num_kde_bins - 1 - num_edge_bins)
        for peak in peaks
    )
    properties = {
        'npeaks': npeaks,
        'npeaks_clipped': npeaks_clipped,
        'peaks': peaks,
        'prominences': props['prominences'],
    }
    # If there is at least a peak in the edge bins,
    # then the signal is probably clipped or distorted
    trace_clipped = npeaks_clipped > 0
    if debug:
        _plot_clipping_analysis(
            trace, max_data, density_points, density, density_weight,
            peaks=peaks, num_edge_bins=num_edge_bins,
            num_kde_bins=num_kde_bins, trace_clipped=trace_clipped)
    return trace_clipped, properties


def get_clipping_score(trace, remove_baseline=False, debug=False):
    """
    Compute a trace clipping score based on kernel density estimation.

    The algorithm is based on the following steps:

    1.  The trace is detrended and demeaned. Optionally, the trace baseline
        can be removed.

    2.  A kernel density estimation is performed on the trace amplitude values.

    3.  Two weighted kernel density functions are computed:

          - a full weighted kernel density, where the kernel density is
            weighted by the distance from the zero mean value, using a 8th
            order power function between 1 and 100.
          - a weighted kernel density without the central peak, where the
            kernel density is weighted by the distance from the zero mean
            value, using a 8th order power function between 0 and 100.

        In both cases, the weight gives more importance to samples far from
        the zero mean value. In the second case, the central peak is ignored.

    4.  The score, ranging from 0 to 100, is the sum of the squared weighted
        kernel density without the central peak, normalized by the sum of
        the squared full weighted kernel density. The score is 0 if there is no
        additional peak beyond the central peak.

    Parameters
    ----------
    trace : :class:`~obspy.core.trace.Trace`
        Trace to check for clipping.
    remove_baseline : bool, optional
        If set, remove the trace baseline before computing the score.
    debug : bool, optional
        If set, plot the trace, the samples histogram, the kernel density
        (unweighted and weighted), the kernel baseline model, and the misfit.
        Default is False.

    Returns
    -------
    clipping_score : float
        Clipping score, in percentage.

    Note
    ----
    Distorted traces (e.g., signals with strong baselines) can also get a
    high clipping score. To avoid this, the trace baseline can be removed.

    A debug mode is available to plot the trace, the samples histogram, the
    kernel density (unweighted and weighted), the kernel baseline model,
    and the misfit.
    """
    trace, trace_baseline, max_data = _preprocess(
        trace, remove_linear_trend=True, remove_baseline=remove_baseline)
    density, density_points = _get_kernel_density(
        trace, -max_data, max_data, 101)
    # Distance weight for full weighted kernel density,
    # 8th order, between 1 and 100
    dist_weight_full = _get_distance_weight(
        density_points, order=8, min_weight=1, max_weight=100)
    # Compute full weighted kernel density, including the central part
    density_weight_full = density*dist_weight_full
    # Distance weight for weighted kernel density without central peak,
    # 8th order, between 0 and 100
    dist_weight_no_central = _get_distance_weight(
        density_points, order=8, min_weight=0, max_weight=100)
    # Compute weighted kernel density, excluding the central part
    density_weight_no_central = density*dist_weight_no_central
    # Compute the clipping score
    clipping_score = 100 *\
        np.sum(density_weight_no_central**2)/np.sum(density_weight_full**2)
    if debug:
        _plot_clipping_analysis(
            trace, max_data, density_points, density, density_weight_full,
            trace_baseline=trace_baseline,
            density_weight_no_central=density_weight_no_central,
            clipping_score=clipping_score)
    return clipping_score


def _get_plotting_axes():
    """Get matplotlib axes for plotting"""
    # Force loading of a matplotlib GUI backend
    import matplotlib
    mpl_backends = 'macosx', 'qt5agg', 'qt4agg', 'gtk3agg', 'tkagg', 'wxagg'
    for backend in mpl_backends:
        try:
            matplotlib.use(backend, force=True)
            from matplotlib import pyplot  #noqa
            break
        except Exception:
            continue
    import matplotlib.pyplot as plt
    from matplotlib.ticker import ScalarFormatter

    class ScalarFormatterForceFormat(ScalarFormatter):
        def _set_format(self, *args):
            self.format = '%1.1f'

    fig, (ax_trace, ax_density) = plt.subplots(
        1, 2, figsize=(15, 5), sharey=True)
    yfmt = ScalarFormatterForceFormat()
    yfmt.set_powerlimits((0, 0))
    ax_trace.yaxis.set_major_formatter(yfmt)
    ax_trace.grid(True)
    ax_trace.set_xlabel('Time (s)')
    ax_trace.set_ylabel('Amplitude')
    ax_density.grid(True, axis='y')
    ax_density.set_xlabel('Density')
    return plt, ax_trace, ax_density


def _plot_clipping_analysis(
        trace, max_data, density_points, density, density_weight,
        trace_baseline=None, density_weight_no_central=None,
        clipping_score=None,
        peaks=None, num_edge_bins=None, num_kde_bins=None,
        trace_clipped=False):
    """
    Plot trace, samples histogram and kernel densities
    (unweighted and weighted)
    """
    plt, ax_trace, ax_density = _get_plotting_axes()
    # trace plot
    ax_trace.plot(
        trace.times(), trace.data, zorder=10, label='baseline removed')
    if trace_baseline is not None:
        ax_trace.plot(
            trace.times(), trace.data + trace_baseline, zorder=5,
            label='original trace')
        ax_trace.plot(trace.times(), trace_baseline, zorder=20, label='baseline')
        ax_trace.legend()
    ax_trace.set_ylim(-max_data, max_data)
    ax_trace.set_title(trace.id)
    # density plots
    counts, bins, bin_width = _get_histogram(trace.data)
    ax_density.hist(
        bins[:-1] + bin_width/2., bins=len(counts), weights=counts,
        orientation='horizontal', zorder=10)
    ax_density.plot(density, density_points, label='kernel density', zorder=20)
    ax_density.plot(
        density_weight, density_points,
        label='weighted kernel density', zorder=30)
    if density_weight_no_central is not None:
        ax_density.fill_betweenx(
            density_points, density_weight_no_central,
            alpha=0.5, color='gray',
            label='weighted kernel density\nno central peak', zorder=10)
    if peaks is not None:
        ax_density.scatter(
            density_weight[peaks], density_points[peaks],
            s=100, marker='x', color='red')
    ax_density.legend()
    if num_edge_bins is not None and num_kde_bins is not None:
        ax_density.axhspan(
            -max_data, density_points[num_edge_bins],
            alpha=0.5, color='yellow')
        ax_density.axhspan(
            density_points[num_kde_bins-1-num_edge_bins], max_data,
            alpha=0.5, color='yellow')
    if clipping_score is not None:
        ax_density.set_title(f'Clipping score: {clipping_score:.2f}%')
    if trace_clipped:
        ax_density.set_title('Clipped!')
    plt.tight_layout()
    plt.show()


def parse_arguments():
    """Parse command line arguments"""
    import sys
    import argparse
    parser = argparse.ArgumentParser(
        description='Check trace for clipping using kernel density '
        'estimation. Two methods are available: (1) check if trace is '
        'clipped, based on counting peaks in kernel density; (2) compute '
        'a trace clipping score based on kernel density estimation.')
    subparser = parser.add_subparsers(dest='command')
    common_parser = argparse.ArgumentParser(add_help=False)
    common_parser.add_argument(
        'infile', nargs='+', help='Input file(s) in any format supported by '
        'ObsPy (e.g., miniSEED, SAC, GSE2, SU, etc.)')
    common_parser.add_argument(
        '--traceid', '-t', help='Only process this trace ID. Trace ID can '
        'contain wildcards (e.g., "IU.ANMO.00.BHZ", "IU.ANMO.00.*", '
        '"IU.*.00.BHZ", etc.)')
    common_parser.add_argument(
        '--remove_baseline', '-r', action='store_true',
        help='Remove trace baseline before processing',
        default=False)
    common_parser.add_argument(
        '--debug', '-d', action='store_true',
        help='Plot trace, samples histogram, kernel density, and clipping '
        'parameters', default=False)
    sp_is_clipped = subparser.add_parser(
        'is_clipped', help='Check if trace is clipped, based on counting '
        'peaks in kernel density (method 1)', parents=[common_parser]
    )
    sp_is_clipped.add_argument(
        '-s', '--sensitivity', type=int, default=3,
        help='Sensitivity level, from 1 (least sensitive) '
        'to 5 (most sensitive). Default is %(default)s.')
    sp_is_clipped.add_argument(
        '-p', '--clipping_percentile', type=float, default=10,
        help='Percentile of trace amplitude range (expressed as percentage) '
        'to check for clipping. Default is %(default)s%%.')
    subparser.add_parser(
        'clipping_score', help='Get trace clipping score from kernel density '
        '(method 2)', parents=[common_parser]
    )
    args = parser.parse_args()
    if args.command is None:
        parser.print_usage(sys.stderr)
        sys.stderr.write(
            'Error: at least one positional argument is required.\n')
        sys.exit(2)
    return args


def run():
    from obspy import read, Stream
    args = parse_arguments()
    st = Stream()
    for file in args.infile:
        try:
            st += read(file)
        except Exception as msg:
            print(f'Error reading file {file}: {msg}')
    if args.traceid is not None:
        st = st.select(id=args.traceid)
    if not st:
        print('No traces found')
        return
    for tr in st:
        if args.command == 'is_clipped':
            trace_clipped, properties = is_clipped(
                tr, args.sensitivity, args.clipping_percentile, args.debug)
            print(
                f'{tr.id} - clipped: {trace_clipped}, '
                f'total peaks: {properties["npeaks"]}, '
                f'clipped peaks: {properties["npeaks_clipped"]}'
            )
        elif args.command == 'clipping_score':
            clipping_score = get_clipping_score(
                tr, args.remove_baseline, args.debug)
            print(f'{tr.id} - clipping score: {clipping_score:.2f}%')


def main():
    import sys
    try:
        run()
    except Exception as msg:
        sys.exit(msg)
    except KeyboardInterrupt:
        sys.exit()


if __name__ == '__main__':
    main()
