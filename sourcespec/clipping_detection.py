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


def _get_baseline(signal):
    """Get the signal baseline using a Savitzky-Golay filter."""
    from scipy.signal import savgol_filter
    npts = len(signal)
    wlen = npts // 10
    if wlen % 2 == 0:
        wlen += 1
    return savgol_filter(signal, wlen, 3)


def _plot_clipping_analysis(
        trace, trace_baseline, max_data,
        counts, bins, bin_width,
        density_points, density,
        density_weight_full, density_weight_no_central,
        clipping_score):
    """
    Plot trace, samples histogram and kernel densities
    (unweighted and weighted)
    """
    # Force loading of a matplotlib GUI backend
    import matplotlib
    mpl_backends = 'macosx', 'qt5agg', 'qt4agg', 'gtk3agg', 'tkagg', 'wxagg'
    for backend in mpl_backends:
        try:
            matplotlib.use(backend, force=True)
            from matplotlib import pyplot as plt
            break
        except Exception:
            continue
    import matplotlib.pyplot as plt
    from matplotlib.ticker import ScalarFormatter

    class ScalarFormatterForceFormat(ScalarFormatter):
        def _set_format(self):
            self.format = '%1.1f'

    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(15, 5), sharey=True)
    ax0.plot(
        trace.times(), trace.data, zorder=10, label='baseline removed')
    if trace_baseline is not None:
        ax0.plot(
            trace.times(), trace.data + trace_baseline, zorder=5,
            label='original trace')
        ax0.plot(trace.times(), trace_baseline, zorder=20, label='baseline')
        ax0.legend()
    ax0.set_ylim(-max_data, max_data)
    yfmt = ScalarFormatterForceFormat()
    yfmt.set_powerlimits((0,0))
    ax0.yaxis.set_major_formatter(yfmt)
    ax0.grid(True)
    ax0.set_title(trace.id)
    ax0.set_xlabel('Time (s)')
    ax0.set_ylabel('Amplitude')
    ax1.hist(
        bins[:-1] + bin_width/2., bins=len(counts), weights=counts,
        orientation='horizontal', zorder=10)
    ax1.plot(density, density_points, label='kernel density', zorder=20)
    ax1.plot(
        density_weight_full, density_points,
        label='weighted kernel density', zorder=30)
    ax1.fill_betweenx(
        density_points, density_weight_no_central,
        alpha=0.5, color='gray',
        label='weighted kernel density\nno central peak', zorder=10)
    ax1.grid(True, axis='y')
    ax1.set_xlabel('Density')
    ax1.legend()
    ax1.set_title(f'Clipping score: {clipping_score:.2f}%')
    plt.tight_layout()
    plt.show()


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

    4.  The score, ranging from 0 to 100, is the ratio between the sum of the
        squares of the values of the full weighted kernel density and the
        sum of the squares of the values of the weighted kernel density without
        the central peak. The score is 0 if there is no additional peak beyond
        the central peak.

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

    Notes
    -----
    Distorted traces (e.g., signals with strong baselines) can also get a
    high clipping score. To avoid this, the trace baseline can be removed.

    A debug mode is available to plot the trace, the samples histogram, the
    kernel density (unweighted and weighted), the kernel baseline model,
    and the misfit.
    """
    trace = trace.copy().detrend('linear')
    trace.detrend('demean')
    max_data = np.max(np.abs(trace.data)) * 1.10  # 10% margin
    # Remove trace baseline if requested
    trace_baseline = None
    if remove_baseline:
        trace_baseline = _get_baseline(trace.data)
        trace.data -= trace_baseline
    # Compute data histogram with a number of bins equal to 0.5% of data points
    # or 31, whichever is larger
    # Note: histogram is only used for plotting, the score is computed from
    # the kernel density
    nbins = max(31, int(len(trace.data)*0.005))
    if nbins % 2 == 0:
        nbins += 1
    counts, bins = np.histogram(trace.data, bins=nbins)
    counts = counts/np.max(counts)
    bin_width = bins[1] - bins[0]
    # Compute gaussian kernel density using default bandwidth (Scott’s Rule)
    kde = gaussian_kde(trace.data)
    n_density_points = 101
    density_points = np.linspace(-max_data, max_data, n_density_points)
    density = kde.pdf(density_points)
    density /= np.max(density)
    # Base distance weight for kernel density, 8th order power function
    # between 0 and 1
    dist_weight = np.abs(density_points)**8
    dist_weight /= dist_weight.max()
    # Distance weight between 1 and 100
    dist_weight_full = 1 + 99*dist_weight
    # Compute full weighted kernel density, including the central part
    density_weight_full = density*dist_weight_full
    # Distance weight between 0 and 100
    dist_weight_no_central = 100*dist_weight
    # Compute weighted kernel density, excluding the central part
    density_weight_no_central = density*dist_weight_no_central
    # Compute the clipping score
    clipping_score = 100 *\
        np.sum(density_weight_no_central**2)/np.sum(density_weight_full**2)
    if debug:
        _plot_clipping_analysis(
            trace, trace_baseline, max_data,
            counts, bins, bin_width,
            density_points, density,
            density_weight_full, density_weight_no_central,
            clipping_score)
    return clipping_score


def run():
    import argparse
    from obspy import read, Stream
    parser = argparse.ArgumentParser(
        description='Compute trace clipping score')
    parser.add_argument(
        'infile', nargs='+', help='input file(s) in any format supported by '
        'ObsPy (e.g., miniSEED, SAC, GSE2, SU, etc.)')
    parser.add_argument(
        '--remove_baseline', '-r', action='store_true',
        help='remove trace baseline before computing the score',
        default=False)
    parser.add_argument(
        '--debug', '-d', action='store_true',
        help='plot trace, samples histogram, kernel density '
        'kernel baseline model and misfit', default=False)
    args = parser.parse_args()
    st = Stream()
    for file in args.infile:
        try:
            st += read(file)
        except Exception as msg:
            print(f'Error reading file {file}: {msg}')
    for tr in st:
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
