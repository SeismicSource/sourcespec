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
from scipy.optimize import curve_fit


def _exp_func(x, a, b):
    """Exponential function for fitting kernel density."""
    return a*np.exp(-b*np.abs(x))


def _get_baseline(signal):
    """Get the signal baseline using a Savitzky-Golay filter."""
    from scipy.signal import savgol_filter
    npts = len(signal)
    wlen = npts // 5
    if wlen % 2 == 0:
        wlen += 1
    return savgol_filter(signal, wlen, 3)


def get_clipping_score(trace, remove_baseline=False, debug=False):
    """
    Compute a trace clipping score based on kernel density estimation.

    Trace is demeaned and detrended before computing the score.
    Optionally, the trace baseline can be removed before computing the score.

    The score, ranging from 0 to 100, is calculated as the misfit between the
    kernel density of trace amplitude values (weighted by distance from the
    zero mean value) and an exponential fit to the kernel density.

    Distance weighting (8th order power function between 1 and 100) is used
    to give more weight to samples far from the zero mean value.

    Distorted traces (e.g., signals with strong baselines) will also get a
    high clipping score. To avoid this, the trace baseline can be removed.

    A debug mode is available to plot the trace, the samples histogram and the
    kernel density (unweighted and weighted) and the exponential fit.

    Parameters
    ----------
    trace : :class:`~obspy.core.trace.Trace`
        Trace to check for clipping.
    remove_baseline : bool, optional
        If set, remove the trace baseline before computing the score.
    debug : bool, optional
        If set, plot trace, samples histogram, kernel density and the
        exponential fit. Default is False.

    Returns
    -------
    clipping_score : float
        Clipping score, in percentage.
    """
    trace = trace.copy().detrend('linear')
    trace.detrend('demean')
    max_data = np.max(np.abs(trace.data)) * 1.05  # 5% margin
    baseline = None
    if remove_baseline:
        baseline = _get_baseline(trace.data)
        trace.data -= baseline
    # Compute data histogram with a number of bins equal to 0.5% of data points
    # or 31, whichever is larger
    nbins = max(31, int(len(trace.data)*0.005))
    if nbins % 2 == 0:
        nbins += 1
    counts, bins = np.histogram(trace.data, bins=nbins)
    counts = counts/np.max(counts)
    bin_width = bins[1] - bins[0]
    # Compute gaussian kernel density
    kde = gaussian_kde(trace.data, bw_method=0.1)
    n_density_points = 101
    density_points = np.linspace(-max_data, max_data, n_density_points)
    density = kde.pdf(density_points)
    maxdensity = np.max(density)
    density /= maxdensity
    # Fit an exponential function to the kernel density
    popt, _ = curve_fit(
        _exp_func, density_points/max_data, density, p0=[1, 1])
    density_fit = _exp_func(density_points/max_data, *popt)
    # Distance weight, 8th order power function, between 1 and 100
    dist_weight = np.abs(density_points)**8
    dist_weight *= 99/dist_weight.max()
    dist_weight += 1
    density_weight = density*dist_weight
    # Compute misfit between weighted density and fit
    misfit = density_weight - density_fit
    # Compute clipping score
    clipping_score = 100 * np.sum(misfit**2)/np.sum(density_weight**2)
    if debug:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 2, figsize=(15, 5), sharey=True)
        ax[0].plot(
            trace.times(), trace.data, zorder=10, label='baseline removed')
        if baseline is not None:
            ax[0].plot(
                trace.times(), trace.data + baseline, zorder=5,
                label='original trace')
            ax[0].plot(trace.times(), baseline, zorder=20, label='baseline')
            ax[0].legend()
        ax[0].set_ylim(-max_data, max_data)
        ax[0].ticklabel_format(axis='y', style='sci', scilimits=(-2, 2))
        ax[0].grid(True)
        ax[0].set_title(trace.id)
        ax[0].set_xlabel('Time (s)')
        ax[0].set_ylabel('Amplitude')
        ax[1].hist(
            bins[:-1] + bin_width/2., bins=len(counts), weights=counts,
            orientation='horizontal', zorder=10)
        ax[1].plot(density, density_points, label='kernel density', zorder=20)
        ax[1].plot(
            density_weight, density_points, label='weighted\nkernel density',
            zorder=30)
        ax[1].plot(density_fit, density_points, label='exp. fit', zorder=40)
        ax[1].grid(True, axis='y')
        ax[1].set_title(f'Clipping score: {clipping_score:.2f}%')
        ax[1].set_xlabel('Density')
        ax[1].legend()
        plt.show()
    return clipping_score


def run():
    import argparse
    from obspy import read
    parser = argparse.ArgumentParser(
        description='Compute trace clipping score')
    parser.add_argument(
        'infile', type=str, help='input file name')
    parser.add_argument(
        '--remove_baseline', '-r', action='store_true',
        help='remove trace baseline before computing the score',
        default=False)
    parser.add_argument(
        '--debug', '-d', action='store_true',
        help='plot trace, samples histogram, kernel density '
        'and the exponential fit', default=False)
    args = parser.parse_args()
    st = read(args.infile)
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
