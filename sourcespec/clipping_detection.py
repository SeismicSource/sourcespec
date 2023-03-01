# -*- coding: utf-8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
A trace clipping detector based on kernel density estimation.

:copyright:
    2023 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
import numpy as np
from scipy.stats import gaussian_kde
from scipy.signal import find_peaks


def _plot_clipping_analysis(
        trace, max_data, density, density_points, density_weight,
        peaks, num_edge_bins, num_kde_bins):
    import matplotlib.pyplot as plt
    from matplotlib.ticker import ScalarFormatter

    class ScalarFormatterForceFormat(ScalarFormatter):
        def _set_format(self):
            self.format = '%1.1f'

    fig, ax = plt.subplots(1, 2, figsize=(15, 5), sharey=True)
    ax[0].plot(trace.times(), trace.data)
    ax[0].set_ylim(-max_data, max_data)
    yfmt = ScalarFormatterForceFormat()
    yfmt.set_powerlimits((0, 0))
    ax[0].yaxis.set_major_formatter(yfmt)
    ax[0].grid()
    ax[0].set_title(trace.id)
    ax[0].set_xlabel('Time (s)')
    ax[0].set_ylabel('Amplitude')

    npts = len(trace.data)
    # Compute data histogram with a number of bins equal to 0.5% of data points
    # or 11, whichever is greater
    nbins = max(11, int(npts*0.005))
    if nbins % 2 == 0:
        nbins += 1
    counts, bins = np.histogram(trace.data, bins=nbins)
    counts = counts/np.max(counts)
    bin_width = bins[1] - bins[0]
    ax[1].hist(
        bins[:-1] + bin_width/2., bins=len(counts), weights=counts,
        orientation='horizontal')
    ax[1].plot(density, density_points, label='kernel density')
    ax[1].plot(
        density_weight, density_points, label='weighted\nkernel density')
    ax[1].scatter(
        density_weight[peaks], density_points[peaks],
        s=100, marker='x', color='red')
    ax[1].set_xlabel('Density')
    ax[1].legend()
    ax[1].axhspan(
        -max_data, density_points[num_edge_bins], alpha=0.5, color='yellow')
    ax[1].axhspan(
        density_points[num_kde_bins-1-num_edge_bins], max_data,
        alpha=0.5, color='yellow')
    ax[1].grid(axis='y')
    plt.tight_layout()
    plt.show()



def is_clipped(trace, sensitivity=3, lower_clip_bound=90, debug=False):
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
    lower_clip_bound : int
        Lower bound of amplitude range (expressed as percentage) to consider
        as potentially clipped
        (default: 90)
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
    lower_clip_bound = (min(100, max(0, lower_clip_bound)))
    num_kde_bins = 101
    num_edge_bins = int(np.ceil(
        (num_kde_bins/2.) * (100 - lower_clip_bound) / 100.))
    trace = trace.copy().detrend('demean')
    # Compute gaussian kernel density
    # Note: current value of bw_method is optimized for num_kde_bins = 101
    kde = gaussian_kde(trace.data, bw_method=0.1)
    min_data = np.min(trace.data)
    max_data = np.max(trace.data)
    density_points = np.linspace(min_data, max_data, num_kde_bins)
    density = kde.pdf(density_points)
    maxdensity = np.max(density)
    density /= maxdensity
    # Distance weight, parabolic, between 1 and 5
    dist_weight = np.abs(density_points)**2
    dist_weight *= 4/dist_weight.max()
    dist_weight += 1
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

    if debug:
        print('%2d peaks found:' % len(peaks))
        for peak, prominence in zip(peaks, props['prominences']):
            print('  idx %d: prominence=%G' % (peak, prominence))
        _plot_clipping_analysis(
            trace, max_data, density, density_points, density_weight,
            peaks, num_edge_bins, num_kde_bins)
    # If there is a peak in the edge bins,
    # then the signal is probably clipped or distorted
    return any(
        peak < num_edge_bins or peak > (num_kde_bins - 1 - num_edge_bins)
        for peak in peaks
    )


def main():
    import argparse
    from obspy import read
    parser = argparse.ArgumentParser(
        description='Check if a trace is clipped, '
        'based on kernel density estimation of trace samples.')
    parser.add_argument(
        'infile', type=str,
        help='Input file name in any ObsPy supported format')
    parser.add_argument(
        '-s', '--sensitivity', type=int, default=3,
        help='Sensitivity level, from 1 (least sensitive) '
        'to 5 (most sensitive)')
    parser.add_argument(
        '-l', '--lower_clip_bound', type=int, default=90,
        help='Lower bound of amplitude range (expressed as percentage) '
        'to consider as potentially clipped')
    parser.add_argument(
        '-d', '--debug', action='store_true',
        help='If set, plot trace, samples histogram and kernel density')
    args = parser.parse_args()
    st = read(args.infile)
    for tr in st:
        print(tr.id, is_clipped(
                tr, args.sensitivity, args.lower_clip_bound, args.debug))


if __name__ == '__main__':
    import sys
    try:
        main()
    except Exception as msg:
        sys.exit(msg)
    except KeyboardInterrupt:
        sys.exit()
