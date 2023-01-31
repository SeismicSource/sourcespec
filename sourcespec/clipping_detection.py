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


def is_clipped(trace, sensitivity, debug=False):
    """
    Check if a trace is clipped, based on kernel density estimation.

    Kernel density estimation is used to find the peaks of the histogram of
    the trace data points. The peaks are then weighted by their distance from
    the trace average (which should be the most common value).
    The peaks with the highest weight are then checked for prominence,
    which is a measure of how much higher the peak is than the surrounding
    data. The prominence threshold is determined by the sensitivity parameter.
    If more than one peak is found, the trace is considered clipped or
    distorted.

    Parameters
    ----------
    trace : obspy.core.trace.Trace
        Trace to check.
    sensitivity : int
        Sensitivity level, from 1 (least sensitive) to 5 (most sensitive).
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
    trace = trace.copy().detrend('demean')
    npts = len(trace.data)
    # Compute data histogram with a number of bins equal to 0.5% of data points
    nbins = int(npts*0.005)
    counts, bins = np.histogram(trace.data, bins=nbins)
    counts = counts/np.max(counts)
    # Compute gaussian kernel density
    kde = gaussian_kde(trace.data, bw_method=0.2)
    max_data = np.max(np.abs(trace.data))*1.2
    density_points = np.linspace(-max_data, max_data, 100)
    density = kde.pdf(density_points)
    maxdensity = np.max(density)
    density /= maxdensity
    # Distance weight, parabolic, between 1 and 5
    dist_weight = np.abs(density_points)**2
    dist_weight *= 4/dist_weight.max()
    dist_weight += 1
    density_weight = density*dist_weight
    # find peaks with minimum prominence based on clipping sensitivity
    min_prominence = [0.1, 0.05, 0.03, 0.02, 0.01]
    peaks, _ = find_peaks(
        density_weight,
        prominence=min_prominence[sensitivity-1]
    )
    if debug:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 2, figsize=(15, 5), sharey=True)
        fig.suptitle(trace.id)
        ax[0].plot(trace.times(), trace.data)
        ax[0].set_ylim(-max_data, max_data)
        ax[0].set_xlabel('Time (s)')
        ax[0].set_ylabel('Amplitude')
        ax[1].hist(
            bins[:-1], bins=len(counts), weights=counts,
            orientation='horizontal')
        ax[1].plot(density, density_points, label='kernel density')
        ax[1].plot(
            density_weight, density_points, label='weighted\nkernel density')
        ax[1].scatter(
            density_weight[peaks], density_points[peaks],
            s=100, marker='x', color='red')
        ax[1].set_xlabel('Density')
        ax[1].legend()
        plt.show()
    # If more than one peak, then the signal is probably clipped or distorted
    if len(peaks) > 1:
        return True
    else:
        return False


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
        '-d', '--debug', action='store_true',
        help='If set, plot trace, samples histogram and kernel density')
    args = parser.parse_args()
    st = read(args.infile)
    for tr in st:
        print(tr.id, is_clipped(tr, args.sensitivity, args.debug))


if __name__ == '__main__':
    import sys
    try:
        main()
    except Exception as msg:
        sys.exit(msg)
    except KeyboardInterrupt:
        sys.exit()
