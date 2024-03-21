# -*- coding: utf8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Compute mean station residuals from source_spec output.

:copyright:
    2013-2014 Claudio Satriano <satriano@ipgp.fr>,
              Agnes Chounet <chounet@ipgp.fr>

    2015-2024 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
import sys
import os
from collections import defaultdict
import pickle
from argparse import ArgumentParser
import matplotlib
import matplotlib.pyplot as plt
from sourcespec.ssp_util import moment_to_mag, mag_to_moment
from sourcespec.spectrum import SpectrumStream
matplotlib.use('Agg')  # NOQA


def parse_args():
    """
    Parse command line arguments.
    """
    parser = ArgumentParser(
        description='Compute mean station residuals from source_spec output.')
    parser.add_argument(
        '-m', '--min_spectra', dest='min_spectra', type=int, action='store',
        default='20', help='minimum number of spectra to '
        'compute residuals (default=20)', metavar='NUMBER')
    parser.add_argument(
        '-o', '--outdir', dest='outdir', action='store',
        default='sspec_residuals',
        help='output directory (default="sspec_residuals")')
    parser.add_argument(
        '-p', '--plot', dest='plot', action='store_true',
        default=False, help='save residuals plots to file')
    parser.add_argument(
        'residual_files_dir',
        help='directory containing source_spec residual files '
             'in pickle format. Residual files can be in subdirectories '
             '(e.g., a subdirectory for each event).')
    return parser.parse_args()


def read_residuals(resfiles_dir):
    """
    Read residuals from pickle files in resfiles_dir.

    Parameters
    ----------
    resfiles_dir : str
        Directory containing source_spec residual files in pickle format.
        Residual files can be in subdirectories (e.g., a subdirectory for
        each event).

    Returns
    -------
    residual_dict : dict
        Dictionary containing residuals for each station.
    """
    if not os.path.exists(resfiles_dir):
        sys.exit(f'Error: directory "{resfiles_dir}" does not exist.')
    resfiles = []
    for root, _dirs, files in os.walk(resfiles_dir):
        resfiles.extend(
            os.path.join(root, file)
            for file in files
            if file.endswith('residuals.pickle')
        )
    if not resfiles:
        sys.exit(f'No residual file found in directory: {resfiles_dir}')
    residual_dict = defaultdict(SpectrumStream)
    for resfile in resfiles:
        print(f'Found residual file: {resfile}')
        with open(resfile, 'rb') as fp:
            residual_st = pickle.load(fp)
        for spec in residual_st:
            residual_dict[spec.id].append(spec)
    return residual_dict


def compute_mean_residuals(residual_dict, min_spectra=20):
    """
    Compute mean residuals for each station.

    Parameters
    ----------
    residual_dict : dict
        Dictionary containing residuals for each station.
    min_spectra : int
        Minimum number of spectra to compute residuals (default=20).

    Returns
    -------
    residual_mean : SpectrumStream
        Stream containing mean residuals for each station.
    """
    residual_mean = SpectrumStream()
    for stat_id in sorted(residual_dict.keys()):
        if len(residual_dict[stat_id]) < min_spectra:
            continue
        print(f'Processing station: {stat_id}')

        res = residual_dict[stat_id]

        freqs_min = [spec.freq.min() for spec in res]
        freqs_max = [spec.freq.max() for spec in res]
        freq_min = min(freqs_min)
        freq_max = max(freqs_max)

        spec_mean = None
        for spec in res:
            spec_slice = spec.slice(freq_min, freq_max, pad=True,
                                    fill_value=mag_to_moment(0))
            spec_slice.data_mag = moment_to_mag(spec_slice.data)
            norm = (spec_slice.data_mag != 0).astype(int)
            if spec_mean is None:
                spec_mean = spec_slice
                norm_mean = norm
            else:
                try:
                    spec_mean.data_mag += spec_slice.data_mag
                except ValueError:
                    continue
                norm_mean += norm
        spec_mean.data_mag /= norm_mean
        spec_mean.data = mag_to_moment(spec_mean.data_mag)

        residual_mean.append(spec_mean)
    return residual_mean


def plot_residuals(residual_dict, residual_mean, outdir):
    """
    Plot residuals.

    Parameters
    ----------
    residual_dict : dict
        Dictionary containing residuals for each station.
    residual_mean : SpectrumStream
        SpectrumStream containing mean residuals for each station.
    outdir : str
        Output directory.
    """
    for spec_mean in residual_mean:
        stat_id = spec_mean.id
        res = residual_dict[stat_id]
        figurefile = os.path.join(outdir, f'{stat_id}-res.png')
        fig = plt.figure(dpi=160)
        for spec in res:
            plt.semilogx(spec.freq, spec.data_mag, 'b-')
        plt.semilogx(spec_mean.freq, spec_mean.data_mag, 'r-')
        plt.xlabel('frequency (Hz)')
        plt.ylabel('residual amplitude (obs - synth) in magnitude units')
        plt.title(f'residuals: {stat_id} – {len(res)} records')
        fig.savefig(figurefile, bbox_inches='tight')
        plt.close()
        print(f'Residual plot saved to: {figurefile}')


def main():
    """Main function."""
    args = parse_args()

    residual_dict = read_residuals(args.residual_files_dir)

    outdir = args.outdir
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    min_spectra = int(args.min_spectra)
    residual_mean = compute_mean_residuals(residual_dict, min_spectra)

    if args.plot:
        plot_residuals(residual_dict, residual_mean, outdir)

    # writes the mean residuals (the stations corrections)
    res_mean_file = 'residual_mean.pickle'
    res_mean_file = os.path.join(outdir, res_mean_file)
    with open(res_mean_file, 'wb') as fp:
        pickle.dump(residual_mean, fp)
    print(f'Mean station residuals saved to: {res_mean_file}')
