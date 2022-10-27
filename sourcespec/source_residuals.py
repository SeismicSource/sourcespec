# -*- coding: utf8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Compute mean station residuals from source_spec output.

:copyright:
    2013-2014 Claudio Satriano <satriano@ipgp.fr>,
              Agnes Chounet <chounet@ipgp.fr>

    2015-2022 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
import sys
import os
from collections import defaultdict
import pickle
from argparse import ArgumentParser
from obspy.core import Stream
from sourcespec.ssp_util import moment_to_mag, mag_to_moment
from sourcespec.spectrum import Spectrum
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('Agg')  # NOQA


def main():
    parser = ArgumentParser(
        description='Compute mean station residuals from source_spec output.')
    parser.add_argument(
        '-m', '--min_spectra', dest='min_spectra', type=int, action='store',
        default='20', help='minimum number of spectra to '
        'compute residuals (default=20)', metavar='NUMBER')
    parser.add_argument(
        '-p', '--plot', dest='plot', action='store_true',
        default=False, help='save residuals plots to file')
    parser.add_argument(
        'residual_files_dir',
        help='directory containing source_spec residual files '
             'in pickle format. Residual files can be in subdirectories '
             '(e.g., a subdirectory for each event).')
    args = parser.parse_args()

    min_spectra = int(args.min_spectra)
    resfiles_dir = args.residual_files_dir
    outdir = 'sspec_residuals'

    if not os.path.exists(resfiles_dir):
        print('Error: directory "{}" does not exist.'.format(resfiles_dir))
        sys.exit(1)

    resfiles = []
    for root, dirs, files in os.walk(resfiles_dir):
        for file in files:
            if file.endswith('residuals.pickle'):
                resfiles.append(os.path.join(root, file))
    if not resfiles:
        print('No residual file found in directory: {}'.format(resfiles_dir))
        sys.exit()

    residual_dict = defaultdict(Stream)
    for resfile in resfiles:
        print('Found residual file: {}'.format(resfile))
        with open(resfile, 'rb') as fp:
            residual_st = pickle.load(fp)
        for spec in residual_st:
            residual_dict[spec.id].append(spec)

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    residual_mean = Stream()
    for stat_id in sorted(residual_dict.keys()):
        if len(residual_dict[stat_id]) < min_spectra:
            continue
        print('Processing station: {}'.format(stat_id))

        res = residual_dict[stat_id]

        freqs_min = [spec.get_freq().min() for spec in res]
        freqs_max = [spec.get_freq().max() for spec in res]
        freq_min = min(freqs_min)
        freq_max = max(freqs_max)

        spec_mean = Spectrum()
        spec_mean.id = stat_id
        spec_mean.stats.begin = freq_min
        spec_mean.stats.delta = res[0].stats.delta
        spec_mean.data_mag = None
        for spec in res:
            spec_slice = spec.slice(freq_min, freq_max, pad=True,
                                    fill_value=mag_to_moment(0))
            spec_slice.data_mag = moment_to_mag(spec_slice.data)
            norm = (spec_slice.data_mag != 0).astype(int)
            if spec_mean.data_mag is None:
                spec_mean.data_mag = spec_slice.data_mag
                norm_mean = norm
            else:
                spec_mean.data_mag += spec_slice.data_mag
                norm_mean += norm
        spec_mean.data_mag /= norm_mean
        spec_mean.data = mag_to_moment(spec_mean.data_mag)

        residual_mean.append(spec_mean)

        # plot traces
        if args.plot:
            figurefile = os.path.join(outdir, '{}-res.png'.format(stat_id))
            fig = plt.figure(dpi=160)
            for spec in res:
                plt.semilogx(spec.get_freq(), spec.data_mag, 'b-')
            plt.semilogx(spec_mean.get_freq(), spec_mean.data_mag, 'r-')
            plt.xlabel('frequency (Hz)')
            plt.ylabel('residual amplitude (obs - synth) in magnitude units')
            plt.title(
                'residuals: {} – {} records'.format(stat_id, len(res))
            )
            fig.savefig(figurefile, bbox_inches='tight')
            plt.close()
            print('Residual plot saved to: {}'.format(figurefile))

    # writes the mean residuals (the stations corrections)
    res_mean_file = 'residual_mean.pickle'
    res_mean_file = os.path.join(outdir, res_mean_file)
    with open(res_mean_file, 'wb') as fp:
        pickle.dump(residual_mean, fp)
    print('Mean station residuals saved to: {}'.format(res_mean_file))
