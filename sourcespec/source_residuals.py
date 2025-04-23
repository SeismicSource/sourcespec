# -*- coding: utf8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Compute mean station residuals from source_spec output.

:copyright:
    2013-2014 Claudio Satriano <satriano@ipgp.fr>,
              Agnes Chounet <chounet@ipgp.fr>

    2015-2025 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
import sys
import os
from collections import defaultdict
from argparse import ArgumentParser
import numpy as np
from scipy.interpolate import interp1d
import matplotlib
import matplotlib.pyplot as plt
from sourcespec._version import get_versions
from sourcespec.spectrum import read_spectra
from sourcespec.ssp_util import mag_to_moment
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
        'compute residuals (default=20)', metavar='NUMBER'
    )
    parser.add_argument(
        '-r', '--runid', dest='runid', action='store',
        default=None, metavar='RUNID',
        help='Select a specific run when multiple runs exist for the same '
             'event. If omitted, residuals will be computed using all '
             'available runs. You can provide a specific RUNID or use '
             '"latest" or "earliest" to select the most recent or the '
             'earliest run, based on alphabetical or numerical order.'
    )
    parser.add_argument(
        '-o', '--outdir', dest='outdir', action='store',
        default='sspec_residuals',
        help='output directory (default="sspec_residuals")'
    )
    parser.add_argument(
        '-e', '--exclude', dest='exclude_subdirs', action="append",
        default=[],
        help='subfolder to exclude (repeat for multiple subfolders)'
    )
    parser.add_argument(
        '-p', '--plot', dest='plot', action='store_true',
        default=False, help='save residuals plots to file'
    )
    parser.add_argument(
        'residual_files_dir',
        help='directory containing source_spec residual files '
             'in HDF5 format. Residual files can be in subdirectories '
             '(e.g., a subdirectory for each event).')
    return parser.parse_args()


def _filter_by_runid(residual_dict, runid='latest'):
    """
    Filter the residual dictionary by the specified runid criteria.

    Parameters
    ----------
    residual_dict : dict
        Dictionary containing residuals for each station.
    runid : str
        Criteria to select the runid. Can be 'latest', 'earliest',
        or a specific runid.

    Returns
    -------
    filt_residual_dict : dict
        Dictionary containing residuals filtered by the specified criteria.
    """
    if runid is None:
        return residual_dict
    if not isinstance(runid, str):
        raise ValueError("runid must be a string.")
    filt_residual_dict = defaultdict(SpectrumStream)
    if runid not in ['latest', 'earliest']:
        for spec_st in residual_dict.values():
            # Filter by specific runid
            filt_residual_dict[spec_st[0].id].extend(
                [spec for spec in spec_st if spec.stats.runid == runid]
            )
        return filt_residual_dict
    # Filter by latest or earliest runid
    for spec_st in residual_dict.values():
        evids_runids = [
            (spec.stats.event.event_id, spec.stats.runid)
            for spec in spec_st
        ]
        # Dictionary to store the selected runid for each evid
        selected_runids = {}
        for _evid, _runid in evids_runids:
            if runid == 'latest':
                selected_runids[_evid] = max(
                    selected_runids.get(_evid, '0'), _runid, key=int)
            elif runid == 'earliest':
                selected_runids[_evid] = min(
                    selected_runids.get(_evid, '999999'), _runid, key=int)
        # Convert back to list of tuples
        filtered_evids_runids = list(selected_runids.items())
        filt_residual_dict[spec_st[0].id].extend(
            [spec for spec in spec_st
                if (spec.stats.event.event_id, spec.stats.runid)
                in filtered_evids_runids]
        )
    return filt_residual_dict


def read_residuals(resfiles_dir, runid=None, exclude_subdirs=[]):
    """
    Read residuals from HDF5 files in resfiles_dir.

    Parameters
    ----------
    resfiles_dir : str
        Directory containing source_spec residual files in HDF5 format.
        Residual files can be in subdirectories (e.g., a subdirectory for
        each event).
    exclude_subdirs: list of str
        List of subdirectories (potentially containing residual files)
        that should not be parsed
        (default: [])

    Returns
    -------
    residual_dict : dict
        Dictionary containing residuals for each station.
    """
    if not os.path.exists(resfiles_dir):
        sys.exit(f'Error: directory "{resfiles_dir}" does not exist.')
    resfiles = []
    for root, _dirs, files in os.walk(resfiles_dir):
        _dirs[:] = [d for d in _dirs if d not in exclude_subdirs]
        resfiles.extend(
            os.path.join(root, file)
            for file in files
            if file.endswith('residuals.hdf5')
        )
    if not resfiles:
        sys.exit(f'No residual file found in directory: {resfiles_dir}')
    residual_dict = defaultdict(SpectrumStream)
    for resfile in sorted(resfiles):
        print(f'Found residual file: {resfile}')
        residual_st = read_spectra(resfile)
        for spec in residual_st:
            residual_dict[spec.id].append(spec)
    return _filter_by_runid(residual_dict, runid)


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
        delta_f_min = min(spec.stats.delta for spec in res)
        freq_min = min(freqs_min)
        freq_max = max(freqs_max)
        freq_array = np.arange(freq_min, freq_max + delta_f_min, delta_f_min)

        spec_mean = None
        for spec in res:
            spec_interp = spec.copy()
            # interpolate to the new frequency array
            spec_interp.freq = freq_array
            # spec_slice.data must exist, so we create it as zeros
            spec_interp.data = np.zeros_like(freq_array)
            # interpolate data_mag to the new frequency array
            f = interp1d(spec.freq, spec.data_mag, bounds_error=False)
            spec_interp.data_mag = f(freq_array)
            # norm is 1 where interpolated data_mag is not nan, 0 otherwise
            norm = (~np.isnan(spec_interp.data_mag)).astype(float)
            # Replace nan data_mag values with zeros for summation
            spec_interp.data_mag[norm == 0] = 0
            if spec_mean is None:
                spec_mean = spec_interp
                norm_mean = norm
            else:
                spec_mean.data_mag += spec_interp.data_mag
                norm_mean += norm
        # Make sure to avoid division by zero
        norm_mean[norm_mean == 0] = np.nan
        spec_mean.data_mag /= norm_mean
        spec_mean.data = mag_to_moment(spec_mean.data_mag)
        spec_mean.stats.software = 'SourceSpec'
        spec_mean.stats.software_version = get_versions()['version']
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
        figurefile = os.path.join(outdir, f'{stat_id}.res.png')
        fig, ax = plt.subplots(dpi=160)
        for spec in res:
            ax.semilogx(spec.freq, spec.data_mag, 'b-', alpha=0.5)
        ax.semilogx(spec_mean.freq, spec_mean.data_mag, 'r-', linewidth=2)
        ax.set_xlabel('frequency (Hz)')
        ax.set_ylabel('residual amplitude (obs - synth) in magnitude units')
        ax.set_title(f'residuals: {stat_id} – {len(res)} records')
        fig.savefig(figurefile, bbox_inches='tight')
        plt.close(fig)
        print(f'Residual plot saved to: {figurefile}')


def main():
    """Main function."""
    args = parse_args()

    exclude_subdirs = args.exclude_subdirs
    residual_dict = read_residuals(args.residual_files_dir, args.runid,
                                    exclude_subdirs=exclude_subdirs)

    outdir = args.outdir
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    min_spectra = int(args.min_spectra)
    residual_mean = compute_mean_residuals(residual_dict, min_spectra)

    if args.plot:
        plot_residuals(residual_dict, residual_mean, outdir)

    # write the mean residuals (the stations corrections)
    res_mean_file = 'residual_mean.hdf5'
    res_mean_file = os.path.join(outdir, res_mean_file)
    if not residual_mean:
        print(
            'Mean residuals not computed: not enough spectra.\n'
            'Minimum number of spectra ("--min_spectra") currently set to: '
            f'{min_spectra}'
        )
        sys.exit(0)
    residual_mean.write(res_mean_file, format='HDF5')
    print(f'Mean station residuals saved to: {res_mean_file}')



if __name__ == '__main__':
    main()
