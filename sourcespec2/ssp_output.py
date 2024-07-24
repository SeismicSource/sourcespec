# -*- coding: utf-8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Output functions for source_spec.

:copyright:
    2012 Claudio Satriano <satriano@ipgp.fr>

    2013-2014 Claudio Satriano <satriano@ipgp.fr>,
              Emanuela Matrullo <matrullo@geologie.ens.fr>

    2015-2025 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
import os
import shutil
import contextlib
import logging
from collections.abc import Mapping
from datetime import datetime
from tzlocal import get_localzone
import numpy as np
from .setup import config
from .ssp_qml_output import write_qml
from .ssp_sqlite_output import write_sqlite
from ._version import get_versions
logger = logging.getLogger(__name__.rsplit('.', maxsplit=1)[-1])
# reduce logging level for tzlocal
logging.getLogger('tzlocal').setLevel(logging.WARNING)


def _write_author_and_agency_to_parfile(parfile):
    author_str = empty_author_str = '\n*** Author:'
    if config.author_name is not None:
        author_str += f' {config.author_name}'
    if config.author_email is not None:
        if author_str != empty_author_str:
            author_str += f' <{config.author_email}>'
        else:
            author_str += f' {config.author_email}'
    if author_str != empty_author_str:
        parfile.write(author_str)
    agency_str = empty_agency_str = '\n*** Agency:'
    if config.agency_full_name is not None:
        agency_str += f' {config.agency_full_name}'
    if config.agency_short_name is not None:
        if agency_str != empty_agency_str:
            agency_str += f' ({config.agency_short_name})'
        else:
            agency_str += f' {config.agency_short_name}'
    if config.agency_url is not None:
        if agency_str != empty_agency_str:
            agency_str += ' -'
        agency_str += f' {config.agency_url}'
    if agency_str != empty_agency_str:
        parfile.write(agency_str)


def _value_error_str(value, error, fmt):
    if error[0] == error[1]:
        s = f'{fmt} +/- {fmt}'
        s = s.format(value, error[0])
    else:
        s = f'{fmt} /- {fmt} /+ {fmt}'
        s = s.format(value, error[0], error[1])
    return s


def _write_parfile(sspec_output):
    """
    Write station source parameters to file.

    Note: this format is deprecated and will not evolve anymore
    (e.g., including new parameters or new way of computing statistics).
    """
    if not os.path.exists(config.options.outdir):
        os.makedirs(config.options.outdir)
    evid = config.event.event_id
    parfilename = os.path.join(
        config.options.outdir, f'{evid}.ssp.out')
    with open(parfilename, 'w', encoding='utf-8') as parfile:
        parfile.write(
            '*** Note: this file is deprecated and might not contain all the '
            'output information. ***\n')
        parfile.write(
            '*** It will be removed in future versions of SourceSpec. ***\n')
        parfile.write(
            '*** Please look at the YAML file in this directory. ***\n')

        evid = config.event.event_id
        hypo = config.event.hypocenter
        evlo = hypo.longitude.value_in_deg
        evla = hypo.latitude.value_in_deg
        evdp = hypo.depth.value_in_km
        parfile.write(
            f'{evid} lon {evlo:8.3f} lat {evla:7.3f} depth {evdp:5.1f} km '
            f'orig_time {hypo.origin_time}\n\n')
        parfile.write('*** Station source parameters ***\n')
        parfile.write(
            '*** Note: outliers are prepended by a star (*) symbol ***\n')
        parkeys = (
            'Mw', 'fc', 't_star', 'Qo', 'Mo',
            'ssd', 'ra', 'hyp_dist', 'az', 'Er'
        )
        formats = {
            'Mo': '{:.3e} ',
            'Er': '{:.3e} ',
            'hyp_dist': '{:7.3f} ',
            'az': '{:7.3f} ',
            'Mw': '{:6.3f} ',
            'fc': '{:6.3f} ',
            'ssd': '{:.3e} ',
            'ra': '{:8.3f} ',
            't_star': '{:6.3f} ',
            'Qo': '{:7.1f} ',
            'Ml': '{:6.3f} '
        }
        formats_none = {
            'Mo': '{:>9} ',
            'Er': '{:>9} ',
            'hyp_dist': '{:>7} ',
            'az': '{:>7} ',
            'Mw': '{:>6} ',
            'fc': '{:>6} ',
            'ssd': '{:>9} ',
            'ra': '{:>8} ',
            't_star': '{:>6} ',
            'Qo': '{:>7} ',
            'Ml': '{:>6} '
        }
        stationpar = sspec_output.station_parameters
        for statId in sorted(stationpar.keys()):
            par = stationpar[statId]
            parfile.write(f'{statId:>15} {par.instrument_type:>6}\t')
            for key in parkeys:
                if key == 'az':
                    _pkey = key
                    val = par['azimuth']
                    outl = False
                elif key == 'hyp_dist':
                    _pkey = key
                    val = par['hypo_dist_in_km']
                    outl = False
                elif key == 'ra':
                    _pkey = 'radius'
                    val = par[_pkey].value
                    outl = par[_pkey].outlier
                else:
                    _pkey = key
                    val = par[_pkey].value
                    outl = par[_pkey].outlier
                space = ' *' if outl else '  '
                parfile.write(f'{space}{key} ')
                if val is not None and ~np.isnan(val):
                    parfile.write(formats[key].format(val))
                else:
                    parfile.write(formats_none[key].format('nan'))
            parfile.write('\n')
            parfile.write(f'{"--- errmin":>22}\t')
            for key in parkeys:
                _pkey = 'radius' if key == 'ra' else key
                if key in ['hyp_dist', 'az']:
                    outl = False
                    err = None
                else:
                    outl = par[_pkey].outlier
                    err = par[_pkey].lower_uncertainty
                    if err is None:
                        err = par[_pkey].uncertainty
                space = ' *' if outl else '  '
                parfile.write(f'{space}{key} ')
                if err is not None:
                    parfile.write(formats[key].format(err))
                else:
                    parfile.write(formats_none[key].format('nan'))
            parfile.write('\n')
            parfile.write(f'{"--- errmax":>22}\t')
            for key in parkeys:
                _pkey = 'radius' if key == 'ra' else key
                if key in ['hyp_dist', 'az']:
                    outl = False
                    err = None
                else:
                    outl = par[_pkey].outlier
                    err = par[_pkey].upper_uncertainty
                    if err is None:
                        err = par[_pkey].uncertainty
                space = ' *' if outl else '  '
                parfile.write(f'{space}{key} ')
                if err is not None:
                    parfile.write(formats[key].format(err))
                else:
                    parfile.write(formats_none[key].format('nan'))
            parfile.write('\n')

        means = sspec_output.mean_values()
        errors = sspec_output.mean_uncertainties()
        means_weight = sspec_output.weighted_mean_values()
        errors_weight = sspec_output.weighted_mean_uncertainties()

        parfile.write('\n*** Average source parameters ***\n')
        parfile.write(
            '*** Note: averages computed after removing outliers ****\n')

        Mw_mean = means['Mw']
        Mw_error = errors['Mw']
        s = _value_error_str(Mw_mean, Mw_error, '{:.2f}')
        parfile.write(f'Mw: {s}\n')
        Mw_mean_weight = means_weight['Mw']
        Mw_error_weight = errors_weight['Mw']
        s = _value_error_str(Mw_mean_weight, Mw_error_weight, '{:.2f}')
        parfile.write(f'Mw (weighted): {s}\n')

        Mo_mean = means['Mo']
        Mo_error = errors['Mo']
        s = _value_error_str(Mo_mean, Mo_error, '{:.3e}')
        parfile.write(f'Mo: {s} N·m\n')
        Mo_mean_weight = means_weight['Mo']
        Mo_error_weight = errors_weight['Mo']
        s = _value_error_str(Mo_mean_weight, Mo_error_weight, '{:.3e}')
        parfile.write(f'Mo (weighted): {s} N·m\n')

        fc_mean = means['fc']
        fc_error = errors['fc']
        s = _value_error_str(fc_mean, fc_error, '{:.3f}')
        parfile.write(f'fc: {s} Hz\n')
        fc_mean_weight = means_weight['fc']
        fc_error_weight = errors_weight['fc']
        s = _value_error_str(fc_mean_weight, fc_error_weight, '{:.3f}')
        parfile.write(f'fc (weighted): {s} Hz\n')

        t_star_mean = means['t_star']
        t_star_error = errors['t_star']
        s = _value_error_str(t_star_mean, t_star_error, '{:.3f}')
        parfile.write(f't_star: {s} s\n')
        t_star_mean_weight = means_weight['t_star']
        t_star_error_weight = errors_weight['t_star']
        s = _value_error_str(t_star_mean_weight, t_star_error_weight, '{:.3f}')
        parfile.write(f't_star (weighted): {s} s\n')

        Qo_mean = means['Qo']
        Qo_error = errors['Qo']
        s = _value_error_str(Qo_mean, Qo_error, '{:.1f}')
        parfile.write(f'Qo: {s}\n')
        try:
            Qo_mean_weight = means_weight['Qo']
            Qo_error_weight = errors_weight['Qo']
        except KeyError:
            # weighted Qo might not be computed
            Qo_mean_weight = np.nan
            Qo_error_weight = [np.nan, np.nan]
        s = _value_error_str(Qo_mean_weight, Qo_error_weight, '{:.1f}')
        parfile.write(f'Qo (weighted): {s}\n')

        ra_mean = means['radius']
        ra_error = errors['radius']
        s = _value_error_str(ra_mean, ra_error, '{:.3f}')
        parfile.write(f'Source radius: {s} m\n')
        ra_mean_weight = means_weight['radius']
        ra_error_weight = errors_weight['radius']
        s = _value_error_str(ra_mean_weight, ra_error_weight, '{:.3f}')
        parfile.write(f'Source radius (weighted): {s} m\n')

        ssd_mean = means['ssd']
        ssd_error = errors['ssd']
        s = _value_error_str(ssd_mean, ssd_error, '{:.3e}')
        parfile.write(f'Static stress drop: {s} MPa\n')
        ssd_mean_weight = means_weight['ssd']
        ssd_error_weight = errors_weight['ssd']
        s = _value_error_str(ssd_mean_weight, ssd_error_weight, '{:.3e}')
        parfile.write(f'Static stress drop (weighted): {s} MPa\n')

        Ml_mean = means.get('Ml', None)
        Ml_error = errors.get('Ml', None)
        if Ml_mean is not None:
            s = _value_error_str(Ml_mean, Ml_error, '{:.3f}')
            parfile.write(f'Ml: {s}\n')

        Er_mean = means['Er']
        Er_error = errors['Er']
        s = _value_error_str(Er_mean, Er_error, '{:.3e}')
        parfile.write(f'Er: {s} N·m\n')

        parfile.write(f'\n*** SourceSpec: {get_versions()["version"]}')
        parfile.write(
            f'\n*** Run completed on: {config.end_of_run} '
            f'{config.end_of_run_tz}')
        if getattr(config.options, 'run_id', None):
            parfile.write(f'\n*** Run ID: {config.options.run_id}')
        _write_author_and_agency_to_parfile(parfile)

    logger.info(f'Output written to file: {parfilename}')


def _dict2yaml(dict_like, level=0):
    """
    Serialize a dict-like object into YAML format.

    :param dict_like: Dict-like object to serialize.
    :type dict_like: dict
    :param level: Indentation level.
    :type level: int

    :return: YAML-formatted string.
    :rtype: str
    """
    if not isinstance(dict_like, Mapping):
        raise TypeError('dict_like must be a dict-like object')
    comments = dict_like.get('_comments', {})
    fmt = dict_like.get('_format', None)
    value_uncertainty_keys = (
        'value', 'uncertainty', 'lower_uncertainty', 'upper_uncertainty')
    target_dict = {
        key: (fmt.format(value)
              if fmt is not None and key in value_uncertainty_keys
              else value)
        for key, value in dict_like.items()
        if not key.startswith('_') and value is not None
    }
    indent = ' ' * 2 * level
    # use oneliners for dict-like objects containing value and uncertainty keys
    if set(target_dict.keys()).intersection(set(value_uncertainty_keys)):
        # make sure all values are represented as strings
        target_dict = {key: str(value) for key, value in target_dict.items()}
        oneliner = str(target_dict).replace("'", "")
        return f'{indent}{oneliner}\n'
    lines = ''
    for key, value in target_dict.items():
        if key.startswith('_') or value is None:
            continue
        if isinstance(value, Mapping):
            if level == 0:
                lines += '\n'
            with contextlib.suppress(KeyError):
                comment = comments[key]
                for line in comment.split('\n'):
                    lines += f'# {line}\n'
            lines += f'{indent}{key}:\n'
            lines += _dict2yaml(value, level + 1)
        else:
            lines += f'{indent}{key}: {value}\n'
    return lines


def _write_yaml(sspec_output):
    """
    Write sspec output in a YAML file.

    :param sspec_output: Output of the spectral inversion.
    :type sspec_output: :class:`~sourcespec.ssp_data_types.SourceSpecOutput`
    """
    if not os.path.exists(config.options.outdir):
        os.makedirs(config.options.outdir)
    evid = config.event.event_id
    yamlfilename = os.path.join(config.options.outdir, f'{evid}.ssp.yaml')
    lines = _dict2yaml(sspec_output)
    comments = sspec_output.get('_comments', {})
    with open(yamlfilename, 'w', encoding='utf-8') as fp:
        begin_comment = comments.get('begin', None)
        if begin_comment is not None:
            for line in begin_comment.split('\n'):
                fp.write(f'# {line}\n')
        fp.write(lines)
    logger.info(f'Output written to file: {yamlfilename}')


def _write_hypo71(sspec_output):
    """
    Write source parameters to hypo71 file.

    :param sspec_output: Output of the spectral inversion.
    :type sspec_output: :class:`~sourcespec.ssp_data_types.SourceSpecOutput`
    """
    if not getattr(config.options, 'hypo_file', None):
        return
    if config.hypo_file_format != 'hypo71':
        return
    with open(config.options.hypo_file, 'r', encoding='ascii') as fp:
        line = fp.readline()
        # Check if first 10 digits of the line contain characters
        if any(c.isalpha() for c in line[:10]):
            line1 = line
            line = fp.readline()
        line = list(line)
    summary_values = sspec_output.reference_values()
    mw_str = f'{summary_values["Mw"]:03.2f}'
    Ml = summary_values.get('Ml', None)
    ml_str = f'{Ml:03.2f}' if Ml is not None and ~np.isnan(Ml) else ' ' * 4
    for i in range(4):
        line[49 + i] = mw_str[0 + i]
        # line[45 + i] = mw_str[0 + i]
        line[69 + i] = ml_str[0 + i]
    outline = ''.join(line)
    evid = config.event.event_id
    hypo_file_out = os.path.join(config.options.outdir, f'{evid}.ssp.h')
    with open(hypo_file_out, 'w', encoding='ascii') as fp:
        with contextlib.suppress(Exception):
            fp.write(line1)
        fp.write(outline)
    logger.info(f'Hypo file written to: {hypo_file_out}')


def _make_symlinks():
    """Make symlinks to input files into output directory."""
    # Windows does not support symlinks
    if os.name == 'nt':
        return
    outdir = config.options.outdir
    out_data_dir = os.path.join(outdir, 'input_files')
    rel_path = os.path.relpath(config.workdir, out_data_dir)
    os.makedirs(out_data_dir, exist_ok=True)
    filelist = []
    with contextlib.suppress(AttributeError, TypeError):
        filelist += list(config.options.trace_path)
    with contextlib.suppress(AttributeError):
        filelist += [config.options.station_metadata]
    with contextlib.suppress(AttributeError):
        filelist += [config.options.hypo_file]
    with contextlib.suppress(AttributeError):
        filelist += [config.options.pick_file]
    with contextlib.suppress(AttributeError):
        filelist += [config.options.qml_file]
    with contextlib.suppress(AttributeError, TypeError):
        filelist += list(config.options.asdf_path)
    for filename in filelist:
        if filename is None or not os.path.exists(filename):
            continue
        basename = os.path.basename(filename)
        linkname = os.path.join(out_data_dir, basename)
        if os.path.isfile(linkname) or os.path.islink(linkname):
            os.remove(linkname)
        elif os.path.isdir(linkname):
            shutil.rmtree(linkname, ignore_errors=True)
        filename = os.path.join(rel_path, filename)
        os.symlink(filename, linkname)


def write_output(sspec_output):
    """
    Write results into different formats.

    :param sspec_output: Output of the spectral inversion.
    :type sspec_output: :class:`~sourcespec.ssp_data_types.SourceSpecOutput`
    """
    # Add run info to output object
    run_info = sspec_output.run_info
    run_info.SourceSpec_version = get_versions()['version']
    config.end_of_run = datetime.now()
    tz = get_localzone()
    config.end_of_run_tz = tz.tzname(config.end_of_run)
    run_info.run_completed = f'{config.end_of_run} {config.end_of_run_tz}'
    if getattr(config.options, 'run_id', None):
        run_info.run_id = config.options.run_id
    run_info.author_name = config.author_name
    run_info.author_email = config.author_email
    run_info.agency_full_name = config.agency_full_name
    run_info.agency_short_name = config.agency_short_name
    run_info.agency_url = config.agency_url
    # Symlink input files into output directory
    _make_symlinks()
    # Write to parfile (deprecated)
    _write_parfile(sspec_output)
    # Write to YAML file
    _write_yaml(sspec_output)
    # Write to SQLite database, if requested
    write_sqlite(sspec_output)
    # Write to hypo file, if requested
    _write_hypo71(sspec_output)
    # Write to quakeml file, if requested
    write_qml(sspec_output)


def save_spectra(spec_st):
    """
    Save spectra to file.

    :param spec_st: Stream of spectra.
    :type spec_st: :class:`~sourcespec.spectrum.SpectrumStream`
    """
    if not config.save_spectra:
        return
    outfile = os.path.join(
        config.options.outdir,
        f'{config.event.event_id}.spectra.hdf5'
    )
    spec_st.sort()
    for spec in spec_st:
        spec.stats.software = 'SourceSpec'
        spec.stats.software_version = get_versions()['version']
    spec_st.write(outfile)
    logger.info(f'Spectra saved to: {outfile}')
