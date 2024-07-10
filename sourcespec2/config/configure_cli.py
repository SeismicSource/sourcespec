# -*- coding: utf8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Configure SourceSpec from command line arguments and config file.

:copyright:
    2013-2024 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
import os
import sys
import types
import shutil
import uuid
import json
from copy import copy
from datetime import datetime
from .config import config
from .configobj_helpers import (
    read_config_file, parse_configspec, get_default_config_obj
)
from .library_versions import library_versions
from .configobj import ConfigObj
from .configobj.validate import Validator
from ..ssp_update_db import update_db_file


def _write_sample_config(configspec, progname):
    """
    Write a sample configuration file.

    :param configspec: The configuration specification
    :type configspec: ConfigObj
    :param progname: The name of the program
    :type progname: str
    """
    config_obj = get_default_config_obj(configspec)
    configfile = f'{progname}.conf'
    write_file = True
    if os.path.exists(configfile):
        ans = input(
            f'{configfile} already exists. Do you want to overwrite it? [y/N] '
        )
        write_file = ans in ['y', 'Y']
    if write_file:
        with open(configfile, 'wb') as fp:
            config_obj.write(fp)
        print(f'Sample config file written to: {configfile}')
        note = """
Note that the default config parameters are suited for a M<5 earthquake
recorded within ~100 km. Adjust `win_length`, `noise_pre_time`, and the
frequency bands (`bp_freqmin_*`, `bp_freqmax_*`, `freq1_*`, `freq2_*`)
according to your setup."""
        print(note)


def _update_config_file(config_file, configspec):
    """
    Update a configuration file to the latest version.

    :param config_file: The path to the configuration file
    :type config_file: str
    :param configspec: The configuration specification
    :type configspec: ConfigObj
    """
    config_obj = read_config_file(config_file, configspec)
    val = Validator()
    config_obj.validate(val)
    mod_time = datetime.fromtimestamp(os.path.getmtime(config_file))
    mod_time_str = mod_time.strftime('%Y%m%d_%H%M%S')
    config_file_old = f'{config_file}.{mod_time_str}'
    ans = input(
        f'Ok to update {config_file}? [y/N]\n'
        f'(Old file will be saved as {config_file_old}) '
    )
    if ans not in ['y', 'Y']:
        sys.exit(0)
    config_new = ConfigObj(configspec=configspec, default_encoding='utf8')
    config_new = read_config_file(None, configspec)
    config_new.validate(val)
    config_new.defaults = []
    config_new.comments = configspec.comments
    config_new.initial_comment = config_obj.initial_comment
    config_new.final_comment = configspec.final_comment
    for k, v in config_obj.items():
        if k not in config_new:
            continue
        # Fix for force_list(default=None)
        if v == ['None', ]:
            v = None
        config_new[k] = v
    migrate_options = {
        's_win_length': 'win_length',
        'traceids': 'traceid_mapping_file',
        'ignore_stations': 'ignore_traceids',
        'use_stations': 'use_traceids',
        'dataless': 'station_metadata',
        'clip_nmax': 'clip_max_percent',
        'PLOT_SHOW': 'plot_show',
        'PLOT_SAVE': 'plot_save',
        'PLOT_SAVE_FORMAT': 'plot_save_format',
        'vp': 'vp_source',
        'vs': 'vs_source',
        'rho': 'rho_source',
        'pre_p_time': 'noise_pre_time',
        'pre_s_time': 'signal_pre_time',
        'rps_from_focal_mechanism': 'rp_from_focal_mechanism',
        'paz': 'station_metadata',
        'pi_bsd_min_max': 'pi_ssd_min_max',
        'max_epi_dist': 'epi_dist_ranges'
    }
    for old_opt, new_opt in migrate_options.items():
        if old_opt in config_obj and config_obj[old_opt] != 'None':
            # max_epi_dist needs to be converted to a list
            if old_opt == 'max_epi_dist':
                config_new[new_opt] = [0, config_obj[old_opt]]
            else:
                config_new[new_opt] = config_obj[old_opt]
    shutil.copyfile(config_file, config_file_old)
    with open(config_file, 'wb') as fp:
        config_new.write(fp)
        print(f'{config_file}: updated')


def _write_config(config_obj, progname, outdir):
    """
    Write the configuration file to the output directory.

    :param config_obj: The configuration object
    :type config_obj: ConfigObj
    :param progname: The name of the program
    :type progname: str
    :param outdir: The output directory
    :type outdir: str
    """
    if progname != 'source_spec':
        return
    configfile = f'{progname}.conf'
    configfile = os.path.join(outdir, configfile)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    with open(configfile, 'wb') as fp:
        # create a copy of config_obj and remove the basemap API key
        _tmp_config_obj = copy(config_obj)
        _tmp_config_obj['plot_map_api_key'] = None
        _tmp_config_obj.write(fp)


def _init_plotting():
    """
    Initialize plotting backend.
    """
    # pylint: disable=import-outside-toplevel
    import matplotlib.pyplot as plt
    if not config.plot_show:
        plt.switch_backend('Agg')


def _init_traceid_map():
    """
    Initialize trace ID map from file.
    """
    config.TRACEID_MAP = None
    if config.traceid_mapping_file is None:
        return
    try:
        with open(config.traceid_mapping_file, 'r', encoding='utf-8') as fp:
            config.TRACEID_MAP = json.loads(fp.read())
    except Exception:
        sys.exit(
            f'traceid mapping file "{config.traceid_map_file}" not found '
            'or not in json format.\n')


def _write_sample_ssp_event_file():
    """
    Write a sample SourceSpec Event file.
    """
    ssp_event_file = 'ssp_event.yaml'
    src_path = os.path.join(
        os.path.dirname(__file__), 'ssp_event.yaml')
    dest_path = os.path.join('.', ssp_event_file)
    write_file = True
    if os.path.exists(dest_path):
        ans = input(
            f'{ssp_event_file} already exists. '
            'Do you want to overwrite it? [y/N] '
        )
        write_file = ans in ['y', 'Y']
    if write_file:
        shutil.copyfile(src_path, dest_path)
        print(f'Sample SourceSpec Event File written to: {ssp_event_file}')


def _fix_and_expand_path(path):
    """
    Fix any path issues and expand it.

    :param path: Path specification
    :type path: str
    :return: The fixed and expanded path
    :rtype: str
    """
    fixed_path = os.path.normpath(path).split(os.sep)
    fixed_path = os.path.join(*fixed_path)
    if path.startswith(os.sep):
        fixed_path = os.path.join(os.sep, fixed_path)
    elif path.startswith('~'):
        fixed_path = os.path.expanduser(fixed_path)
    return fixed_path


def configure_cli(options=None, progname='source_spec', config_overrides=None):
    """
    Configure SourceSpec from command line arguments and config file.

    The global config object is updated with the configuration parameters
    read from the configuration file and the command line options.

    :param options: An object containing command line options
    :type options: A generic object
    :param progname: The name of the program
    :type progname: str
    :param config_overrides: A dictionary with parameters that override or
        extend those defined in the config file
    :type config_overrides: dict
    """
    if options is None:
        # create an empty object to support the following getattr() calls
        options = types.SimpleNamespace()
    configspec = parse_configspec()
    if getattr(options, 'sampleconf', None):
        _write_sample_config(configspec, progname)
        sys.exit(0)
    if getattr(options, 'updateconf', None):
        _update_config_file(options.updateconf, configspec)
        sys.exit(0)
    if getattr(options, 'updatedb', None):
        update_db_file(options.updatedb)
        sys.exit(0)
    if getattr(options, 'samplesspevent', None):
        _write_sample_ssp_event_file()
        sys.exit(0)

    if getattr(options, 'config_file', None):
        options.config_file = _fix_and_expand_path(options.config_file)
        config_obj = read_config_file(options.config_file, configspec)
        # Apply overrides
        if config_overrides is not None:
            try:
                for key, value in config_overrides.items():
                    config_obj[key] = value
            except AttributeError as e:
                raise ValueError(
                    '"config_override" must be a dict-like.'
                ) from e

    # TODO: we should allow outdir to be None and not producing any output
    options.outdir = getattr(options, 'outdir', 'sspec_out')
    # Fix and expand paths in options
    options.outdir = _fix_and_expand_path(options.outdir)
    if getattr(options, 'trace_path', None):
        # trace_path is a list
        options.trace_path = [
            _fix_and_expand_path(path) for path in options.trace_path]
    if getattr(options, 'qml_file', None):
        options.qml_file = _fix_and_expand_path(options.qml_file)
    if getattr(options, 'hypo_file', None):
        options.hypo_file = _fix_and_expand_path(options.hypo_file)
    if getattr(options, 'pick_file', None):
        options.pick_file = _fix_and_expand_path(options.pick_file)

    if getattr(options, 'config_file', None):
        # Create a 'no_evid_' subdir into outdir.
        # The random hex string will make it sure that this name is unique
        # It will be then renamed once an evid is available
        hexstr = uuid.uuid4().hex
        options.outdir = os.path.join(options.outdir, f'no_evid_{hexstr}')
        _write_config(config_obj, progname, options.outdir)
        # Update config object with the contents of the config file
        config.update(config_obj.dict())
        config.running_from_command_line = True

    try:
        config.validate()
    except ValueError as msg:
        sys.exit(msg)

    # Add options to config:
    config.options = options

    # Override station_metadata config option with command line option
    if getattr(options, 'station_metadata', None):
        config.station_metadata = options.station_metadata

    # Fix and expand paths in config
    if config.database_file:
        config.database_file = _fix_and_expand_path(config.database_file)
    if config.traceid_mapping_file:
        config.traceid_mapping_file = _fix_and_expand_path(
            config.traceid_mapping_file)
    if config.station_metadata:
        config.station_metadata = _fix_and_expand_path(config.station_metadata)
    if config.residuals_filepath:
        config.residuals_filepath = _fix_and_expand_path(
            config.residuals_filepath)

    if config.plot_station_map:
        try:
            library_versions.check_cartopy_version()
            library_versions.check_pyproj_version()
            if config.plot_map_style == 'geotiff':
                library_versions.check_rasterio_version()
        except ImportError as err:
            for msg in config.warnings:
                print(msg)
            sys.exit(err)

    if config.NLL_time_dir is not None or config.NLL_model_dir is not None:
        try:
            library_versions.check_nllgrid_version()
        except ImportError as err:
            sys.exit(err)

    _init_traceid_map()
    _init_plotting()
