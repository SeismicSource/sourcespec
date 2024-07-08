# -*- coding: utf8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Config class for sourcespec.

:copyright:
    2013-2025 Claudio Satriano <satriano@ipgp.fr>
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
import contextlib
from copy import copy
from datetime import datetime
from collections import defaultdict
from .library_versions import library_versions
from .mandatory_deprecated import (
    mandatory_config_params, deprecated_config_params
)
from .configobj import ConfigObj
from .configobj.validate import Validator
from ..ssp_update_db import update_db_file

# TODO: remove these when the global config object will be implemented
INSTR_CODES_VEL = []
INSTR_CODES_ACC = []


class Config(dict):
    """Config class for sourcespec."""
    def __init__(self):
        # Additional config values. Tey must be defined using the dict syntax.
        self['running_from_command_line'] = False
        self['vertical_channel_codes'] = ['Z']
        self['horizontal_channel_codes_1'] = ['N', 'R']
        self['horizontal_channel_codes_2'] = ['E', 'T']
        # Empty options object, for compatibility with the command line version
        self['options'] = types.SimpleNamespace()
        # A list of warnings to be issued when logger is set up
        self['warnings'] = []
        # Create a dict to store figure paths
        self['figures'] = defaultdict(list)
        # store the absolute path of the current working directory
        self['workdir'] = os.getcwd()
        # SEED standard instrument codes:
        # https://ds.iris.edu/ds/nodes/dmc/data/formats/seed-channel-naming/
        self['INSTR_CODES_VEL'] = ['H', 'L']
        self['INSTR_CODES_ACC'] = ['N', ]
        # Initialize config object to the default values
        configspec = _parse_configspec()
        config_obj = _get_default_config_obj(configspec)
        self.update(config_obj.dict())

    def __setitem__(self, key, value):
        """Make Config keys accessible as attributes."""
        super().__setattr__(key, value)
        super().__setitem__(key, value)

    def __getattr__(self, key):
        """Make Config keys accessible as attributes."""
        try:
            return self.__getitem__(key)
        except KeyError as err:
            raise AttributeError(err) from err

    __setattr__ = __setitem__

    def update(self, other):
        """
        Update the configuration with the values from another dictionary.

        :param dict other: The dictionary with the new values

        :raises ValueError: If an error occurs while parsing the options
        """
        for key, value in other.items():
            self[key] = value
        # Set to None all the 'None' strings
        for key, value in self.items():
            if value == 'None':
                self[key] = None
        # Make sure that self['figures'] is still a defaultdict
        self['figures'] = defaultdict(list, self['figures'])
        self._update_channel_codes()
        self._update_instrument_codes()

    def _update_channel_codes(self):
        """
        Update channel codes with mis-oriented channels.
        """
        msc = self.mis_oriented_channels
        if msc is None:
            return
        self['vertical_channel_codes'].append(msc[0])
        self['horizontal_channel_codes_1'].append(msc[1])
        self['horizontal_channel_codes_2'].append(msc[2])
        self['vertical_channel_codes'] =\
            list(set(self['vertical_channel_codes']))
        self['horizontal_channel_codes_1'] =\
            list(set(self['horizontal_channel_codes_1']))
        self['horizontal_channel_codes_2'] =\
            list(set(self['horizontal_channel_codes_2']))

    def _update_instrument_codes(self):
        """
        Update instrument codes from user-defined values.
        """
        # User-defined instrument codes:
        instr_code_acc_user = self['instrument_code_acceleration']
        instr_code_vel_user = self['instrument_code_velocity']
        # Remove user-defined instrument codes if they conflict
        # with another instrument
        with contextlib.suppress(ValueError):
            self['INSTR_CODES_VEL'].remove(instr_code_acc_user)
        with contextlib.suppress(ValueError):
            self['INSTR_CODES_ACC'].remove(instr_code_vel_user)
        # Add user-defined instrument codes
        if instr_code_vel_user is not None:
            self['INSTR_CODES_VEL'].append(instr_code_vel_user)
        if instr_code_acc_user is not None:
            self['INSTR_CODES_ACC'].append(instr_code_acc_user)
        self['INSTR_CODES_VEL'] = list(set(self['INSTR_CODES_VEL']))
        self['INSTR_CODES_ACC'] = list(set(self['INSTR_CODES_ACC']))

    def validate(self):
        """
        Validate the configuration.

        :raises ValueError: If an error occurs while validating the options
        """
        config_obj = ConfigObj(self, configspec=_parse_configspec())
        val = Validator()
        test = config_obj.validate(val)
        # The variable "test" is:
        # - True if everything is ok
        # - False if no config value is provided
        # - A dict if invalid values are present,
        #   with the invalid values as False
        msg = ''
        if isinstance(test, dict):
            for entry in [e for e in test if not test[e]]:
                msg += f'\nInvalid value for "{entry}": "{config_obj[entry]}"'
            raise ValueError(msg)
        if not test:
            raise ValueError('No configuration value present!')
        # Reupdate the config object with the validated values, which include
        # type conversion from string to numeric values
        self.update(config_obj.dict())
        # Check deprecated and mandatory parameters
        self._check_deprecated_config_params()
        self._check_mandatory_config_params()
        self._check_force_list()
        self._check_list_lengths()
        self._check_Er_freq_range()
        self._check_html_report()

    def _check_deprecated_config_params(self):
        deprecation_msgs = []
        for param, msgs in deprecated_config_params.items():
            if param in self:
                # add two spaces before each line of the message
                msg = ''.join(f'  {line}\n' for line in msgs)
                # replace first character with '>'
                msg = f'>{msg[1:]}'
                deprecation_msgs.append(msg)
        if not deprecation_msgs:
            return
        msg = ''
        if self['running_from_command_line']:
            msg += (
                'Error: your config file contains deprecated parameters:\n\n'
            )
        msg += ''.join(deprecation_msgs)
        if self['running_from_command_line']:
            msg += (
                '\nPlease upgrade your config file manually or '
                'via the "-U" option.\n'
            )
        raise ValueError(msg)

    def _check_mandatory_config_params(self):
        messages = []
        for par in mandatory_config_params:
            if self[par] is None:
                msg = f'"{par}" is mandatory and cannot be None'
                messages.append(msg)
        if messages:
            msg = '\n'.join(messages)
            raise ValueError(msg)

    def _check_force_list(self):
        """
        Check the force_list options and convert them to lists of floats.

        :raises ValueError: If an error occurs while parsing the options
        """
        try:
            for param in [
                'vp_source', 'vs_source', 'rho_source', 'layer_top_depths'
            ]:
                self[param] = _float_list(self[param])
        except ValueError as msg:
            raise ValueError(
                f'Error parsing parameter "{param}": {msg}'
            ) from msg

    def _check_list_lengths(self):
        """
        Check that the lists describing the source model have the same length.
        """
        n_vp_source = _none_lenght(self.vp_source)
        n_vs_source = _none_lenght(self.vs_source)
        n_rho_source = _none_lenght(self.rho_source)
        n_layer_top_depths = _none_lenght(self.layer_top_depths)
        try:
            assert n_vp_source == n_vs_source == n_rho_source \
                == n_layer_top_depths
        except AssertionError as err:
            raise ValueError(
                'Error: "vp_source", "vs_source", "rho_source", and '
                '"layer_top_depths" must have the same length.'
            ) from err

    def _check_Er_freq_range(self):
        """
        Check the Er_freq_range option.

        :raises ValueError: If an error occurs while parsing the options
        """
        if self['Er_freq_range'] is None:
            self['Er_freq_range'] = [None, None]
        try:
            self['Er_freq_range'] = _float_list(
                self['Er_freq_range'], max_length=2,
                accepted_values=[None, 'noise']
            )
        except ValueError as msg:
            raise ValueError(
                f'Error parsing parameter "Er_freq_range": {msg}'
            ) from msg

    def _check_html_report(self):
        """
        Check the html_report option.
        """
        if self['html_report']:
            if not self['plot_save']:
                self['warnings'].append(
                    'The "html_report" option is selected but "plot_save" '
                    'is "False". HTML report will have no plots.'
                )
            if self['plot_save_format'] not in ['png', 'svg']:
                self['warnings'].append(
                    'The "html_report" option is selected but '
                    '"plot_save_format" is not "png" or "svg". '
                    'HTML report will have no plots.'
                )


def _read_config_file(config_file, configspec=None):
    kwargs = {
        'configspec': configspec,
        'file_error': True,
        'default_encoding': 'utf8'
    }
    if configspec is None:
        kwargs.update({
            'interpolation': False,
            'list_values': False,
            '_inspec': True
        })
    try:
        config_obj = ConfigObj(config_file, **kwargs)
    except IOError as err:
        sys.stderr.write(f'{err}\n')
        sys.exit(1)
    except Exception as err:
        sys.stderr.write(f'Unable to read "{config_file}": {err}\n')
        sys.exit(1)
    return config_obj


def _parse_configspec():
    configspec_file = os.path.join(
        os.path.dirname(__file__), 'configspec.conf')
    return _read_config_file(configspec_file)


def _get_default_config_obj(configspec):
    config_obj = ConfigObj(configspec=configspec, default_encoding='utf8')
    val = Validator()
    config_obj.validate(val)
    config_obj.defaults = []
    config_obj.initial_comment = configspec.initial_comment
    config_obj.comments = configspec.comments
    config_obj.final_comment = configspec.final_comment
    return config_obj


def _write_sample_config(configspec, progname):
    config_obj = _get_default_config_obj(configspec)
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
    config_obj = _read_config_file(config_file, configspec)
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
    config_new = _read_config_file(None, configspec)
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


def _init_plotting(plot_show):
    # pylint: disable=import-outside-toplevel
    import matplotlib.pyplot as plt
    if not plot_show:
        plt.switch_backend('Agg')


def _init_traceid_map(config):
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

    :param str path: Path specification
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


def _float_list(input_list, max_length=None, accepted_values=None):
    """
    Convert an input list to a list of floats.

    :param list input_list: Input list or None
    :return: A list of floats or None
    :rtype: list
    """
    if input_list is None:
        return None
    if accepted_values is None:
        accepted_values = []

    def _parse_float(val):
        val = None if val == 'None' else val
        return val if val in accepted_values else float(val)

    try:
        return [_parse_float(val) for val in input_list[:max_length]]
    except ValueError as e:
        raise ValueError('Cannot parse all values in list') from e


def _none_lenght(input_list):
    """
    Return the length of input list, or 1 if input list is None

    :param list input_list: Input list or None
    :return: List length or 1
    :rtype: int
    """
    return 1 if input_list is None else len(input_list)


# Global config object, initialized with default values
# API users should use this object to access configuration parameters
# and update them as needed
config = Config()


def configure(options=None, progname='source_spec', config_overrides=None):
    """
    Parse command line arguments and read config file.

    :param object options: An object containing command line options
    :param str progname: The name of the program
    :param dict config_overrides: A dictionary with parameters that override or
        extend those defined in the config file
    :return: A ``Config`` object with both command line and config options.
    """
    if options is None:
        # create an empty object to support the following getattr() calls
        options = types.SimpleNamespace()
    configspec = _parse_configspec()
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
        config_obj = _read_config_file(options.config_file, configspec)
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

    # TODO: remove these when the global config object will be implemented
    global INSTR_CODES_VEL
    global INSTR_CODES_ACC
    INSTR_CODES_VEL = config.INSTR_CODES_VEL
    INSTR_CODES_ACC = config.INSTR_CODES_ACC

    _init_traceid_map(config)
    _init_plotting(config.plot_show)
    return config
