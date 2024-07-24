# -*- coding: utf8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Config class for sourcespec.

:copyright:
    2013-2024 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
import os
import contextlib
import warnings
from collections import defaultdict
from .configobj_helpers import parse_configspec, get_default_config_obj
from .mandatory_deprecated import (
    mandatory_config_params, deprecated_config_params
)
from .configobj import ConfigObj
from .configobj.validate import Validator


# ---- Helper functions ----
def _float_list(input_list, max_length=None, accepted_values=None):
    """
    Convert an input list to a list of floats.

    :param input_list: Input list or None
    :type input_list: list or None
    :param max_length: Maximum length of the list
    :type max_length: int or None
    :param accepted_values: List of accepted values
    :type accepted_values: list or None

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

    :param input_list: Input list or None
    :type input_list: list or None

    :return: List length or 1
    :rtype: int
    """
    return 1 if input_list is None else len(input_list)
# ---- End Helper functions ----


class _Options(dict):
    """
    Options class for sourcespec, with builtin checks for API users.
    """
    def __setitem__(self, key, value):
        """
        Make Config keys accessible as attributes.

        Perform specific checks for some parameters.
        """
        if (
            key in ['trace_path', 'asdf_path'] and
            (value is not None and not isinstance(value, (list, tuple)))
        ):
            warnings.warn(
                f'"{key}" must be a list. Converting to a list with one '
                'element.'
            )
            value = [value]
        super().__setattr__(key, value)
        super().__setitem__(key, value)

    def __getattr__(self, key):
        """Make Config keys accessible as attributes."""
        try:
            return self.__getitem__(key)
        except KeyError as err:
            raise AttributeError(err) from err

    __setattr__ = __setitem__


class _Config(dict):
    """
    Config class for sourcespec.

    This class stores the configuration parameters for sourcespec.
    Parameters can be accessed as dictionary keys (e.g. config['param'])
    or as attributes (e.g. config.param).
    The class is initialized with default values.

    .. note::

        This class is private and should not be used directly.
        Import the global config object instead.
    """
    def __init__(self):
        # Additional config values that must exist for the code to run without
        # errors. They must be defined using the dict syntax.
        self['running_from_command_line'] = False
        self['vertical_channel_codes'] = ['Z']
        self['horizontal_channel_codes_1'] = ['N', 'R']
        self['horizontal_channel_codes_2'] = ['E', 'T']
        self['TRACEID_MAP'] = None
        # options object with some built-in checks for API users
        self['options'] = _Options()
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
        configspec = parse_configspec()
        config_obj = get_default_config_obj(configspec)
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

    def clear(self):
        """Clear the configuration and the options."""
        self['options'].clear()
        super().clear()

    def update(self, other):
        """
        Update the configuration with the values from another dictionary.

        :param other: The dictionary with the new values
        :type other: dict

        :raises ValueError: If an error occurs while parsing the parameters
        """
        for key, value in other.items():
            self[key] = value
        # Set to None all the 'None' strings
        for key, value in self.items():
            if value == 'None':
                self[key] = None
            if isinstance(value, list):
                self[key] = [None if val == 'None' else val for val in value]
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

        :raises ValueError: If an error occurs while validating the parameters
        """
        config_obj = ConfigObj(self, configspec=parse_configspec())
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
        """
        Check the deprecated configuration parameters.

        :raises ValueError: If an error occurs while parsing the parameters
        """
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
        """
        Check the mandatory configuration parameters.

        :raises ValueError: If an error occurs while parsing the parameters
        """
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
        Check the force_list parameters and convert them to lists of floats.

        :raises ValueError: If an error occurs while parsing the parameters
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

        :raises ValueError: If an error occurs while parsing the parameters
        """
        n_vp_source = _none_lenght(self['vp_source'])
        n_vs_source = _none_lenght(self['vs_source'])
        n_rho_source = _none_lenght(self['rho_source'])
        n_layer_top_depths = _none_lenght(self['layer_top_depths'])
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
        Check the Er_freq_range parameter.

        :raises ValueError: If an error occurs while parsing the parameters
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
        Check the html_report parameter.
        """
        if self['html_report']:
            if not self['plot_save']:
                self['warnings'].append(
                    'The "html_report" parameter is selected but "plot_save" '
                    'is "False". HTML report will have no plots.'
                )
            if self['plot_save_format'] not in ['png', 'svg']:
                self['warnings'].append(
                    'The "html_report" parameter is selected but '
                    '"plot_save_format" is not "png" or "svg". '
                    'HTML report will have no plots.'
                )


# Global config object, initialized with default values
# API users should use this object to access configuration parameters
# and update them as needed
config = _Config()
