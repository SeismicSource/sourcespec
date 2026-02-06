# -*- coding: utf8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Config class for sourcespec.

:copyright:
    2013-2026 Claudio Satriano <satriano@ipgp.fr>
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
    mandatory_config_params, check_deprecated_config_params
)
from .configobj import ConfigObj
from .configobj.validate import Validator
from ..ssp_parse_arguments import (
    _get_description, _init_parser, _update_parser
)


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
    if input_list is None or input_list == ['None', ] or input_list == [None]:
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


def _none_length(input_list):
    """
    Return the length of input list, or 1 if input list is None

    :param input_list: Input list or None
    :type input_list: list or None

    :return: List length or 1
    :rtype: int
    """
    return 1 if input_list is None else len(input_list)


def _format_value_list(value):
    """
    Format a list or tuple value for display.

    :param value: List or tuple to format
    :return: Formatted string representation
    """
    if len(value) == 0:
        return '[]'
    formatted_items = [str(item) for item in value[:3]]
    result = '[' + ', '.join(formatted_items)
    if len(value) > 3:
        result += f', ... ({len(value)} total)'
    result += ']'
    return result


def _format_value(value, max_length=100):
    """
    Format a configuration value for display.

    :param value: The value to format
    :param max_length: Maximum length before truncating
    :return: Formatted string representation
    """
    if value is None:
        return '<None>'
    if isinstance(value, str):
        return (
            f'{value[:max_length]}...'
            if len(value) > max_length
            else value
        )
    if isinstance(value, (list, tuple)):
        return _format_value_list(value)
    if isinstance(value, dict):
        return f'<dict with {len(value)} items>'
    return str(value)


def _generate_html_repr(
    obj, title='Configuration', key_order=None, internal_keys=None,
    help_texts=None
):
    """
    Generate HTML representation for a dict-like object.

    :param obj: Dictionary-like object to display
    :param title: Title for the HTML display
    :param key_order: Optional list of keys in desired order
    :param internal_keys: Optional list of internal keys to separate
    :param help_texts: Optional dict mapping keys to help text strings
    :return: HTML string
    """
    if help_texts is None:
        help_texts = {}
    if internal_keys is None:
        internal_keys = []

    # Generate unique ID for this instance
    import random  # pylint: disable=import-outside-toplevel
    import string  # pylint: disable=import-outside-toplevel
    unique_id = ''.join(
        random.choices(string.ascii_lowercase, k=8)
    )

    style_div = (
        'font-family: monospace; border: 1px solid #ddd; '
        'padding: 10px; border-radius: 5px; '
        'background-color: #f9f9f9; color: #000;'
    )
    html = f'<div style="{style_div}">'
    html += (
        f'<h4 style="margin-top: 0; color: #000;">'
        f'{title}</h4>'
    )

    # Add search input
    search_style = (
        'width: 100%; padding: 8px; margin-bottom: 10px; '
        'border: 1px solid #ccc; border-radius: 4px; '
        'box-sizing: border-box; font-size: 14px;'
    )
    html += (
        f'<input type="text" id="search_{unique_id}" '
        f'placeholder="Search parameters..." '
        f'style="{search_style}">'
    )

    # Scrollable table container
    style_scroll = (
        'max-height: 400px; overflow-y: auto; '
        'border: 1px solid #ddd; border-radius: 3px;'
    )
    html += f'<div style="{style_scroll}">'
    html += (
        f'<table id="table_{unique_id}" '
        f'style="border-collapse: collapse; width: auto;">'
    )

    # Use provided key_order, or fall back to sorted keys
    keys_to_display = (
        key_order
        if key_order is not None
        else sorted(obj.keys())
    )
    # Track whether we've added the section headers
    internal_section_added = False
    config_section_added = False

    for key in keys_to_display:
        # Add section header before transitioning from internal to config keys
        if internal_keys:
            is_internal = key in internal_keys
            if is_internal and not internal_section_added:
                # Add "Internal Configuration" header
                style = (
                    'background-color: #f0f0f0; '
                    'border-bottom: 2px solid #999;'
                )
                cell_style = (
                    'padding: 8px; font-weight: bold; '
                    'color: #333; background-color: #f0f0f0;'
                )
                html += (
                    f'<tr style="{style}">'
                    f'<td colspan="2" style="{cell_style}">'
                    'Internal Configuration</td></tr>'
                )
                internal_section_added = True
            elif (
                not is_internal and not config_section_added
                and internal_section_added
            ):
                # Add "Configuration Parameters" header
                style = (
                    'background-color: #f0f0f0; '
                    'border-bottom: 2px solid #999;'
                )
                cell_style = (
                    'padding: 8px; font-weight: bold; '
                    'color: #333; background-color: #f0f0f0;'
                )
                html += (
                    f'<tr style="{style}">'
                    f'<td colspan="2" style="{cell_style}">'
                    'Configuration Parameters</td></tr>'
                )
                config_section_added = True

        if key not in obj:
            continue

        value = obj[key]
        formatted = _format_value(value, max_length=200)

        # Escape HTML special characters
        key_html = (
            key.replace('&', '&amp;')
            .replace('<', '&lt;')
            .replace('>', '&gt;')
        )
        value_html = (
            formatted.replace('&', '&amp;')
            .replace('<', '&lt;')
            .replace('>', '&gt;')
        )

        # Use slightly different styling for internal vs config keys
        if key in internal_keys:
            bg_color = '#fafafa'
            key_color = '#555'
        else:
            bg_color = '#fff'
            key_color = '#000'

        tr_style = (
            f'border-bottom: 1px solid #eee; '
            f'background-color: {bg_color};'
        )
        html += f'<tr class="data-row" style="{tr_style}">'
        style_key = (
            'padding: 8px; font-weight: bold; width: 30%; '
            f'vertical-align: top; color: {key_color};'
        )
        html += f'<td style="{style_key}">{key_html}</td>'
        style_val = (
            f'padding: 8px; color: {key_color}; text-align: left;'
        )
        html += f'<td style="{style_val}">{value_html}'

        if help_text := help_texts.get(key):
            # Replace lines containing only "#" with empty lines
            help_lines = [
                '' if line.strip() == '#' else line
                for line in help_text.split('\n')
            ]
            processed_help = '\n'.join(help_lines)
            help_html = (
                processed_help.replace('&', '&amp;')
                .replace('<', '&lt;')
                .replace('>', '&gt;')
                .replace('\n', '<br>')
                .replace('  ', '&nbsp;&nbsp;')
            )
            help_style = (
                'display: block; margin-top: 4px; '
                'font-size: 0.85em; color: #666; '
                'font-style: italic; white-space: pre-wrap;'
            )
            html += (
                f'<div style="{help_style}">{help_html}</div>'
            )

        html += '</td>'
        html += '</tr>'

    html += '</table>'
    html += '</div>'

    # Add JavaScript for search functionality
    html += f'''
<script>
(function() {{
    const searchInput = document.getElementById('search_{unique_id}');
    const table = document.getElementById('table_{unique_id}');
    const rows = table.getElementsByClassName('data-row');
    if (searchInput) {{
        searchInput.addEventListener('input', function() {{
            const searchText = this.value.toLowerCase();
            for (let i = 0; i < rows.length; i++) {{
                const row = rows[i];
                const text = row.textContent.toLowerCase();
                if (searchText === '' || text.includes(searchText)) {{
                    row.style.display = '';
                }} else {{
                    row.style.display = 'none';
                }}
            }}
        }});
    }}
}})();
</script>
'''

    html += '</div>'
    return html

# ---- End Helper functions ----


class _Options(dict):
    """
    Options class for sourcespec, with builtin checks for API users.
    """
    def __init__(self):
        """Initialize the Options object with default values."""
        super().__init__()
        self.progname = 'source_spec'
        self.set_defaults(self.progname)

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

    def _get_parser_actions(self):
        """Get the argparse actions for the given program name."""
        description, epilog, nargs = _get_description(self.progname)
        parser = _init_parser(description, epilog, nargs)
        _update_parser(parser, self.progname)
        return [
            action
            for action in getattr(parser, '_actions', [])
            if action.dest not in ('help', 'version')
        ]

    def set_defaults(self, progname='source_spec'):
        """
        Populate Options object with all possible keys, set to None, False
        or default values

        :param progname: name of program, 'source_spec' or 'source_model'
        :type progname: str
        """
        if progname not in ('source_spec', 'source_model'):
            raise ValueError(
                'progname must be "source_spec" or "source_model"'
            )
        self.progname = progname
        actions = self._get_parser_actions()
        for action in actions:
            self[action.dest] = action.default

    def get_help(self, option=None):
        """
        Print information about given option

        :param option: name of option or None to list all valid options
        :type option: str
        """
        actions = self._get_parser_actions()
        # get a list of valid options
        valid_options = [action.dest for action in actions]
        for action in actions:
            if action.dest == option:
                print(action.help, end='')
                return
        if option is not None:
            print(f'Unknown option "{option}"')
        valid_options_str = '\n'.join(valid_options)
        print(f'Valid options:\n{valid_options_str}', end='')

    def __repr__(self):
        """Return a detailed string representation of options."""
        lines = ['_Options(']
        # Get keys in the order they appear in argparse actions
        actions = self._get_parser_actions()
        action_dests = [action.dest for action in actions]
        # Use action order, then sort any remaining keys
        keys_to_display = [
            key for key in action_dests if key in self
        ] + [
            key for key in sorted(self.keys())
            if key not in action_dests
        ]
        for key in keys_to_display:
            value = self[key]
            formatted = _format_value(value)
            lines.append(f'  {key}: {formatted}')
        lines.append(')')
        return '\n'.join(lines)

    def __str__(self):
        """Return a readable string representation of options."""
        return self.__repr__()

    def _repr_html_(self):
        """
        Return HTML representation for Jupyter notebooks.

        This method is automatically called by Jupyter when displaying
        the options object in a notebook cell.
        """
        # Get keys in the order they appear in argparse actions
        actions = self._get_parser_actions()
        action_dests = [action.dest for action in actions]
        key_order = [
            key for key in action_dests if key in self
        ] + [
            key for key in sorted(self.keys())
            if key not in action_dests
        ]
        return _generate_html_repr(
            self, title='SourceSpec Options', key_order=key_order
        )


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
        """Initialize the config object with default values."""
        super().__init__()
        self._set_defaults()

    def _set_defaults(self):
        """Set the default values for the config object."""
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
        self['INSTR_CODES_VEL'] = ['H', 'L', 'P']
        self['INSTR_CODES_ACC'] = ['N', ]
        # Initialize config object to the default values
        configspec = parse_configspec()
        config_obj = get_default_config_obj(configspec)
        self.update(config_obj.dict())
        # Store the key order from configspec for use in repr methods
        # Internal keys go first, then configspec keys in order
        self._internal_keys = [
            'running_from_command_line', 'vertical_channel_codes',
            'horizontal_channel_codes_1', 'horizontal_channel_codes_2',
            'TRACEID_MAP', 'options', 'warnings', 'figures', 'workdir',
            'INSTR_CODES_VEL', 'INSTR_CODES_ACC'
        ]
        self._key_order = self._internal_keys + list(config_obj.keys())

    def reset(self):
        """Reset the config object to the default values."""
        super().clear()
        self._set_defaults()

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

    def _get_help_text(self, parameter):
        """
        Get help text for a configuration parameter.

        :param parameter: name of parameter
        :type parameter: str
        :return: Help text string or empty string if not found
        :rtype: str
        """
        configspec = parse_configspec()
        config_obj = get_default_config_obj(configspec)
        if parameter not in config_obj:
            return ''
        comment = config_obj.comments.get(parameter)
        if len(comment) == 0 or comment == ['']:
            comment = config_obj.inline_comments.get(parameter) or []
        # Clean up comment lines by removing leading '# ', empty lines,
        # and lines containing '----'. Join the remaining lines into a
        # single string.
        comment = '\n'.join(
            line.replace('# ', '') for line in comment
            if '----' not in line and line.strip()
        )
        return comment

    def get_help(self, parameter=None):
        """
        Print information about given configuration parameter
        (i.e., associated comment in template configuration file)

        :param parameter: name of parameter
            or None to list all valid parameters
        :type parameter: str
        """
        configspec = parse_configspec()
        config_obj = get_default_config_obj(configspec)
        if parameter not in config_obj:
            if parameter is not None:
                print(f'Unknown parameter "{parameter}"')
            valid_params = '\n'.join(config_obj.comments.keys())
            print(
                'Please specify one of the following parameters:\n'
                f'{valid_params}',
                end=''
            )
            return
        comment = self._get_help_text(parameter)
        print(comment)

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
        # Make sure that self.options is still an _Options object
        if not isinstance(self['options'], _Options):
            _options = self['options']
            self['options'] = _Options()
            self['options'].update(_options)
        # Check deprecated and mandatory parameters
        self._check_deprecated_config_params()
        self._check_mandatory_config_params()
        self._check_force_list()
        self._check_list_lengths()
        self._check_Er_freq_range()
        self._parse_free_surface_amplification()
        self._check_html_report()

    def _check_deprecated_config_params(self):
        """
        Check the deprecated configuration parameters.

        :raises ValueError: If an error occurs while parsing the parameters
        """
        deprecation_msgs = check_deprecated_config_params(self)
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
                'vp',
                'vs',
                'rho',
                'layer_top_depths',
                'vp_source',
                'vs_source',
                'rho_source',
                'layer_top_depths_source',
            ]:
                self[param] = _float_list(self[param])
        except ValueError as msg:
            raise ValueError(
                f'Error parsing parameter "{param}": {msg}'
            ) from msg

    def _check_list_lengths(self):
        """
        Check that the lists describing the velocity and density models
        have the same length.

        :raises ValueError: If an error occurs while parsing the parameters
        """
        # Check that the tt velocity models have the same length
        n_vp = _none_length(config.vp)
        n_vs = _none_length(config.vs)
        n_rho = _none_length(config.rho)
        n_layer_top_depths = _none_length(config.layer_top_depths)
        try:
            assert n_vp == n_vs == n_rho == n_layer_top_depths
        except AssertionError as err:
            raise ValueError(
                'Error: "vp", "vs", "rho", and "layer_top_depths" '
                'must have the same length.'
            ) from err
        # Check that the velocity and density models have the same length
        n_vp_source = _none_length(config.vp_source)
        n_vs_source = _none_length(config.vs_source)
        n_rho_source = _none_length(config.rho_source)
        n_layer_top_depths_source = _none_length(
            config.layer_top_depths_source)
        try:
            assert (
                n_vp_source == n_vs_source == n_rho_source ==
                n_layer_top_depths_source
            )
        except AssertionError as err:
            raise ValueError(
                'Error: "vp_source", "vs_source", "rho_source", and '
                '"layer_top_depths_source" must have the same length.'
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

    def _parse_free_surface_amplification(self):
        """
        Parse free_surface_amplification config parameter.
        """
        fsa = self['free_surface_amplification']
        fsa_items = ()
        # First, check if values are in the form:
        # ((STATID1, FLOAT), (STATID2, FLOAT), ...)
        if (
            isinstance(fsa, (list, tuple))
            and all(
                isinstance(item, (list, tuple)) and len(item) == 2
                for item in fsa
            )
        ):
            with contextlib.suppress(ValueError, TypeError):
                fsa_items = tuple(
                    (fsa_item[0], float(fsa_item[1]))
                    for fsa_item in fsa
                )
        # If the above failed, check if single FLOAT value is given
        if not fsa_items and len(fsa) == 1:
            with contextlib.suppress(ValueError, TypeError):
                fsa_items = (('*', float(fsa[0])),)
        # If the above failed, or if there are multiple elements,
        # try to parse as STATID: FLOAT pairs
        if not fsa_items:
            # Values must be in the form: STATID1: FLOAT, STATID2: FLOAT, ...
            try:
                if any(':' not in item for item in fsa):
                    # empty ValueError to jump to except clause
                    raise ValueError()
                fsa_items = tuple(
                    (key.strip(), float(val)) for key, val in
                    (item.split(':') for item in fsa)
                )
            except (ValueError, TypeError) as e:
                raise ValueError(
                    'Invalid format. '
                    'Expected: STATID1: FLOAT, STATID2: FLOAT, ...'
                ) from e
        # sanity check
        if not fsa_items:
            raise ValueError(
                'Invalid format. '
                'Expected: STATID1: FLOAT, STATID2: FLOAT, ...'
            )
        # Check that all values are positive
        for _, val in fsa_items:
            if val <= 0:
                raise ValueError('All values must be positive.')
        self['free_surface_amplification'] = fsa_items

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

    def _format_value(self, value, max_length=100):
        """
        Format a configuration value for display.

        :param value: The value to format
        :param max_length: Maximum length before truncating
        :return: Formatted string representation
        """
        return _format_value(value, max_length)

    def _format_value_list(self, value):
        """Format a list or tuple value for display."""
        return _format_value_list(value)

    def __repr__(self):
        """Return a detailed string representation of the config."""
        lines = ['_Config(']
        # Use stored key order, or fall back to sorted keys
        key_order = getattr(self, '_key_order', sorted(self.keys()))
        for key in key_order:
            if key in self:
                value = self[key]
                formatted = _format_value(value)
                lines.append(f'  {key}: {formatted}')
        lines.append(')')
        return '\n'.join(lines)

    def __str__(self):
        """Return a readable string representation of the config."""
        return self.__repr__()

    def _repr_html_(self):
        """
        Return HTML representation for Jupyter notebooks.

        This method is automatically called by Jupyter when displaying
        the config object in a notebook cell.
        """
        # Use stored key order for consistent display
        key_order = getattr(self, '_key_order', None)
        if key_order:
            # Filter to only include keys that exist in the config
            key_order = [k for k in key_order if k in self]
        # Get internal keys for visual separation
        internal_keys = getattr(self, '_internal_keys', [])
        # Gather help texts for all config parameters
        help_texts = {}
        for key in self.keys():
            if key not in internal_keys:
                if help_text := self._get_help_text(key):
                    help_texts[key] = help_text
        return _generate_html_repr(
            self, title='SourceSpec Configuration',
            key_order=key_order, internal_keys=internal_keys,
            help_texts=help_texts
        )


# Global config object, initialized with default values
# API users should use this object to access configuration parameters
# and update them as needed
config = _Config()
