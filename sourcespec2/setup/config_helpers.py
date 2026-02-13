# -*- coding: utf8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Helper functions for the config and option classes.

:copyright:
    2013-2026 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""


def float_list(input_list, max_length=None, accepted_values=None):
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


def none_length(input_list):
    """
    Return the length of input list, or 1 if input list is None

    :param input_list: Input list or None
    :type input_list: list or None

    :return: List length or 1
    :rtype: int
    """
    return 1 if input_list is None else len(input_list)


def format_value_list(value):
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


def format_value(value, max_length=100):
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
        return format_value_list(value)
    if isinstance(value, dict):
        return f'<dict with {len(value)} items>'
    return str(value)


def generate_html_repr(
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
        formatted = format_value(value, max_length=200)

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
            .replace('$', '<span>$</span>')  # avoid math rendering
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
                .replace('$', '<span>$</span>')  # avoid math rendering
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
        // Prevent notebook keyboard shortcuts from interfering
        // with search input
        searchInput.addEventListener('keydown', function(event) {{
            // Don't stop propagation for Ctrl/Cmd combinations
            // (e.g., Ctrl+A, Ctrl+C)
            // These should work normally in the input field
            if (event.ctrlKey || event.metaKey) {{
                return;
            }}
            // Stop all other keys from reaching the notebook
            event.stopPropagation();
        }});
        searchInput.addEventListener('keyup', function(event) {{
            event.stopPropagation();
        }});
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
