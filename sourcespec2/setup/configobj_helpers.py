# -*- coding: utf8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Helper functions for using ConfigObj in SourceSpec.

:copyright:
    2013-2025 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
import os
import sys
from .configobj import ConfigObj
from .configobj.validate import Validator


def read_config_file(config_file, configspec=None):
    """
    Read a configuration file and return a ConfigObj object.

    :param config_file: path to the configuration file
    :type config_file: str
    :param configspec: path to the configuration specification file
    :type configspec: str

    :return: ConfigObj object
    :rtype: ConfigObj
    """
    kwargs = {
        'configspec': configspec,
        'file_error': True,
        'default_encoding': 'utf8'
    }
    if configspec is None:
        kwargs |= {
            'interpolation': False,
            'list_values': False,
            '_inspec': True,
        }
    try:
        config_obj = ConfigObj(config_file, **kwargs)
    except IOError as err:
        sys.stderr.write(f'{err}\n')
        sys.exit(1)
    except Exception as err:
        sys.stderr.write(f'Unable to read "{config_file}": {err}\n')
        sys.exit(1)
    return config_obj


def parse_configspec():
    """
    Parse the configuration specification file and return a ConfigObj object.

    :return: ConfigObj object
    :rtype: ConfigObj
    """
    configspec_file = os.path.join(
        os.path.dirname(__file__), 'configspec.conf')
    return read_config_file(configspec_file)


def get_default_config_obj(configspec):
    """
    Return a ConfigObj object with default values.

    :param configspec: ConfigObj object with the configuration specification
    :type configspec: ConfigObj

    :return: ConfigObj object
    :rtype: ConfigObj
    """
    config_obj = ConfigObj(configspec=configspec, default_encoding='utf8')
    val = Validator()
    config_obj.validate(val)
    config_obj.defaults = []
    config_obj.initial_comment = configspec.initial_comment
    config_obj.comments = configspec.comments
    config_obj.final_comment = configspec.final_comment
    return config_obj
