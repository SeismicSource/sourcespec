# -*- coding: utf8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Config class for sourcespec.

:copyright:
    2013-2022 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""


class Config(dict):
    """Config class for sourcespec."""

    def __setitem__(self, key, value):
        """Make Config keys accessible as attributes."""
        super(Config, self).__setattr__(key, value)
        super(Config, self).__setitem__(key, value)

    def __getattr__(self, key):
        """Make Config keys accessible as attributes."""
        try:
            return self.__getitem__(key)
        except KeyError as err:
            raise AttributeError(err)

    __setattr__ = __setitem__
