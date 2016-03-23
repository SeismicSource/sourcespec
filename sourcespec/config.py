# -*- coding: utf8 -*-
#
# (c) 2013-2016 Claudio Satriano <satriano@ipgp.fr>


class Config(dict):
    """Config class for source_spec."""

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
