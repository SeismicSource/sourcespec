# -*- coding: utf8 -*-
#
# (c) 2013-2014 Claudio Satriano <satriano@ipgp.fr>

class Config(dict):
    '''
    Config class for source_spec.
    '''
    def __setitem__(self, key, value):
        '''
        Makes Config keys accessible as attributes
        '''
        super(Config, self).__setattr__(key, value)
        super(Config, self).__setitem__(key, value)
    def __getattr__(self, key):
        try:
            return self.__getitem__(key)
        except KeyError, message:
            raise AttributeError, message
    __setattr__ = __setitem__
