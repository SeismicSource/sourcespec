# -*- coding: utf8 -*- 
class Config(dict):
    #The following is to make Config keys accessible as attributes
    def __setitem__(self, key, value):
	    super(Config, self).__setattr__(key, value)
	    super(Config, self).__setitem__(key, value)
    def __getattr__(self, key):
        try:
            return self.__getitem__(key)
        except KeyError, message:
            raise AttributeError, message
    __setattr__ = __setitem__
