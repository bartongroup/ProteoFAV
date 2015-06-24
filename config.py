#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
---------
config.py
---------

This defines the methods that load and validate user defined
parameters.
Usage
>>> from config import defaults
>>> print(defaults.api_pdbe)
http://www.ebi.ac.uk/pdbe/api/
>>> from config import Defaults
>>> local_defaults = Defaults("config_local.txt")
>>> print(local_defaults.db_pdb)
/Users/tbrittoborges/Downloads
>>> print(local_defaults.contact_email)
tbrittoborges@dundee.ac.uk
>>> print(local_defaults.email)
Traceback (most recent call last):
...
AttributeError: 'Defaults' object has no attribute 'email'"""

from __future__ import print_function
from ConfigParser import ConfigParser
import logging
from os import path

__all__ = ["defaults", "Defaults"]

logger = logging.getLogger(__name__)


class Defaults(object):
    def __init__(self, config_file=None):
        default_config = path.join(path.dirname(__file__), "config.txt")
        config = ConfigParser()
        config_default = config_file or default_config
        config.read(config_default)
        self.__config = config
        self.__populate_attributes()

    def __populate_attributes(self):
        for section in self.__config.sections():
            for var_name, var_par in self.__config.items(section):
                if var_par == "...":
                    var_par = self._manual_filling(var_name)
                setattr(self, var_name, var_par)

    def _manual_filling(self, var_name):
        """Offers user to write in config.txt"""
        # raise NotImplementedError
        return "..."

defaults = Defaults()
