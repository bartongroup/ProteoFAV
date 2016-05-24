#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
---------
config.py
---------

This defines the methods that load and validate user defined
parameters.
Usage
>>> from proteofav.config import defaults
>>> print(defaults.api_pdbe)
http://www.ebi.ac.uk/pdbe/api/
>>> from proteofav.config import Defaults
>>> local_defaults = Defaults("config.txt")
>>> print(local_defaults.api_uniprot)
http://www.uniprot.org/uniprot/
>>> print(local_defaults.sifts_extension)
.xml.gz
>>> print(local_defaults.email)
Traceback (most recent call last):
...
AttributeError: 'Defaults' object has no attribute 'email'"""

from __future__ import print_function

import tempfile
import logging
from os import path
from ConfigParser import ConfigParser

__all__ = ["defaults", "Defaults"]

logger = logging.getLogger(__name__)


class Defaults(object):
    def __init__(self, config_file='config.txt'):
        config_file = path.join(path.dirname(__file__), config_file)
        config = ConfigParser()
        if path.isfile(config_file):
            config.read(config_file)
            self.__config = config
            self.__populate_attributes()
            self.config_file = config_file
        else:
            raise IOError('Config file {} not available.'.format(config_file))

    def __populate_attributes(self):
        for var_name, var_par in self:
            setattr(self, var_name, var_par)

    def __iter__(self):
        logged_header = False
        for section in self.__config.sections():
            for var_name, var_par in self.__config.items(section):
                if var_par == "...":
                    if not logged_header:
                        logger.warning(" Some parameters in the config.txt file need "
                                       "to be updated...")
                        logged_header = True
                    logger.warning(" Update the value for parameter {}"
                                   "...".format(var_name))
                    var_par = tempfile.gettempdir()
                yield var_name, var_par

    def update(self, config_file):
        return self.__init__(config_file)


defaults = Defaults()
