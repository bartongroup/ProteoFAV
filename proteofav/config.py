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
>>> print(local_defaults.http_uniprot)
http://www.uniprot.org/uniprot/
>>> print(local_defaults.sifts_extension)
.xml.gz
>>> print(local_defaults.email)
Traceback (most recent call last):
...
AttributeError: 'Defaults' object has no attribute 'email'"""

from __future__ import print_function

import logging
from os import path
from ConfigParser import ConfigParser

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
        logged_header = False
        for var_name, var_par in self:
            if var_par == "...":
                if not logged_header:
                    logger.warning(" Some parameters in the config.txt file need "
                                   "to be provided...")
                    logged_header = True
                logger.warning(" Update the value for parameter {}"
                               "...".format(var_name))
            setattr(self, var_name, var_par)

    def __iter__(self):
        for section in self.__config.sections():
            for var_name, var_par in self.__config.items(section):
                yield var_name, var_par

defaults = Defaults()
