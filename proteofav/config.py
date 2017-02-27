#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Defines the methods that load and validate user configuration parameters, such
    as data resources web address or local and remote file paths.

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
AttributeError: 'Defaults' object has no attribute 'email'

"""

from __future__ import print_function

import tempfile
import logging
from os import path

import sys

try:
    # python 2.7
    from ConfigParser import ConfigParser
except ImportError:
    from configparser import ConfigParser

__all__ = ["defaults", "Defaults"]
log = logging.getLogger(__name__)
logging.captureWarnings(True)
logging.basicConfig(stream=sys.stderr, level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s ')


class Defaults(object):
    """
    Container for configuration parameters. The object attributes receive
        the values parsed from the config file. Can also be set on the run.

    >>> from proteofav.config import defaults
    >>> print(defaults.api_pdbe)
    http://www.ebi.ac.uk/pdbe/api/
    >>> defaults.api_pdbe = 'test'
    >>> print(defaults.api_pdbe)
    test #

    """
    def __init__(self, config_file='config.txt'):
        config_file = path.join(path.dirname(__file__), config_file)
        config = ConfigParser()
        if path.isfile(config_file):
            config.read(config_file)
            self.__config = config
            self.__populate_attributes()
            self.config_file = config_file
        else:     # pragma: no cover
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
                        log.warning(" Some parameters in the config.txt file need "
                                       "to be updated...")
                        logged_header = True
                    log.warning(" Update the value for parameter {}"
                                   "...".format(var_name))
                    var_par = tempfile.gettempdir()
                yield var_name, var_par

    def update(self, config_file):
        return self.__init__(config_file)


defaults = Defaults()
