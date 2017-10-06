# -*- coding: utf-8 -*-

"""
Defines the methods that load and validate user configuration parameters, such
    as data resources web address or local and remote file paths.

>>> from proteofav.config import defaults
>>> print(defaults.api_pdbe)
http://www.ebi.ac.uk/pdbe/api/
>>> from proteofav.config import Defaults
>>> local_defaults = Defaults('config.ini')
>>> print(local_defaults.api_uniprot)
http://www.uniprot.org/uniprot/
>>> print(local_defaults.email)
Traceback (most recent call last):
...
AttributeError: 'Defaults' object has no attribute 'email'

"""

import os
import sys
import click
import logging
import tempfile

try:
    # python 2.7
    from ConfigParser import ConfigParser
except ImportError:
    from configparser import ConfigParser

log = logging.getLogger(__name__)
logging.captureWarnings(True)
logging.basicConfig(stream=sys.stderr, level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s ')

__all__ = ["defaults", "Defaults"]


class Defaults(object):
    """
    Container for configuration parameters. The object attributes receive
        the values parsed from the config file. Can also be set on the run.

    >>> from proteofav.config import defaults
    >>> print(defaults.api_pdbe)
    http://www.ebi.ac.uk/pdbe/api/
    >>> defaults.api_pdbe = 'test'
    >>> print(defaults.api_pdbe)
    test
    """

    def __init__(self, config_file=None):
        if config_file:
            # user provided config
            pass
        elif os.path.isfile(os.path.join(click.get_app_dir('proteofav'), 'config.ini')):
            # os config
            config_file = os.path.join(click.get_app_dir('proteofav'), 'config.ini')
        else:
            # proteofav default config
            config_file = os.path.join(os.path.dirname(__file__), 'config.ini')
        config = ConfigParser()
        if os.path.isfile(config_file):
            config.read(config_file)
            self.__config = config
            self.populate_attributes()
            self.config_file = config_file
        else:  # pragma: no cover
            raise IOError('Config file {} not available.'.format(config_file))

    def populate_attributes(self):
        for var_name, var_par in self:
            setattr(self, var_name, var_par)

    def __iter__(self):
        for section in self.__config.sections():
            for var_name, var_par in self.__config.items(section):
                if var_name.startswith('db') and var_par == "...":
                    var_par = tempfile.gettempdir()
                yield var_name, var_par

    def update(self, config_file):
        return self.__init__(config_file)

    def commit_configuration(self):
        pass

    def write(self, file_path=None):
        file_path = file_path or os.path.join(click.get_app_dir('proteofav'), 'config.ini')

        if not os.path.exists(os.path.dirname(file_path)):
            os.makedirs(os.path.dirname(file_path))
        # self.reverse_attributes() # TODO reverse populate from __dict__ to self.config
        with open(file_path, 'w') as f:
            self.__config.write(f)


defaults = Defaults()
