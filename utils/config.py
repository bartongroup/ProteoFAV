#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-

#!/sw/opt/pypy-2.4.0-src/pypy/goal/pypy-c


"""
---------
config.py
---------

This defines the methods that load and validate user defined
parameters.

.. moduleauthor:: Fabio Madeira

:module_version: 1.0
:created_on: 26-02-2015

"""

from __future__ import print_function

import os
import sys
from ConfigParser import SafeConfigParser

from utils import flash

# edit this name if you rename your config file
CONFIG_FILE = "config.txt"


def load_config(option):
    """
    Loads and tests input parameters.

    :param option: option name
    :return: returns a valid config value for the inputed option
    """

    parser = SafeConfigParser()
    try:
        filename = "{}/{}".format(os.getcwd(), CONFIG_FILE)
        parser.read(filename)

        # global paths
        if option == "tmp_dir":
            return parser.get('Global', 'tmp_dir')
        if option == "db_root":
            return parser.get('Global', 'db_root')
        elif option == "db_data":
            return parser.get('Global', 'db_data')
        elif option == "db_pdb":
            return parser.get('Global', 'db_pdb')
        elif option == "db_cif":
            return parser.get('Global', 'db_cif')
        elif option == "db_sifts":
            return parser.get('Global', 'db_sifts')
        elif option == "db_dssp":
            return parser.get('Global', 'db_dssp')
        elif option == "db_uniprot":
            return parser.get('Global', 'db_uniprot')
        elif option == "db_results":
            return parser.get('Global', 'db_results')
        elif option == "db_biolip":
            return parser.get('Global', 'db_biolip')
        elif option == "db_cath":
            return parser.get('Global', 'db_cath')

        # addresses for fetching/downloading data
        elif option == "api_pdbe":
            return parser.get('Addresses', 'api_pdbe')
        elif option == "api_cath":
            return parser.get('Addresses', 'api_cath')
        elif option == "web_cath":
            return parser.get('Addresses', 'web_cath')
        elif option == "api_rcsb":
            return parser.get('Addresses', 'api_rcsb')
        elif option == "rsync_pdb":
            return parser.get('Addresses', 'rsync_pdb')
        elif option == "rsync_cif":
            return parser.get('Addresses', 'rsync_cif')
        elif option == "rsync_sifts":
            return parser.get('Addresses', 'rsync_sifts')
        elif option == "rsync_dssp":
            return parser.get('Addresses', 'rsync_dssp')
        elif option == "ftp_sifts":
            return parser.get('Addresses', 'ftp_sifts')
        elif option == "ftp_obsolete":
            return parser.get('Addresses', 'ftp_obsolete')
        elif option == "http_uniprot":
            return parser.get('Addresses', 'http_uniprot')
        else:
            raise ValueError("ERROR: Invalid option")
    except IOError:
        flash("ERROR: Invalid Config File")


if __name__ == "__main__":
    print("Testing...")
