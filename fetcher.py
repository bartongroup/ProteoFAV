#!/usr/bin/env python
# -*- coding: utf-8

__author__ = 'tbrittoborges'
"""
Created on 11:41 24/07/15 2015 

"""
import urllib
from config import defaults

def fetch_files(identifier, dir=defaults.temp, sources=("cif", "dssp", "sifts")):
    """
    Use defaults to fetch files from servers. Defaults server are defined in
     the default config.txt. Returns None since it has side effects.
    :param identifier: protein identifier as PDB identifier
    :param dir: path to download
    :param sources: where to fetch the data. Must be in the config file.
    :return: None
    """
    for source in sources:
        url = getattr(defaults, 'fetch_' + source)
        file_name =  getattr(defaults, 'extension_' + source)
        urllib.urlretrieve(url, dir + file_name)
    # TODO assert raise meaninfull error
    # TODO test fetching from test dir


if __name__ == '__main__':
    fetch_files("3mn5", "/Users/tbrittoborges/Downloads/")
