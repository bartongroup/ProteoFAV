#!/usr/bin/env python
# -*- coding: utf-8

__author__ = 'tbrittoborges'
"""
Created on 11:41 24/07/15 2015 

"""
import os
import urllib
from config import defaults

def fetch_files(identifier, directory=None, sources=("cif", "dssp", "sifts")):
    """
    Use defaults to fetch files from servers. Defaults server are defined in
     the default config.txt. Returns None since it has side effects.
    :param identifier: protein identifier as PDB identifier
    :param directory: path to download. The default downloads all files to
     default.temp folder. If is a string and valid path downloads all files to
     that folder. If its a iterable and all valid path downloads to the
     respective folders
    :param sources: where to fetch the data. Must be in the config file.
    :return: None
    """
    if not directory:
        directory = defaults.temp
    elif isinstance(directory, str):
        if not os.path.isdir(directory):
            raise IOError(directory + " is not a valid path.")
        directory = [directory] * len(sources)
    elif hasattr(directory, "__iter__"):
        if len(directory) != len((sources)):
            raise IOError(directory + " you need one diretory for each source,"
                                      "or a path for all.")
        for d in directory:
            if not os.path.isdir(d):
                raise IOError(d + " is not a valid path.")
    #else: FAIL?

    for source, destination in zip(sources, directory):
        # test if destination is writtable?
        url = getattr(defaults, 'fetch_' + source)
        file_name =  getattr(defaults, 'extension_' + source)
        urllib.urlretrieve(url, destination + file_name)
    # TODO assert raise meaninfull error
    # TODO test fetching from test dir

if __name__ == '__main__':
    fetch_files("3mn5", "/Users/tbrittoborges/Downloads/")
