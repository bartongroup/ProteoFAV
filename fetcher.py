#!/usr/bin/env python
# -*- coding: utf-8
"""
Created on 11:41 24/07/15 2015
Small routine to fetch data files from their respective repositries.
"""
import logging
import gzip
import shutil
import os
import urllib
from config import defaults
import socket

socket.setdefaulttimeout(15)
log = logging.getLogger(__name__)
__author__ = 'tbrittoborges'


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
    if isinstance(sources, str):
        sources = [sources]
    if not directory:
        directory = defaults.temp
    elif isinstance(directory, str):
        if not os.path.isdir(directory):
            raise IOError(directory + " is not a valid path.")
        directory = [directory] * len(sources)
    elif hasattr(directory, "__iter__"):
        if len(directory) != len(sources):
            raise IOError(directory + " you need one directory for each source,"
                                      "or a path for all.")
        for d in directory:
            if not os.path.isdir(d):
                raise IOError(d + " is not a valid path.")
    else:
        raise TypeError('Unidentified source|directory combination.')

    result = []
    for source, destination in zip(sources, directory):
        filename = identifier + getattr(defaults, source + '_extension')
        url = getattr(defaults, source + '_fetch') + filename
        try:
            urllib.urlretrieve(url, destination + filename)
        except IOError as e:
            log.error('Unable to retrieve {} for {}'.format(url, str(e)))
            raise
        if filename.endswith('.gz'):
            with gzip.open(destination + filename, 'rb') as input_f, \
                    open(destination + filename.replace('.gz', ''),
                         'wb') as output_f:
                shutil.copyfileobj(input_f, output_f)
                filename = filename.replace('.gz', '')
        result.append(destination + filename)
    return result


if __name__ == '__main__':
    pass
