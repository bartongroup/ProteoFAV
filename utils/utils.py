#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""
--------
utils.py
--------

Defines a number of general purpose routines. Some need to be removed or
moved to more appropriated locations.

.. moduleauthor:: Fabio Madeira

:module_version: 1.0
:created_on: 01-03-2015

"""

from __future__ import print_function

import os
import re
import sys
sys.path.insert(0, '../')
import random
import time
from datetime import datetime

import requests

from config import get_config


def current_date():
    """
    Gets the current date.

    :return: outputs formatted date as 'Day/Month/Year'
    """

    date = datetime.now()
    month = date.strftime("%b")
    day = date.strftime("%d")
    year = date.strftime("%Y")

    return "{}/{}/{}".format(day, month, year)


def current_time():
    """
    Gets current date and time.

    :return: outputs formatted time as 'Day/Month/Year H:M:S'
    """

    date = datetime.now()
    year = date.strftime("%Y")
    month = date.strftime("%m")
    day = date.strftime("%d")
    hour = date.strftime("%H")
    minute = date.strftime("%M")
    second = date.strftime("%S")

    return "{}/{}/{} {}:{}:{}".format(day, month, year, hour, minute, second)


def create_directory(directory):
    """
    Creates a directory structure if it does not exist.

    :param directory: directory name (expects full path)
    :return: creates a directory if it does not exist yet
    """

    return os.makedirs(directory)


def flash(message):
    """
    Flashes a message out.

    :param message: input message str()
    """

    print(str(message))
    sys.stdout.flush()
    return


def string_split(s):
    return filter(None, re.split(r'(\d+)', s))



def write_log(message, output_file):
    """
    Appends a message to a log file.

    @param message: message to be printed
    @param output_file: path to the output file
    """

    with open(output_file, "a") as outlog:
        outlog.write(message + "\n")

    return


def get_random_string_with_n_digits(n):
    """
    gets a random number with n digits.

    @param n: number of digits
    @return: returns a random integer (int)
    """

    range_start = 10 ** (n - 1)
    range_end = (10 ** n) - 1

    return random.randint(range_start, range_end)


def request_info_url(url, params=None, verbose=False):
    """
    Gets content from the provided url.

    If it fails, tries to access it a few more times to make sure it
    was not a temporary server hang.

    @param url: input URL
    @param params: parameters to be encoded and requested
    @param verbose: boolean
    @return: returns a data object from *requests*
    """

    if params is None:
        params = {}

    # limiting ourselves to up to 10 requests per second
    time.sleep(0.1)

    try:
        request = requests.get(url, params=params)
    except Exception as e:
        message = "{}\t{}".format(url, e)
        flash(message)
        raise Exception(message)

    if verbose:
        message = "{}\t{}".format(url, request.status_code)
        flash(message)

    if request.ok:
        return request
    # Ensembl status for exceeding number of connections per second
    elif request.status_code == 429:
        time.sleep(0.5)
        request_info_url(url, params=params, verbose=verbose)
    else:
        message = "{}\t{}".format(url, request.status_code)
        flash(message)
        raise Exception(message)


def isvalid_uniprot(identifier):
    """
    'Quickly' checks if a UniProt id is valid.

    :param identifier: testing ID
    :return: boolean
    """

    config = get_config('http_uniprot')
    try:
        if identifier != '':
            request_info_url("{}{}.fasta".format(config.http_uniprot,
                                                 str(identifier)))
            return True
        else:
            raise Exception
    except Exception:
        return False


def isvalid_pdb(identifier):
    """
    'Quickly' checks if a PDB id is valid.

    :param identifier: testing ID
    :return: boolean
    """

    config = get_config('http_pdbe')
    try:
        if identifier != '':
            request_info_url("{}{}".format(config.http_pdbe,
                                           str(identifier)))
            return True
        else:
            raise Exception
    except Exception:
        return False

if __name__ == '__main__':
    # testing routines
    pass