#!/usr/bin/env python2.7
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
import random
from datetime import datetime


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

    if not os.path.exists(directory):
        os.makedirs(directory)
    return


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



if __name__ == '__main__':
    # testing routines
    print(current_time())
    pass