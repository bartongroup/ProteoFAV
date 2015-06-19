#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-


"""
--------------
test_main_parser.py
--------------

This module contains a set of commands initializing standard
:py:class:`argparse.ArgumentParser` objects with standard sets of pre-defined
options. The idea here is that I have a standard basic parser with set syntax
but can also have a 'cluster parser' with a set of pre-defined cluster oriented
options that can be added to the standard basic parser for scripts that need
cluster support, etc etc

.. moduleauthor:: Nick Schurch, modified by Fabio Madeira

:module_version: 1.2
:created_on: 2013-04-08
:modified_on: 28-02-2015

----------------
"""

from __future__ import print_function

import argparse
import tempfile


import os
import main_callable as cc


def standard_parser(ver, prog=None, usage=None, description=None,
                    epilog=None, tmpdir=True, infile=True, infiletype="txt",
                    outfile=True):
    """Set up a command line parser with standard options.

    Depending on the options supplied to the function the standard options
    include an input file, an output file, a log file, a temporary directory
    a verbosity switch and standard :py:class:`argparse.ArgumentParser` version
    and help switches.

    Returns a :py:class:`argparse.ArgumentParser` instance and two lists; the
    first is a list of the positional arguments and their defaults (which
    should be None), the second is a list of keyword arguments and their
    defaults. These lists are used by :ref:`standard_logging.py
    <standard_logging_autodocs>` to give clarity to the logged script options.
    Note, the destination in the resulting namespace (and in the first value of
    tuples in the keyword arguments lists) is the same as the long option for
    all optional arguments; lets try to keep it that way!
    """

    # use these to store some metadata on the options added to the parser.
    # These will be lists of tuples with the argument name, and default value.
    pos_args = []
    kw_args = []

    # setup the argparse parser
    formatter = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(prog=prog,
                                     usage=usage,
                                     description=description,
                                     epilog=epilog,
                                     formatter_class=formatter,
                                     add_help=False)

    reqarggroup = parser.add_argument_group('Standard Options (required)')

    ########## input file #########
    if infile:
        infilehelp = "Specify an input %s file (inc. path if different from " \
                     "the current working directory) to be consumed by this " \
                     "script." % infiletype

        reqarggroup.add_argument('infile',
                                 action='store',
                                 type=cc.input_file,
                                 help=infilehelp)
        pos_args.append(('infile', None))

    ########## output file #########
    if outfile:
        outfilehelp = "Specify an output file (inc. path if different from " \
                      "the current working directory) to be generated by " \
                      "this script."

        reqarggroup.add_argument('outfile',
                                 action='store',
                                 type=cc.output_file,
                                 help=outfilehelp)
        pos_args.append(('outfile', None))

    ########## log file #########
    loghelp = "Specify a log file (inc. path if different from the current " \
              "working directory) of the log file generated by this script."

    reqarggroup.add_argument('-l', '--log',
                             action='store',
                             dest='log',
                             required=True,
                             type=cc.output_file,
                             help=loghelp)

    optarggroup = parser.add_argument_group('Standard Options (optional)')
    kw_args.append(('log', None))

    ########## tmpdir #########

    if tmpdir:
        tmpdirhelp = "Specify a directory to use as a temp dir for this " \
                     "script. If --tmpdir is listed without a directory, or " \
                     "is omitted entirely, then a system-generated tmpdir " \
                     "will be used."

        tmp_dir = tempfile.mkdtemp()
        optarggroup.add_argument('--tmpdir',
                                 nargs='?',
                                 action='store',
                                 dest='tmpdir',
                                 const=tmp_dir,
                                 default=tmp_dir,
                                 type=cc.output_path,
                                 help=tmpdirhelp)

        kw_args.append(('tmpdir', tmp_dir))

    ########## version, verbose, help #########
    optarggroup.add_argument('--version',
                             action='version',
                             dest='version',
                             version='%(prog)s' + ' %s' % str(ver),
                             help='Show program\'s version number and exit')

    optarggroup.add_argument('-v', '--verbose',
                             action='store_true',
                             dest='verbose',
                             help='Verbosity switch for logging and warnings')

    kw_args.append(('verbose', False))

    optarggroup.add_argument('-h', '--help',
                             action='help',
                             help="Show this help message and exit")

    return parser, pos_args, kw_args


def get_std_req_group(parser):
    """ Returns the 'Standard Options (required)' argument group from the
        standard parser. """

    for group in parser._action_groups:
        if group.title == "Standard Options (required)":
            return group

    return None


def get_std_opt_group(parser):
    """ Returns the 'Standard Options (optional)' argument group from the
        standard parser. """

    for group in parser._action_groups:
        if group.title == "Standard Options (optional)":
            return group

    return None


if __name__ == "__main__":
    print("Testing...")