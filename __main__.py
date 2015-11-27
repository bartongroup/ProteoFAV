#!/usr/bin/env python
# -*- coding: utf-8

"""
Command line application of the library.
"""

import argparse
import logging

from main import merge_tables

log = logging.getLogger(__name__)
usage = ''
# TODO select log level
parser = argparse.ArgumentParser(description=__doc__, usage=usage)
parser.add_argument('--pdb', type=str)
parser.add_argument('--chain')
parser.add_argument('--add_annotation', action='store_true')
parser.add_argument('output')

args = parser.parse_args()
table = merge_tables(pdb_id=args.pdb, chain=args.chain,
                     add_annotation=args.add_annotation)
table.to_csv(args.output)

