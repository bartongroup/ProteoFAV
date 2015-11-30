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
parser = argparse.ArgumentParser(description=__doc__, usage=usage)
parser.add_argument('--pdb', type=str)
parser.add_argument('--chain')
parser.add_argument('--sp_accession')
parser.add_argument('output', required=True)

parser.add_argument('-v', '--verbose', action='store_true')
parser.add_argument('-q', '--quiet', action='store_true')
parser.add_argument('-l', '--log') # set log to file.

parser.add_argument('--add_annotation', action='store_true')
parser.add_argument('--add_validation', action='store_true')
parser.add_argument('--add_variants', action='store_true')
parser.add_argument('--add_annotations', action='store_true')

args = parser.parse_args()
table = merge_tables(pdb_id=args.pdb, chain=args.chain,
                     add_annotation=args.add_annotation)
table.to_csv(args.output)

