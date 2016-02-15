#!/usr/bin/env python
# -*- coding: utf-8

"""
Command line application of the library.
"""

import logging
import sys
import argparse

from .main import merge_tables

log = logging.getLogger(__name__)
usage = ''
parser = argparse.ArgumentParser(description=__doc__, usage=usage)
parser.add_argument('--pdb', type=str)
parser.add_argument('--chain')
parser.add_argument('--sp_accession')
parser.add_argument('output')

parser.add_argument('-v', '--verbose', action='store_true')
parser.add_argument('-q', '--quiet', action='store_true')
parser.add_argument('-l', '--log', default=sys.stderr)

parser.add_argument('--add_annotation', action='store_true')
parser.add_argument('--add_validation', action='store_true')
parser.add_argument('--add_variants', action='store_true')
parser.add_argument('--remove_redundant', action='store_true')

args = parser.parse_args()
logger_level = 10 if args.verbose else 20
logger_level = 50 if args.quite else 20
logging.basicConfig(stream=args.log, level=logger_level)


log.setLevel()
table = merge_tables(pdb_id=args.pdb, chain=args.chain,
                     add_annotation=args.add_annotation,
                     add_validation=args.add_validation,
                     add_variants=args.add_variants,
                     remove_redundant=args.remove_redundant)

table.to_csv(args.output)

