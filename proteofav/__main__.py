#!/usr/bin/env python
# -*- coding: utf-8

"""
Command line application of the library.
"""
import argparse
import logging
import sys

from .main import merge_tables

log = logging.getLogger(__name__)
usage = """Default logger is is set to show warnings."""
parser = argparse.ArgumentParser(description=__doc__, usage=usage)
parser.add_argument('--pdb', type=str)
parser.add_argument('--chain')
parser.add_argument('--sp_accession')
parser.add_argument('output')
parser.add_argument('--type', choices=['csv', 'json', 'tab'], default='csv', dest="out_type")
# TODO add JALVIEW
# TODO add CHIMERA
parser.add_argument('-v', '--verbose',
                    action='store_true',
                    const=logging.DEBUG,
                    default=logging.WARNING,
                    help="Set logger to debug messages")
parser.add_argument('-q', '--quiet',
                    action='store_true',
                    const=logging.INFO,
                    dest="loglevel",
                    help="Set logger to error only messages")
parser.add_argument('-l', '--log',  default=sys.stderr)
parser.add_argument('--add_annotation', action='store_true')
parser.add_argument('--add_validation', action='store_true')
parser.add_argument('--add_variants', action='store_true')
parser.add_argument('--remove_redundant', action='store_true')

args = parser.parse_args()

logging.basicConfig(stream=args.log, level=args.loglevel,
                    format='%(asctime)s - %(levelname)s - %(message)s ')

table = merge_tables(pdb_id=args.pdb, chain=args.chain,
                     add_annotation=args.add_annotation,
                     add_validation=args.add_validation,
                     add_ensembl_variants=args.add_variants,
                     drop_empty_cols=args.remove_redundant)

if args.out_type == 'json':
    table.to_json(args.output)
elif args.out_type == 'tab':
    table.to_table(args.output)
else:
    table.to_csv(args.output)