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
parser.add_argument('--pdb_acc', type=str)
parser.add_argument('--chain')
parser.add_argument('output')

args = parser.parse_args()
table = merge_tables(pdb_id=args.pdb_acc, chain=args.chain)
table.to_csv(args.output)