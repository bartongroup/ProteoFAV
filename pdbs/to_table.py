#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Created on 03/06/2015

"""

import sys
sys.path.insert(0, '../')
import os
import logging
from StringIO import StringIO

import pandas as pd

logger = logging.getLogger(__name__)



def _dssp_to_table(filename):
    """
    Loads and parses DSSP files generating a pandas dataframe.

    :param filename: input SIFTS xml file
    :return: pandas table dataframe
    """

    if not os.path.isfile(filename):
        raise IOError('File {} not found or unavailable.'.format(filename))

    # column width descriptors
    cols_widths = ((0, 5), (6, 10), (11, 12), (13, 14), (16, 17), (35, 38),
                   (103, 109), (109, 115))
    # simplified headers for the table
    dssp_header = ("dssp_index", "icode", "chain_id", "aa", "ss", "acc", "phi",
                   "psi")
    return pd.read_fwf(filename, skiprows=28, names=dssp_header,
                       colspecs=cols_widths, index_col=0, compression=None)


def _mmcif_atom_to_table(filename, delimiter=None):
    """
    Testing a loader of mmCIF ATOM lines with pandas.

    :param filename: input CIF file
    :return: pandas table dataframe
    """

    if not os.path.isfile(filename):
        raise IOError('File {} not found or unavailable.'.format(filename))

    _header_mmcif = []
    lines = []
    with open(filename) as inlines:
        for line in inlines:
            if line.startswith("_atom_site."):
                _header_mmcif.append(line.split('.')[1].rstrip())
            elif line.startswith("ATOM"):
                lines.append(line)
            elif line.startswith("HETATM"):
                lines.append(line)
    lines = "".join(lines)

    if delimiter is None:
        return pd.read_table(StringIO(lines), delim_whitespace=True,
                             names=_header_mmcif, compression=None)
    else:
        return pd.read_table(StringIO(lines), sep=str(delimiter),
                             names=_header_mmcif, compression=None)


if __name__ == '__main__':
    # testing routines
    pass
