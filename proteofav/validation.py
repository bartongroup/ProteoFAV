#!/usr/bin/env python
# -*- coding: utf-8 -*-


import os
import logging
import numpy as np
import pandas as pd
from lxml import etree

from proteofav.config import defaults
from proteofav.utils import (row_selector, fetch_files,
                             constrain_column_types, exclude_columns)
from proteofav.library import validation_types


log = logging.getLogger('proteofav.config')

__all__ = ['parse_validation_residues', 'select_validation']


def parse_validation_residues(filename, excluded_cols=None,
                              global_parameters=False,
                              fix_label_alt_id=True, fix_ins_code=True,):
    """
    Parse the PDB's validation file to a pandas DataFrame.

    :param filename: path to the Validation file
    :param excluded_cols: list of columns to be excluded
    :param global_parameters: bool
    :param fix_label_alt_id: boolean
    :param fix_ins_code: boolean
    :return: returns a pandas DataFrame
    """

    tree = etree.parse(filename)
    root = tree.getroot()
    if global_parameters:
        global_parameters = root.find('Entry').attrib
        log.info(global_parameters)
    rows = []
    header = set()
    for i, elem in enumerate(root.iterfind('ModelledSubgroup')):
        rows.append(dict(elem.attrib))
        header.update(rows[-1].keys())
    for row in rows:
        not_in = {k: None for k in header.difference(row.keys())}
        row.update(not_in)

    table = pd.DataFrame(rows, columns=header)

    # fixes the 'icode'
    if fix_ins_code:
        table = _fix_pdb_ins_code(table)
    # fixes the 'altcode
    if fix_label_alt_id:
        table = _fix_label_alt_id(table)

    # exclude columns
    table = exclude_columns(table, excluded=excluded_cols)

    # enforce some specific column types
    table = constrain_column_types(table, col_type_dict=validation_types)
    # replace 'NoneTypes' by NaNs
    table.fillna(value=np.nan, inplace=True)

    if table.empty:
        log.error('Validation file {} resulted in a empty Dataframe'.format(filename))
        raise ValueError('Validation file {} resulted in a empty Dataframe'.format(
            filename))
    return table


def _fix_pdb_ins_code(table):
    """
    Utility that fixes the 'pdbx_PDB_ins_code' column to match what is expected
    in the mmCIF format.

    :param table: pandas DataFrame object
    :return: returns a modified pandas DataFrame
    """

    table['icode'] = table['icode'].str.replace('\ |', '?')
    table['icode'] = table['icode'].fillna('?').astype(str)
    return table


def _fix_label_alt_id(table):
    """
    Utility that fixes the 'label_alt_id' column to match what is
    expected in the mmCIF format.

    :param table: pandas DataFrame object
    :return: returns a modified pandas DataFrame
    """

    table['altcode'] = table['altcode'].str.replace('\ |\?', '.')
    table['altcode'] = table['altcode'].fillna('.').astype(str)
    return table


def select_validation(pdb_id, chains=None):
    """
    Produces table from PDB validation XML file.

    :param pdb_id: PDB identifier
    :param chains: PDB protein chain
    :return: pandas dataframe
    """
    val_path = os.path.join(defaults.db_validation, pdb_id + defaults.validation_extension)
    try:
        val_table = parse_validation_residues(val_path)
    except IOError:
        val_path = fetch_files(pdb_id, sources='validation', directory=defaults.db_pdb)[0]
        val_table = parse_validation_residues(val_path)

    if chains:
        val_table = row_selector(val_table, 'chain', chains)
    if not val_table.empty:
        val_table.columns = ["validation_" + name for name in val_table.columns]
        return val_table
    raise ValueError('Parsing {} resulted in a empty table'.format(val_path))
