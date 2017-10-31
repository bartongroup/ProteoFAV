# -*- coding: utf-8 -*-

import os
import logging
import numpy as np
import pandas as pd
from lxml import etree

from proteofav.config import defaults
from proteofav.utils import (row_selector, constrain_column_types,
                             exclude_columns, Downloader, GenericInputs)
from proteofav.library import validation_types

log = logging.getLogger('proteofav.config')

__all__ = ['parse_validation_residues', 'select_validation',
           'filter_validation', 'download_validation', 'Validation']


def parse_validation_residues(filename, excluded_cols=None, global_parameters=False,
                              fix_label_alt_id=True, fix_ins_code=True):
    """
    Parse the PDB's validation file to a pandas DataFrame.

    :param filename: path to the Validation file
    :param excluded_cols: list of columns to be excluded
    :param global_parameters: bool
    :param fix_label_alt_id: boolean
    :param fix_ins_code: boolean
    :return: returns a pandas DataFrame
    """

    log.debug("Parsing Validation records from file...")

    tree = etree.parse(filename)
    root = tree.getroot()
    if global_parameters:
        global_parameters = root.find('Entry').attrib
        log.debug(global_parameters)
    rows = []
    header = set()
    for i, elem in enumerate(root.iterfind('ModelledSubgroup')):
        rows.append(dict(elem.attrib))
        header.update(rows[-1].keys())
    for row in rows:
        not_in = {k: None for k in header.difference(row.keys())}
        row.update(not_in)

    table = pd.DataFrame(rows, columns=header)

    # column renaming
    table.columns = ["validation_" + name for name in table.columns]

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
        raise ValueError('Validation file {} resulted in a empty Dataframe'
                         ''.format(filename))
    return table


def _fix_pdb_ins_code(table):
    """
    Utility that fixes the 'pdbx_PDB_ins_code' column to match what is expected
    in the mmCIF format.

    :param table: pandas DataFrame object
    :return: returns a modified pandas DataFrame
    """

    table['validation_icode'] = table['validation_icode'].str.replace('\ |', '?')
    table['validation_icode'] = table['validation_icode'].fillna('?').astype(str)
    return table


def _fix_label_alt_id(table):
    """
    Utility that fixes the 'label_alt_id' column to match what is
    expected in the mmCIF format.

    :param table: pandas DataFrame object
    :return: returns a modified pandas DataFrame
    """

    table['validation_altcode'] = table['validation_altcode'].str.replace('\ |\?', '.')
    table['validation_altcode'] = table['validation_altcode'].fillna('.').astype(str)
    return table


def _add_validation_res_full(table):
    """
    Utility that adds a new column to the table.
    Adds a new column with the 'full res' (i.e. seq_id + ins_code).

    :param table: pandas DataFrame object
    :return: returns a modified pandas DataFrame
    """

    # adds both 'label' and 'auth' entries
    if 'validation_resnum' in table and 'validation_icode' in table:
        table['validation_resnum_full'] = (table['validation_resnum'] +
                                           table['validation_icode'].str.replace('?', ''))
    return table


def select_validation(identifier, excluded_cols=None, overwrite=False, **kwargs):
    """
    Produces table from PDB validation XML file.

    :param identifier: PDB/mmCIF accession ID
    :param excluded_cols: option to exclude mmCIF columns
    :param overwrite: boolean
    :return: returns a pandas DataFrame
    """
    filename = os.path.join(defaults.db_validation,
                            "{}_validation.xml".format(identifier))

    download_validation(identifier=identifier, filename=filename,
                        overwrite=overwrite)

    table = parse_validation_residues(filename=filename, excluded_cols=excluded_cols)

    table = filter_validation(table, excluded_cols=excluded_cols, **kwargs)
    table = constrain_column_types(table, col_type_dict=validation_types)
    return table


def filter_validation(table, excluded_cols=None, chains=None, res=None,
                      add_res_full=True):
    """
    Filter for Validation Pandas Dataframes.

    :param table: pandas DataFrame object
    :param excluded_cols: option to exclude Validation columns
    :param chains: (tuple) chain IDs or None
    :param res: (tuple) res IDs or None
    :param add_res_full: option to extend the table with 'res_full'
    :return: returns a pandas DataFrame
    """

    # selections / filtering
    # excluding columns
    table = exclude_columns(table, excluded=excluded_cols)

    # table modular extensions or selections
    if add_res_full:
        table = _add_validation_res_full(table)
        log.debug("Validation added full res (res + ins_code)...")

    # excluding rows
    if chains is not None:
        table = row_selector(table, 'validation_chain', chains)
        log.debug("Validation table filtered by CHAIN...")

    if res is not None:
        table = row_selector(table, 'validation_resnum', res)
        log.debug("Validation table filtered by RES...")

    if table.empty:
        raise ValueError("The filters resulted in an empty DataFrame...")
    return table


def download_validation(identifier, filename, overwrite=False):
    """
    Downloads a Validation Data XML from the PDBe.

    :param identifier: (str) PDB accession ID
    :param filename: path to the Validation file
    :param overwrite: (boolean)
    :return: (side effects) output file path
    """

    url_root = defaults.validation_fetch
    url_endpoint = "{}_validation.xml".format(identifier)
    url = url_root + url_endpoint
    Downloader(url=url, filename=filename,
               decompress=False, overwrite=overwrite)


class Validation(GenericInputs):
    def read(self, filename=None, **kwargs):
        filename = self._get_filename(filename)
        self.table = parse_validation_residues(filename=filename, **kwargs)
        return self.table

    def download(self, identifier=None, filename=None, **kwargs):
        identifier = self._get_identifier(identifier)
        filename = self._get_filename(filename)
        return download_validation(identifier=identifier, filename=filename, **kwargs)

    def select(self, identifier=None, **kwargs):
        identifier = self._get_identifier(identifier)
        self.table = select_validation(identifier=identifier, **kwargs)
        return self.table


Validation = Validation()
