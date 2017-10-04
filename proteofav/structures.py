#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Functions that handle the reading data files and extracting their information
as a pandas.DataFrame. Also include wrapper functions that select and index
the information. Prefers the use o the wrapper instead the private functions
for better error handling. Both levels are covered by test cases.
"""
from __future__ import absolute_import
import logging
try:
    # python 2.7
    from StringIO import StringIO
except ImportError:
    from io import StringIO
from os import path

import pandas as pd
from lxml import etree
from scipy.spatial import cKDTree
from string import ascii_uppercase

from proteofav.config import defaults
from proteofav.utils import (fetch_files, fetch_from_url_or_retry,
                             get_preferred_assembly_id, row_selector,
                             InputFileHandler,
                             constrain_column_types, exclude_columns)
from proteofav.library import pdbx_types

log = logging.getLogger('proteofav.config')
__all__ = ['parse_mmcif_atoms', '_pdb_validation_to_table',
           '_rcsb_description', '_get_contacts_from_table',
           # '_residues_as_centroid', '_import_dssp_chains_ids',
           'select_cif', 'select_validation']

UNIFIED_COL = ['pdbx_PDB_model_num', 'auth_asym_id', 'auth_seq_id']


##############################################################################
# Private methods
##############################################################################
def yield_lines(filename):
    """
    Custom function for iterating over line from filename.
    :param filename: path to filename
    :return None:
    """
    with open(filename) as lines:
        for line in lines:
            yield line


def parse_mmcif_atoms(filename, excluded_cols=None):
    """
    Parse mmCIF ATOM and HETATM lines.

    :param filename: path to the mmCIF file
    :param excluded_cols: list of columns to be excluded
    :return: returns a pandas DataFrame
    """

    log.info("Parsing mmCIF atoms from lines...")

    # example lines with some problems
    """
    _atom_site.pdbx_PDB_model_num
    _atom_site.pdbe_label_seq_id
    _atom_site.orig_label_asym_id
    _atom_site.orig_auth_asym_id
    ATOM 1 N N . VAL A 1 1 ? -7.069 21.943 18.770 1.0 56.51 ? ? ? ? ? ? 118 VAL A N 1 1 A A
    ATOM 2 C CA . VAL A 1 1 ? -7.077 21.688 20.244 1.0 59.09 ? ? ? ? ? ? 118 VAL A CA 1 1 A A
    ATOM 3 C C . VAL A 1 1 ? -5.756 21.077 20.700 1.0 44.63 ? ? ? ? ? ? 118 VAL A C 1 1 A A
    ATOM 4 O O . VAL A 1 1 ? -5.346 20.029 20.204 1.0 59.84 ? ? ? ? ? ? 118 VAL A O 1 1 A A
    """

    InputFileHandler(filename)

    # parsing atom lines
    header = []
    lines = []
    with open(filename) as inlines:
        for line in inlines:
            if line.startswith("_atom_site."):
                header.append(line.split('.')[1].rstrip())
            elif line.startswith("ATOM") or "ATOM" in line[0:6]:
                lines.append(line)
            elif line.startswith("HETATM"):
                lines.append(line)
    lines = "".join(lines)

    all_str = {key: str for key in header}
    table = pd.read_table(StringIO(lines), delim_whitespace=True, low_memory=False,
                          names=header, compression=None, converters=all_str,
                          keep_default_na=False)
    # excluding columns
    if excluded_cols is None:
        excluded_cols = ('Cartn_x_esd', 'Cartn_y_esd', 'Cartn_z_esd',
                         'occupancy_esd', 'B_iso_or_equiv_esd',
                         'pdbx_formal_charge')

    table = exclude_columns(table, excluded=excluded_cols)

    # enforce some specific column types
    table = constrain_column_types(table, col_type_dict=pdbx_types)

    if table.empty:
        log.error('mmCIF file {} resulted in a empty Dataframe'.format(filename))
        raise ValueError('mmCIF file {} resulted in a empty Dataframe'.format(
            filename))
    return table


def parse_pdb_atoms(filename, excluded_cols=None,
                    fix_label_alt_id=True, fix_ins_code=True, fix_type_symbol=True):
    """
    Parse PDB ATOM and HETATM lines. The ATOM lines are imported
    to the dictionary names used in the mmCIF format.

    :param filename: path to the PDB file
    :param excluded_cols: list of columns to be excluded
    :param fix_label_alt_id: boolean
    :param fix_ins_code: boolean
    :param fix_type_symbol: boolean
    :return: returns a pandas DataFrame
    """

    log.info("Parsing PDB atoms from lines...")

    # example lines
    """
    MODEL        1
    ATOM      0  N   SER A  -1     104.083  78.916  -1.349  1.00 61.47           N
    ATOM      1  N   VAL A 118      -7.069  21.943  18.770  1.00 56.51           N  
    ATOM      2  CA  VAL A 118      -7.077  21.688  20.244  1.00 59.09           C  
    ATOM      3  C   VAL A 118      -5.756  21.077  20.700  1.00 44.63           C  
    ATOM      4  O   VAL A 118      -5.346  20.029  20.204  1.00 59.84           O
    """

    InputFileHandler(filename)

    # parsing atom lines, converting it to mmcif-style headers
    lines = []
    modelnumb = '1'
    with open(filename) as inlines:
        for line in inlines:
            line = line.rstrip()
            line = line[0:78]
            if line.startswith("MODEL"):
                modelnumb = line.split()[1]
            elif line.startswith("ATOM"):
                lines.append(line + "%s" % modelnumb)
            elif line.startswith("HETATM"):
                lines.append(line + "%s" % modelnumb)
    lines = "\n".join(lines)

    header = ('group_PDB', 'id', 'label_atom_id', 'label_alt_id', 'label_comp_id',
              'label_asym_id', 'label_seq_id_full', 'label_seq_id',
              'pdbx_PDB_ins_code', 'Cartn_x', 'Cartn_y', 'Cartn_z',
              'occupancy', 'B_iso_or_equiv', 'type_symbol', 'auth_atom_id', 'auth_comp_id',
              'auth_asym_id', 'auth_seq_id_full', 'auth_seq_id', 'pdbx_PDB_model_num')

    # https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html
    widths = ((0, 6), (6, 11), (12, 16), (16, 17), (17, 20), (21, 22), (22, 27), (22, 26), (26, 27),
              (30, 38), (38, 46), (46, 54), (54, 60), (60, 66), (76, 78),  # (72, 76), ('seg_id')
              (12, 16), (17, 20), (21, 22), (22, 27), (22, 26), (78, 79))

    all_str = {key: str for key in header}
    table = pd.read_fwf(StringIO(lines), names=header, colspecs=widths,
                        compression=None, converters=all_str, keep_default_na=False)

    # fixes the 'pdbx_PDB_ins_code'
    if fix_ins_code:
        table = _fix_pdb_ins_code(table)
    # fixes the 'label_alt_id
    if fix_label_alt_id:
        table = _fix_label_alt_id(table)
    # fixes 'type_symbol' if missing
    if fix_type_symbol:
        table = _fix_type_symbol(table)

    # excluding columns
    if excluded_cols is None:
        excluded_cols = ('Cartn_x_esd', 'Cartn_y_esd', 'Cartn_z_esd',
                         'occupancy_esd', 'B_iso_or_equiv_esd',
                         'pdbx_formal_charge')

    table = exclude_columns(table, excluded=excluded_cols)

    # enforce some specific column types
    table = constrain_column_types(table, col_type_dict=pdbx_types)

    if table.empty:
        log.error('PDB file {} resulted in a empty Dataframe'.format(filename))
        raise ValueError('PDB file {} resulted in a empty Dataframe'.format(
            filename))
    return table


def _fix_pdb_ins_code(table):
    """
    Utility that fixes the 'pdbx_PDB_ins_code' column to match what is expected
    in the mmCIF format.

    :param table: pandas DataFrame object
    :return: returns a modified pandas DataFrame
    """

    table['pdbx_PDB_ins_code'] = table['pdbx_PDB_ins_code'].str.replace('\ |', '?')
    table['pdbx_PDB_ins_code'] = table['pdbx_PDB_ins_code'].fillna('?').astype(str)
    return table


def _fix_label_alt_id(table):
    """
    Utility that fixes the 'label_alt_id' column to match what is
    expected in the mmCIF format.

    :param table: pandas DataFrame object
    :return: returns a modified pandas DataFrame
    """

    table['label_alt_id'] = table['label_alt_id'].str.replace('\ |\?', '.')
    table['label_alt_id'] = table['label_alt_id'].fillna('.').astype(str)
    return table


def _fix_type_symbol(table):
    """
    Utility that fixes the 'type_symbol' column to match what is
    expected in the mmCIF format - when missing in the Structure.

    :param table: pandas DataFrame object
    :return: returns a modified pandas DataFrame
    """

    def get_type_symbol(table, key, key_fix):
        # this maybe a bit crude way of assigning this value
        if table[key] != " " and table[key] != "" and len(table[key]):
            return table[key]
        else:
            return ''.join([x for x in table[key_fix] if x in ascii_uppercase])[0]

    table.is_copy = False
    table['type_symbol'] = table.apply(get_type_symbol, axis=1,
                                       args=('type_symbol', 'label_atom_id'))
    return table


def _mmcif_fields(filename, field_name='exptl.',
                  require_index=False):
    """
    Generic method that gets a particular field to pandas table.
    :param filename: input mmCIF file
    :param field_name: name of the field to be parsed
    :param require_index: boolean for when lines need to start with and index ID
      of any sort
    :return: Pandas table
    """
    header = []
    lines = []
    with open(filename, "r+") as handle:
        for line in handle:
            if line.startswith(field_name):
                break
            last_line = line

        if 'loop_' in last_line:
            while line.startswith(field_name):
                header.append(line.replace(field_name, '').rstrip())
                try:
                    line = handle.next()
                except AttributeError:
                    line = next(handle)
            while not line.startswith('#'):
                lines.append(line.replace('"', "'"))
                try:
                    line = handle.next()
                except AttributeError:
                    line = next(handle)
        else:
            while line.startswith(field_name):
                line = line.replace(field_name, '').rstrip()
                head, data = line.split(None, 1)
                header.append(head)
                lines.append(data)
                try:
                    line = handle.next()
                except AttributeError:
                    line = next(handle)
            lines = (' '.join(lines)).replace('"', "'")

    if require_index:
        # requires the lines to start with and index ID
        nlines = []
        for e in lines:
            try:
                int(e[0:2])
                e = e.replace('\n', '')
            except (TypeError, ValueError):
                pass
            nlines.append(e)
        lines = ''.join(nlines)
    else:
        lines = ''.join(lines)

    table = pd.read_table(StringIO(lines),
                          names=header,
                          delim_whitespace=True,
                          quotechar="'",
                          index_col=False)
    return table


def _pdb_validation_to_table(filename, global_parameters=False):
    """
    Parse the PDB's validation validation file to a pandas DataFrame.
    Private method, prefer its higher level wrapper.

    :type global_parameters: bool
    :param filename: path to file
    :return: table with validation information
    :rtype: pandas.DataFrame
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
    df = pd.DataFrame(rows, columns=header)
    return df


def _rcsb_description(pdb_id, tag, key):
    """
    Gets description from RCSB PDB api.

    :param pdb_id: PDB id
    :param tag: name tag as defined in the api
    :param key: key name as defined in the api
    :return: list of values
    """
    api = defaults.api_rcsb
    endpoint = 'describeMol'
    query = '?structureId=' + pdb_id

    url = api + endpoint + query

    tree = etree.fromstring(fetch_from_url_or_retry(url).content)
    values = []
    for i in tree.iter(tag):
        values.append(i.attrib[key])
    return values


def _get_contacts_from_table(df, distance=5, ignore_consecutive=3):
    """
    Just a simple testing distance measure.

    :param df: pd.Dataframe
    :param distance: distance threshold in Angstrom
    :param ignore_consecutive: number of consecutive residues that will be ignored
      (in both directions)
    :return: new pd.Dataframe
    """
    ig = ignore_consecutive

    # using KDTree
    tree = cKDTree(df[['Cartn_y', 'Cartn_y', 'Cartn_z']])
    nearby = []
    for i in df.index:
        query_point = df.loc[i, ['Cartn_y', 'Cartn_y', 'Cartn_z']]
        idx = tree.query_ball_point(query_point, r=distance, p=2)
        idx = df.index[idx]
        # ignoring nearby residues (not likely to be true contacts)
        # TODO: need to assess this
        idx = [j for j in idx if (j <= i - ig or j >= i + ig)]
        nearby.append(idx)
        # for j in idx:
        #     chain_i = str(df.loc[i, 'auth_asym_id'])
        #     chain_j = str(df.loc[j, 'auth_asym_id'])
        #     if chain_i != chain_j:
        #
        #         cont = [df.loc[i, ['auth_asym_id', 'auth_seq_id',
        #                            'pdbx_PDB_ins_code', 'auth_comp_id' ]],
        #                 df.loc[j, ['auth_asym_id', 'auth_seq_id',
        #                            'pdbx_PDB_ins_code', 'auth_comp_id' ]]]
        #         contacts.append(cont)

    df['contacts'] = nearby
    return df


def _residues_as_centroid(table):
    """
    Gets the residues' atoms and their centroids (mean).

    :param table: main dataframe
    :return: compressed dataframe
    """
    columns_to_agg = {col: "first" if table[col].dtype == 'object' else 'mean'
                      for col in table.columns
                      if col not in UNIFIED_COL}
    columns_to_agg['auth_atom_id'] = 'unique'
    return table.groupby(by=UNIFIED_COL, as_index=False).agg(columns_to_agg)


##############################################################################
# Public methods
##############################################################################
def select_cif(pdb_id, models='first', chains=None, lines='ATOM', atoms='CA',
               biounit=False, assembly_id=None):
    """
    Produce table read from mmCIF file.

    :param atoms: Which atom should represent the structure
    :param pdb_id: PDB identifier
    :param models: protein structure entity
    :param chains: protein structure chain
    :param lines: choice of ATOM, HETATMS or both (list).
    :param biounit: boolean to use the preferred biounit available
    :param assembly_id: only applies when biounit is True
    :return: Table read to be joined
    """

    # asymmetric unit or biological unit?
    if biounit:
        if assembly_id is None:
            # get the preferred bio assembly id from the PDBe API
            assembly_id = get_preferred_assembly_id(pdb_id)

        # load the table
        cif_path = path.join(defaults.db_mmcif,
                             pdb_id + '-assembly-' + assembly_id + '.cif')

        try:
            cif_table = parse_mmcif_atoms(cif_path)
        except IOError:
            cif_path = fetch_files(pdb_id + '-assembly-' + assembly_id,
                                   sources='bio', directory=defaults.db_mmcif)[0]
            cif_table = parse_mmcif_atoms(cif_path)
    else:
        # load the table
        cif_path = path.join(defaults.db_mmcif, pdb_id + '.cif')
        try:
            cif_table = parse_mmcif_atoms(cif_path)
        except IOError:
            cif_path = fetch_files(pdb_id, sources='cif', directory=defaults.db_mmcif)[0]
            cif_table = parse_mmcif_atoms(cif_path)

    # select the models
    if models:
        try:
            cif_table = row_selector(cif_table, 'pdbx_PDB_model_num', models)
        except AttributeError:
            err = 'Structure {} has only one model, which was kept'.format
            log.info(err(pdb_id))

    # select chains
    if chains:
        cif_table = row_selector(cif_table, 'auth_asym_id', chains)

    # select lines
    if lines:
        cif_table = row_selector(cif_table, 'group_PDB', lines)

    # select which atom line will represent
    if atoms == 'centroid':
        cif_table = _residues_as_centroid(cif_table)
    elif atoms == 'backbone_centroid':
        cif_table = row_selector(
            cif_table, 'label_atom_id', ('CA', 'N', 'C', 'O'))
        cif_table = _residues_as_centroid(cif_table)
    elif atoms:
        cif_table = row_selector(cif_table, 'label_atom_id', atoms)

    # id is the atom identifier and it is need for all atoms tables.
    if cif_table[UNIFIED_COL + ['id']].duplicated().any():
        log.error('Failed to find unique index for {}'.format(cif_path))

    return cif_table


def select_validation(pdb_id, chains=None):
    """
    Produces table from PDB validation XML file.

    :param pdb_id: PDB identifier
    :param chains: PDB protein chain
    :return: pandas dataframe
    """
    val_path = path.join(defaults.db_validation, pdb_id + defaults.validation_extension)
    try:
        val_table = _pdb_validation_to_table(val_path)
    except IOError:
        val_path = fetch_files(pdb_id, sources='validation', directory=defaults.db_pdb)[0]
        val_table = _pdb_validation_to_table(val_path)

    if chains:
        val_table = row_selector(val_table, 'chain', chains)
    if not val_table.empty:
        val_table.columns = ["validation_" + name for name in val_table.columns]
        return val_table
    raise ValueError('Parsing {} resulted in a empty table'.format(val_path))


if __name__ == '__main__':
    pass
