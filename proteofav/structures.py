#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Functions that handle the reading data files and extracting their information
as a pandas.DataFrame. Also include wrapper functions that select and index
the information. Prefers the use o the wrapper instead the private functions
for better error handling. Both levels are covered by test cases.
"""
from __future__ import absolute_import
import os
import logging
import pandas as pd
from lxml import etree
from io import StringIO
from requests import HTTPError
from scipy.spatial import cKDTree
from proteofav.parsers import parse_mmcif_atoms_from_file
from proteofav.parsers import parse_dssp_from_file
from proteofav.parsers import parse_sifts_residues_from_file

from proteofav.utils import fetch_files
from proteofav.utils import get_url_or_retry
from proteofav.utils import get_preferred_assembly_id
from proteofav.utils import row_selector

from proteofav.config import defaults as config

log = logging.getLogger('proteofav.config')
__all__ = ['_pdb_validation_to_table',
           '_rcsb_description', '_get_contacts_from_table',
           # '_residues_as_centroid', '_import_dssp_chains_ids',
           'select_cif', 'select_dssp', 'select_sifts', 'select_validation', 'sifts_best']

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
    api = config.api_rcsb
    endpoint = 'describeMol'
    query = '?structureId=' + pdb_id

    url = api + endpoint + query

    tree = etree.fromstring(get_url_or_retry(url))
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


def _import_dssp_chains_ids(pdb_id):
    """Imports mmCIF chain identifier to DSSP.

    :param pdb_id:
    :return: DSSP table with corrected chain ids.
    """
    dssp_table = select_dssp(pdb_id)
    cif_table = select_cif(pdb_id)
    cif_seq = cif_table.auth_comp_id.apply(scop_3to1.get)
    dssp_has_seq = dssp_table.aa.isin(scop_3to1.values())
    dssp_seq = dssp_table.aa[dssp_has_seq]
    # Import only if the sequences are identical
    if not (cif_seq == dssp_seq).all():
        err = ('Inconsitent DSSP / mmCIF sequence for {} protein structure cannot be fixed'
               'by import_dssp_chains_ids')
        raise ValueError(err.format(pdb_id))
    dssp_table.loc[dssp_has_seq, 'chain_id'] = cif_table.auth_asym_id
    return dssp_table


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
        cif_path = os.path.join(config.db_mmcif,
                             pdb_id + '-assembly-' + assembly_id + '.cif')

        try:
            cif_table = parse_mmcif_atoms_from_file(cif_path)
        except IOError:
            cif_path = fetch_files(pdb_id + '-assembly-' + assembly_id,
                                   sources='bio', directory=config.db_mmcif)[0]
            cif_table = parse_mmcif_atoms_from_file(cif_path)
    else:
        # load the table
        cif_path = os.pathjoin(config.db_mmcif, pdb_id + '.cif')
        try:
            cif_table = parse_mmcif_atoms_from_file(cif_path)
        except IOError:
            cif_path = fetch_files(pdb_id, sources='cif', directory=config.db_mmcif)[0]
            cif_table = parse_mmcif_atoms_from_file(cif_path)

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


def select_dssp(pdb_id, chains=None):
    """
    Produce table from DSSP file output.

    :param pdb_id: PDB identifier
    :param chains: PDB protein chain
    :return: pandas dataframe
    """
    dssp_path = os.path.join(config.db_dssp, pdb_id + '.dssp')
    try:
        dssp_table = parse_dssp_from_file(dssp_path)
    except IOError:
        dssp_path = fetch_files(pdb_id, sources='dssp', directory=config.db_dssp)[0]
        dssp_table = parse_dssp_from_file(dssp_path)
    except StopIteration:
        raise IOError('{} is unreadable.'.format(dssp_path))

    if chains:
        try:
            dssp_table = row_selector(dssp_table, 'chain_id', chains)
        except ValueError:
            # TODO:
            # Could not find the correct PDB chain. It happens for protein structures with complex
            # chain identifier, as 4v9d.
            # dssp_table = _import_dssp_chains_ids(pdb_id)
            # dssp_table = row_selector(dssp_table, 'chain_id', chains)
            log.error('Error loading DSSP file: Chain {} not in {}'.format(chains, pdb_id))
            return None
    # remove dssp line of transition between chains
    dssp_table = dssp_table[dssp_table.aa != '!']

    dssp_table.reset_index(inplace=True)
    if dssp_table.duplicated(['icode', 'chain_id']).any():
        log.info('DSSP file for {} has not unique index'.format(pdb_id))
    try:
        dssp_table.loc[:, 'icode'] = dssp_table.loc[:, 'icode'].astype(int)
    except ValueError:
        log.warning("{} insertion code detected in the DSSP file.".format(pdb_id))
        dssp_table.loc[:, 'icode'] = dssp_table.loc[:, 'icode'].astype(str)

    return dssp_table


def select_sifts(pdb_id, chains=None):
    """
    Produce table ready from SIFTS XML file.

    :param pdb_id: PDB identifier
    :param chains: Protein structure chain
    :return: table read to be merged
    """
    sifts_path = os.path.join(config.db_sifts, pdb_id + '.xml')

    try:
        sift_table = parse_sifts_residues_from_file(sifts_path)
    except IOError:
        sifts_path = fetch_files(pdb_id, sources='sifts',
                                 directory=config.db_sifts)[0]
        sift_table = parse_sifts_residues_from_file(sifts_path)
        # standardise column types
    for col in sift_table:
        #  bool columns
        if col.startswith('is'):
            sift_table[col].fillna(False)
    if chains is None:
        return sift_table
    else:
        return row_selector(sift_table, 'PDB_dbChainId', chains)


def select_validation(pdb_id, chains=None):
    """
    Produces table from PDB validation XML file.

    :param pdb_id: PDB identifier
    :param chains: PDB protein chain
    :return: pandas dataframe
    """
    val_path = os.path.join(config.db_validation, pdb_id + config.validation_extension)
    try:
        val_table = _pdb_validation_to_table(val_path)
    except IOError:
        val_path = fetch_files(pdb_id, sources='validation', directory=config.db_pdb)[0]
        val_table = _pdb_validation_to_table(val_path)

    if chains:
        val_table = row_selector(val_table, 'chain', chains)
    if not val_table.empty:
        val_table.columns = ["validation_" + name for name in val_table.columns]
        return val_table
    raise ValueError('Parsing {} resulted in a empty table'.format(val_path))


def sifts_best(uniprot_id, first=False):
    """
    Retrieves the best structures from the SIFTS endpoint in the PDBe api.

    :param uniprot_id: Uniprot ID
    :param first: gets the first entry
    :return: url content or url content in json data structure.
    """
    sifts_endpoint = "mappings/best_structures/"
    url = config.api_pdbe + sifts_endpoint + str(uniprot_id)
    try:
        response = get_url_or_retry(url, json=True)
    except HTTPError as e:
        if e.response.status_code == 404:
            logging.error('No SIFTS mapping found for {}'.format(uniprot_id))
            return None
        else:
            raise
    return response if not first else response[uniprot_id][0]


if __name__ == '__main__':
    pass
