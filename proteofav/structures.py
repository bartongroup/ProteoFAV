#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Functions that handle the reading data files and extracting their information
as a pandas.DataFrame. Also include wrapper functions that select and index
the information. Prefers the use o the wrapper instead the private functions
for better error handling. Both levels are covered by test cases.
"""

import logging
import pandas as pd
from os import path
from lxml import etree
from StringIO import StringIO
from requests import HTTPError
from scipy.spatial import cKDTree


from .config import defaults
from .utils import fetch_files, get_url_or_retry
from .library import scop_3to1

log = logging.getLogger(__name__)

UNIFIED_COL = ['pdbx_PDB_model_num', 'auth_asym_id', 'auth_seq_id']


##############################################################################
# Private methods
##############################################################################
def _dssp(filename):
    """
    Parses DSSP file output to a pandas DataFrame.

    :param filename: input SIFTS xml file
    :return: pandas table dataframe
    """
    # column width descriptors
    cols_widths = ((0, 5), (6, 11), (11, 12), (13, 14), (16, 17), (35, 38),
                   (103, 109), (109, 115))
    # simplified headers for the table
    dssp_header = ("dssp_index", "icode", "chain_id", "aa", "ss", "acc", "phi", "psi")
    dssp_table = pd.read_fwf(filename, skiprows=28, names=dssp_header,
                             colspecs=cols_widths, index_col=0, compression=None)
    if dssp_table.empty:
        log.error('DSSP file {} resulted in a empty Dataframe'.format(filename))
        raise ValueError('DSSP file {} resulted in a empty Dataframe'.format(
                filename))
    return dssp_table


def _mmcif_atom(filename, delimiter=None):
    """
    Parse mmCIF ATOM and HETEROATOM lines to a pandas DataFrame.

    :param filename: input CIF file path
    :return: pandas table dataframe
    """
    # parsing atom lines
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
        return pd.read_table(StringIO(lines),
                             delim_whitespace=True, low_memory=False,
                             names=_header_mmcif, compression=None)
    else:
        return pd.read_table(StringIO(lines),
                             sep=str(delimiter), low_memory=False,
                             names=_header_mmcif, compression=None)


def _sifts_residues(filename, cols=None):
    """
    Parses the residue fields of a SIFTS XML file to a pandas DataFrame.

    :param filename: input SIFTS xml file
    :return: pandas table dataframe
    """
    # parsing sifts segments
    tree = etree.parse(filename)
    root = tree.getroot()
    namespace = root.nsmap[None]
    nsmap = {'ns': namespace}
    cross_reference = "{{{}}}crossRefDb".format(namespace)
    residue_detail = "{{{}}}residueDetail".format(namespace)
    rows = []
    reference = root.attrib['dbCoordSys']

    for segment in root.iterfind('.//ns:entity[@type="protein"]',
                                 namespaces=nsmap):
        for list_residue in segment.iterfind('.//ns:listResidue',
                                             namespaces=nsmap):
            for residue in list_residue:
                # get residue annotations
                residue_annotation = {}
                # key, value pairs
                for k, v in residue.attrib.items():
                    # skipping dbSource
                    if k == 'dbSource':
                        continue
                    # renaming all keys with dbSource prefix
                    try:
                        k = "{}_{}".format(residue.attrib["dbSource"], k)
                    except KeyError:
                        k = "{}_{}".format("REF", k)
                    # adding to the dictionary
                    residue_annotation[k] = v
                # parse extra annotations for each residue
                for annotation in residue:
                    for k, v in annotation.attrib.items():
                        # crossRefDb entries
                        if annotation.tag == cross_reference:
                            # skipping dbSource
                            if k == 'dbSource':
                                continue

                            # renaming all keys with dbSource prefix
                            try:
                                k = "{}_{}".format(
                                        annotation.attrib["dbSource"], k)
                            except KeyError:
                                k = "{}_{}".format("REF", k)

                        # residueDetail entries
                        elif annotation.tag == residue_detail:
                            if annotation.attrib["property"] == 'Annotation':
                                k = 'is ' + annotation.text.lower()
                                k = k.replace(' ', '_')
                                v = True
                            else:
                                k = "_".join([annotation.attrib["dbSource"],
                                              annotation.attrib["property"]])
                                # value is the text field in the XML
                                v = annotation.text

                        # adding to the dictionary
                        try:
                            if v in residue_annotation[k]:
                                continue
                            residue_annotation[k].append(v)
                        except KeyError:
                            residue_annotation[k] = v
                        except AttributeError:
                            residue_annotation[k] = [residue_annotation[k]]
                            residue_annotation[k].append(v)
                        except TypeError:
                            # bool column for annotation
                            residue_annotation[k] = v

                rows.append(residue_annotation)
    if cols:
        data = pd.DataFrame(rows, columns=cols)
    else:
        data = pd.DataFrame(rows)
    data.columns = [
        # this should come from reference
        col if not col.startswith("PDBe")
        else col.replace(reference, "REF")
        for col in data.columns]
    return data


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

    df['contacts'] = nearby
    return df


def _table_selector(table, column, value):
    """
    Generic table selector.

    :param table: a pandas DataFrame
    :param column: a column in the DataFrame
    :param value: the query
    :return: the DataFrame filtered for
    """
    if value == 'first':
        value = table[column].iloc[0]
        table = table[table[column] == value]

    elif not hasattr(value, '__iter__') or isinstance(value, str):
        table = table[table[column] == value]

    else:
        table = table[table[column].isin(value)]

    if table.empty:
        raise ValueError('Column {} does not contain {} value(s)'.format(
                column, value))
    return table


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
def select_cif(pdb_id, models='first', chains=None, lines='ATOM', atoms='CA'):
    """
    Produce table read from mmCIF file.

    :param atoms: Which atom should represent the structure
    :param pdb_id: PDB identifier
    :param pdb_id: PDB identifier
    :param models: protein structure entity
    :param chains: protein structure chain
    :param lines: choice of ATOM, HETATMS or both (list).
    :return: Table read to be joined
    """
    # load the table
    cif_path = path.join(defaults.db_mmcif, pdb_id + '.cif')
    try:
        cif_table = _mmcif_atom(cif_path)
    except IOError:
        cif_path = fetch_files(pdb_id, sources='cif',
                               directory=defaults.db_mmcif)[0]
        cif_table = _mmcif_atom(cif_path)

    # select the models
    if models:
        try:
            cif_table = _table_selector(cif_table, 'pdbx_PDB_model_num', models)
        except AttributeError:
            err = 'Structure {} has only one model, which was kept'.format
            log.info(err(pdb_id))

    # select chains
    if chains:
        cif_table = _table_selector(cif_table, 'auth_asym_id', chains)

    # select lines
    if lines:
        cif_table = _table_selector(cif_table, 'group_PDB', lines)

    # select which atom line will represent
    if atoms == 'centroid':
        cif_table = _residues_as_centroid(cif_table)
    elif atoms == 'backbone_centroid':
        cif_table = _table_selector(
                cif_table, 'label_atom_id', ('CA', 'N', 'C', 'O'))
        cif_table = _residues_as_centroid(cif_table)
    elif atoms:
        cif_table = _table_selector(cif_table, 'label_atom_id', atoms)

    # from here we treat the corner cases for duplicated ids
    # in case of alternative id, use the ones with maximum occupancy
    if len(cif_table.label_alt_id.unique()) > 1:
        idx = cif_table.groupby(['auth_seq_id']).occupancy.idxmax()
        cif_table = cif_table.ix[idx]

    # in case of insertion code, add it to the index
    if len(cif_table.pdbx_PDB_ins_code.unique()) > 1:
        cif_table['pdbx_PDB_ins_code'].replace("?", "", inplace=True)
        cif_table['auth_seq_id'] = (cif_table['auth_seq_id'].astype(str) +
                                    cif_table['pdbx_PDB_ins_code'])

    # otherwise try using the pdbe_label_seq_id
    if 'pdbe_label_seq_id' in cif_table:
        if cif_table['pdbe_label_seq_id'].duplicated().any():
            cif_table['auth_seq_id'] = cif_table['pdbe_label_seq_id']

    # id is the atom identifier and it is need for all atoms tables.
    if cif_table[UNIFIED_COL + ['id']].duplicated().any():
        log.error('Failed to find unique index for {}'.format(cif_path))
    return cif_table.set_index(['auth_seq_id'])


def select_dssp(pdb_id, chains=None):
    """
    Produce table from DSSP file output.

    :param pdb_id: PDB identifier
    :param chains: PDB protein chain
    :return: pandas dataframe
    """
    dssp_path = path.join(defaults.db_dssp, pdb_id + '.dssp')
    try:
        dssp_table = _dssp(dssp_path)
    except IOError:
        dssp_path = fetch_files(pdb_id, sources='dssp', directory=defaults.db_dssp)[0]
        dssp_table = _dssp(dssp_path)
    except StopIteration:
        raise IOError('{} is unreadable.'.format(dssp_path))

    if chains:
        try:
            dssp_table = _table_selector(dssp_table, 'chain_id', chains)
        except ValueError:
            # Could not find the correct PDB chain. It happens for protein structures with complex
            # chain identifier, as 4v9d.
            dssp_table = _import_dssp_chains_ids(pdb_id)
            dssp_table = _table_selector(dssp_table, 'chain_id', chains)

    if dssp_table.index.name == 'icode':
        if dssp_table.index.duplicated().any():
            log.info('DSSP file for {} has not unique index'.format(pdb_id))
        return dssp_table

    if chains and dssp_table.icode.duplicated().any():
        log.info('DSSP file for {} has not unique index'.format(pdb_id))
    return dssp_table.set_index(['icode'])


def select_sifts(pdb_id, chains=None):
    """
    Produce table ready from SIFTS XML file.

    :param pdb_id: PDB identifier
    :param chains: Protein structure chain
    :return: table read to be merged
    """
    sift_path = path.join(defaults.db_sifts, pdb_id + '.xml')

    try:
        sift_table = _sifts_residues(sift_path)
    except IOError:
        sift_path = fetch_files(pdb_id, sources='sifts',
                                directory=defaults.db_sifts)[0]
        sift_table = _sifts_residues(sift_path)
        # standardise column types
        for col in sift_table:
            #  bool columns
            if col.startswith('is'):
                sift_table[col].fillna(False)
    if chains is None:
        return sift_table
    else:
        return _table_selector(sift_table, 'PDB_dbChainId', chains)


def select_validation(pdb_id, chains=None):
    """
    Produces table from PDB validation XML file.

    :param pdb_id: PDB identifier
    :param chains: PDB protein chain
    :return: pandas dataframe
    """
    val_path = path.join(defaults.db_pdb, pdb_id + defaults.validation_extension)
    try:
        val_table = _pdb_validation_to_table(val_path)
    except IOError:
        val_path = fetch_files(pdb_id, sources='validation', directory=defaults.db_pdb)[0]
        val_table = _pdb_validation_to_table(val_path)

    if chains:
        val_table = _table_selector(val_table, 'chain', chains)
    if not val_table.empty:
        val_table.columns = ["val_" + name for name in val_table.columns]
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
    url = defaults.api_pdbe + sifts_endpoint + str(uniprot_id)
    try:
        response = get_url_or_retry(url, json=True)
    except HTTPError, e:
        if e.response.status_code == 404:
            logging.error('No SIFTS mapping found for {}'.format(uniprot_id))
            return None
        else:
            raise
    return response if not first else response[uniprot_id][0]


if __name__ == '__main__':
    X = select_dssp('4v9d', chains='BD')