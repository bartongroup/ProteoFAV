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
from requests import HTTPError
from scipy.spatial import cKDTree

from proteofav.config import defaults
from proteofav.library import scop_3to1
from proteofav.utils import fetch_files, fetch_from_url_or_retry, get_preferred_assembly_id

log = logging.getLogger('proteofav.config')
__all__ = ['_dssp', '_mmcif_atom', '_sifts_residues_regions', '_pdb_validation_to_table',
           '_rcsb_description', '_get_contacts_from_table', '_table_selector',
           # '_residues_as_centroid', '_import_dssp_chains_ids',
           'select_cif', 'select_dssp', 'select_sifts',  'select_validation', 'sifts_best']

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
    dssp_table = pd.read_fwf(filename,
                             skiprows=28,
                             names=dssp_header,
                             colspecs=cols_widths,
                             index_col=0,
                             compression=None)
    if dssp_table.empty:
        log.error('DSSP file {} resulted in a empty Dataframe'.format(filename))
        raise ValueError('DSSP file {} resulted in a empty Dataframe'.format(
            filename))
    return dssp_table


def yield_lines(filename):
    """
    Custom function for iterating over line from filename.
    :param filename: path to filename
    :return None:
    """
    with open(filename) as lines:
        for line in lines:
            yield line


def _mmcif_atom(filename, delimiter=None):
    """
    Parse mmCIF ATOM and HETEROATOM lines to a pandas DataFrame.

    :param filename: input CIF file path
    :return: pandas table dataframe
    """
    # parsing atom lines
    _header_mmcif = []
    lines = []
    for line in yield_lines(filename):
        if line.startswith("_atom_site."):
            _header_mmcif.append(line.split('.')[1].rstrip())
        elif line.startswith("ATOM") or "ATOM" in line[0:6]:
            lines.append(line)
        elif line.startswith("HETATM"):
            lines.append(line)
    lines = "".join(lines)

    if delimiter is None:
        table = pd.read_table(StringIO(lines),
                              delim_whitespace=True,
                              low_memory=False,
                              names=_header_mmcif,
                              compression=None)
    else:
        table = pd.read_table(StringIO(lines),
                              sep=str(delimiter),
                              low_memory=False,
                              names=_header_mmcif,
                              compression=None)

    # drop the 'esd' entries
    excluded = [x for x in table.columns.values if x.endswith('_esd')]
    table = table.drop(excluded, axis=1)
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


def _sifts_residues_regions(filename, cols=None,
                            sources=('CATH', 'SCOP', 'Pfam', 'InterPro')):
    """
    Parses the residue fields of a SIFTS XML file to a pandas DataFrame.

    :param filename: input SIFTS xml file
    :param cols: option to select a columns (post-parsing)
    :param sources: option to select SIFTS dbSources
    :return: pandas table dataframe
    """

    # parsing sifts segments
    tree = etree.parse(filename)
    root = tree.getroot()
    namespace = root.nsmap[None]
    namespace_map = {'ns': namespace}
    cross_reference = "{{{}}}crossRefDb".format(namespace)
    residue_detail = "{{{}}}residueDetail".format(namespace)
    db_reference = "{{{}}}db".format(namespace)
    # db_detail = "{{{}}}dbDetail".format(namespace)
    rows = []

    sources += ('PDB', 'UniProt')

    for segment_list in root.iterfind('.//ns:entity[@type="protein"]',
                                      namespaces=namespace_map):

        entity_id = segment_list.attrib['entityId']
        for segment in segment_list:
            # 1st parse the regions found for this segment
            regions_full = {}
            for region_list in segment.iterfind('.//ns:listMapRegion',
                                                namespaces=namespace_map):
                region_source = {}
                for region in region_list:
                    # get region annotations
                    # parse extra annotations for each region
                    for annotation in region:
                        for k, v in annotation.attrib.items():
                            # db entries
                            if annotation.tag == db_reference:
                                if k == 'dbSource' and v in sources:
                                    source = v
                                else:
                                    continue

                                if source not in region_source:
                                    region_source[source] = []
                                try:
                                    coord = annotation.attrib['dbCoordSys']
                                except KeyError:
                                    coord = '-'
                                region_source[source].append([annotation.attrib['dbAccessionId'],
                                                              region.attrib['start'],
                                                              region.attrib['end'],
                                                              coord])
                regions_full[entity_id] = region_source

            # 2nd parse each residue and
            for list_residue in segment.iterfind('.//ns:listResidue',
                                                 namespaces=namespace_map):
                for residue in list_residue:
                    # get residue annotations
                    residue_annotation = {}
                    # key, value pairs
                    for k, v in residue.attrib.items():
                        # skipping dbSource
                        if k == 'dbSource' or k == 'dbCoordSys' or k == 'dbResName':
                            continue

                        k = "PDB_index"
                        # adding to the dictionary
                        residue_annotation[k] = v
                        resnum = int(v)

                    # parse extra annotations for each residue
                    for annotation in residue:
                        for k, v in annotation.attrib.items():
                            # crossRefDb entries
                            if annotation.tag == cross_reference:
                                if annotation.attrib["dbSource"] in sources:
                                    # skipping dbSource
                                    if k == 'dbSource' or k == 'dbCoordSys':
                                        continue
                                    if (annotation.attrib["dbSource"] != "PDB" and
                                            annotation.attrib["dbSource"] != "UniProt"):
                                        if k == 'dbResName' or k == 'dbResNum' or k == 'dbChainId':
                                            continue
                                    if annotation.attrib["dbSource"] == "PDB" and k == "dbAccessionId":
                                        continue

                                    # adding a new column with the regionId from the 'regions'
                                    if k == "dbAccessionId":
                                        if annotation.attrib["dbSource"] in regions_full[entity_id]:
                                            for c, entry in enumerate(regions_full[entity_id][annotation.attrib["dbSource"]]):
                                                if v == entry[0]:
                                                    start = int(entry[1])
                                                    end = int(entry[2])
                                                    if resnum in range(start, end + 1, 1):
                                                        nk = "{}_regionId".format(annotation.attrib["dbSource"])
                                                        nv = str(c + 1)
                                                        residue_annotation[nk] = nv

                                    # renaming all keys with dbSource prefix
                                    k = "{}_{}".format(
                                            annotation.attrib["dbSource"], k)

                            if annotation.tag == residue_detail:
                                k = "PDB_{}".format(annotation.attrib["property"])
                                # value is the text field in the XML
                                v = annotation.text

                            # adding to the dictionary
                            if "_" in k:
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
            cif_table = _mmcif_atom(cif_path)
        except IOError:
            cif_path = fetch_files(pdb_id + '-assembly-' + assembly_id,
                                   sources='bio', directory=defaults.db_mmcif)[0]
            cif_table = _mmcif_atom(cif_path)
    else:
        # load the table
        cif_path = path.join(defaults.db_mmcif, pdb_id + '.cif')
        try:
            cif_table = _mmcif_atom(cif_path)
        except IOError:
            cif_path = fetch_files(pdb_id, sources='cif', directory=defaults.db_mmcif)[0]
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
            # TODO:
            # Could not find the correct PDB chain. It happens for protein structures with complex
            # chain identifier, as 4v9d.
            # dssp_table = _import_dssp_chains_ids(pdb_id)
            # dssp_table = _table_selector(dssp_table, 'chain_id', chains)
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
    sifts_path = path.join(defaults.db_sifts, pdb_id + '.xml')

    try:
        sift_table = _sifts_residues_regions(sifts_path,
                                             sources=('CATH', 'SCOP', 'Pfam'))
    except IOError:
        sifts_path = fetch_files(pdb_id, sources='sifts',
                                 directory=defaults.db_sifts)[0]
        sift_table = _sifts_residues_regions(sifts_path)
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
    val_path = path.join(defaults.db_validation, pdb_id + defaults.validation_extension)
    try:
        val_table = _pdb_validation_to_table(val_path)
    except IOError:
        val_path = fetch_files(pdb_id, sources='validation', directory=defaults.db_pdb)[0]
        val_table = _pdb_validation_to_table(val_path)

    if chains:
        val_table = _table_selector(val_table, 'chain', chains)
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
    url = defaults.api_pdbe + sifts_endpoint + str(uniprot_id)
    try:
        response = fetch_from_url_or_retry(url, json=True).json()
    except HTTPError as e:
        if e.response.status_code == 404:
            logging.error('No SIFTS mapping found for {}'.format(uniprot_id))
            return None
        else:
            raise
    return response if not first else response[uniprot_id][0]


if __name__ == '__main__':
    pass
