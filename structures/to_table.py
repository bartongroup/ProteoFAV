#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Functions that handle the reading data files and extracting their information
as a pandas.DataFrame. Also include wrapper functions that select and index
the information. Prefers the use o the wrapper instead the private functions
for better error handling. Both levels are convered by test cases.
"""

import logging
from StringIO import StringIO
from os import path

import pandas as pd
from lxml import etree
from scipy.spatial import cKDTree

from config import defaults
from utils import fetch_files
from utils import get_url_or_retry
from utils import is_valid

log = logging.getLogger(__name__)

__all__ = ["select_cif", "select_sifts", "select_dssp", "select_validation",
           "sifts_best", "_rcsb_description"]
UNIFIED_COL = ['auth_seq_id', 'pdbx_PDB_model_num', 'auth_asym_id']

##############################################################################
# Private methods
##############################################################################
def _to_unique(series):
    """Lambda-like expression for returning unique elements of a Series.
    :param series: pandas.Series
    :return: pandas.Series
    """
    return series.unique()


def _dssp(filename):
    """Parses DSSP file output to a pandas DataFrame.

    :param filename: input SIFTS xml file
    :return: pandas table dataframe
    """
    # column width descriptors
    cols_widths = ((0, 5), (6, 11), (11, 12), (13, 14), (16, 17), (35, 38),
                   (103, 109), (109, 115))
    # simplified headers for the table
    dssp_header = ("dssp_index", "icode", "chain_id", "aa", "ss", "acc", "phi",
                   "psi")
    dssp_table = pd.read_fwf(filename, skiprows=28, names=dssp_header,
                             colspecs=cols_widths, index_col=0,
                             compression=None)
    if dssp_table.empty:
        log.error('DSSP file {} resulted in a empty Dataframe'.format(filename))
        raise ValueError('DSSP file {} resulted in a empty Dataframe'.format(
                filename))
    return dssp_table


def _mmcif_atom(filename, delimiter=None):
    """Parse mmCIF ATOM and HETEROATOM lines to a pandas DataFrame.

    :param filename: input CIF file path
    :return: pandas table dataframe
    """

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
                             delim_whitespace=True,
                             low_memory=False,
                             names=_header_mmcif,
                             compression=None)
    else:
        return pd.read_table(StringIO(lines),
                             sep=str(delimiter),
                             low_memory=False,
                             names=_header_mmcif,
                             compression=None)


def _sifts_residues(filename, cols=None):
    """Parses the residue fields of a SIFTS XML file to a pandas DataFrame

    :param filename: input SIFTS xml file
    :return: pandas table dataframe
    """

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
        col if not col.startswith("PDBe")  # this should come from reference
        else col.replace(reference, "REF")
        for col in data.columns]
    return data


def _sifts_regions(filename):
    """Parse ther region field of the SIFTS XML file to a pandas dataframe

    :param filename: input SIFTS xml file path
    :return: pandas table dataframe
    """

    if not path.isfile(filename):
        raise IOError('File {} not found or unavailable.'.format(filename))

    tree = etree.parse(filename)
    root = tree.getroot()
    namespace = 'http://www.ebi.ac.uk/pdbe/docs/sifts/eFamily.xsd'
    namespace_map = {'ns': namespace}
    db_reference = "{{{}}}db".format(namespace)
    db_detail = "{{{}}}dbDetail".format(namespace)
    rows = []
    regions = {}

    for segment in root.find('.//ns:entity[@type="protein"]',
                             namespaces=namespace_map):
        for region in segment.find('.//ns:listMapRegion',
                                   namespaces=namespace_map):
            # get region annotations
            region_annotation = {}

            # parse extra annotations for each region
            for annotation in region:
                for k, v in annotation.attrib.items():
                    # db entries
                    if annotation.tag == db_reference:
                        # skipping dbSource
                        if k == 'dbSource':
                            continue

                        start = region.attrib['start']
                        end = region.attrib['end']
                        coord = annotation.attrib.get('dbCoordSys', '')
                        region_annotation['Start'] = start
                        region_annotation['End'] = end

                        # region id
                        r = (start, end, coord)

                        # renaming all keys with dbSource prefix
                        k = "{}_{}".format(annotation.attrib["dbSource"], k)

                    # dbDetail entries
                    elif annotation.tag == db_detail:
                        # joining dbSource and property keys
                        k = "_".join([annotation.attrib["dbSource"],
                                      annotation.attrib["property"]])
                        # value is the text field in the XML
                        v = annotation.text

                    # adding to the dictionary
                    try:
                        if v in region_annotation[k]:
                            continue
                        region_annotation[k].append(v)
                    except KeyError:
                        region_annotation[k] = v
                    except AttributeError:
                        region_annotation[k] = [region_annotation[k]]
                        region_annotation[k].append(v)

                    if r not in regions:
                        regions[r] = [region_annotation]
                    else:
                        regions[r].append(region_annotation)

        # group regions together
        for region in regions:
            region_annotation = {}
            for region_annot in regions[region]:
                for k in region_annot:
                    v = region_annot[k]
                    try:
                        if v in region_annotation[k]:
                            continue
                        region_annotation[k].append(v)
                    except KeyError:
                        region_annotation[k] = v
                    except AttributeError:
                        region_annotation[k] = [region_annotation[k]]
                        region_annotation[k].append(v)

            rows.append(region_annotation)
    return pd.DataFrame(rows)


def _pdb_uniprot_sifts_mapping(identifier):
    """Queries the PDBe API for SIFTS mapping between PDB - UniProt. One to many
     relationship expected.

    :param identifier: PDB id
    :return: pandas table dataframe
    """

    if not is_valid(identifier, 'pdbe'):
        raise ValueError(
                "{} is not a valid PDB identifier.".format(identifier))

    sifts_endpoint = "mappings/uniprot/"
    url = defaults.api_pdbe + sifts_endpoint + identifier
    information = get_url_or_retry(url, json=True)

    rows = []
    for uniprot in information[identifier]['UniProt']:
        uniprots = {'uniprot_id': uniprot}
        rows.append(uniprots)
    return pd.DataFrame(rows)


def _uniprot_pdb_sifts_mapping(identifier):
    """Queries the PDBe API for SIFTS mapping between UniProt - PDB entries.
    One to many relationship expected.

    :param identifier: UniProt ID
    :return: pandas table dataframe
    """
    sifts_endpoint = "mappings/best_structures/"
    url = defaults.api_pdbe + sifts_endpoint + str(identifier)
    information = get_url_or_retry(url, json=True)

    rows = []
    for entry in information[identifier]:
        rows.append(entry)
    return pd.DataFrame(rows)


def _pdb_validation_to_table(filename, global_parameters=False):
    """Parse the PDB's validation validation file to a pandas DataFrame. Private
     method, prefer its higher level wrapper.

    :type global_parameters: bool
    :param filename: path to file
    :return: table with validation information
    :rtype: pandas.DataFrame
    """
    if not path.isfile(filename):
        raise IOError('File {} not found or unavailable.'.format(filename))

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


##############################################################################
# Public methods
##############################################################################
def table_selector(table, column, value):
    """Generic table selector

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


def select_cif(pdb_id, models='first', chains=None, lines='ATOM', atoms='CA'):
    """Produce table read from mmCIF file.

    :param atoms: Which atom should represent the structure
    :param pdb_id: PDB identifier
    :param pdb_id: PDB identifier
    :param models: protein structure entity
    :param chains: protein structure chain
    :param lines: choice of ATOM, HETEROATOMS or both.
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
            cif_table = table_selector(cif_table, 'pdbx_PDB_model_num', models)
        except AttributeError:
            err = 'Structure {} has only one model, which was kept'.format
            log.info(err(pdb_id))

    # select chains
    if chains:
        cif_table = table_selector(cif_table, 'auth_asym_id', chains)

    # select lines
    if lines:
        cif_table = table_selector(cif_table, 'group_PDB', lines)

    # select which atom line will represent
    if atoms == 'centroid':
        cif_table = residues_as_cetroid(cif_table)
    elif atoms == 'backbone_centroid':
        cif_table = table_selector(
                cif_table, 'label_atom_id', ('CA', 'N', 'C', 'O'))
        cif_table = residues_as_cetroid(cif_table)
    elif atoms:
        cif_table = table_selector(cif_table, 'label_atom_id', atoms)

    # majority case
    if not cif_table[UNIFIED_COL].duplicated().any():
        return cif_table.set_index(['auth_seq_id'])

    # from here we treat the corner cases for duplicated ids
    # in case of alternative id, use the ones with maximum occupancy
    elif len(cif_table.label_alt_id.unique()) > 1:
        idx = cif_table.groupby(['auth_seq_id']).occupancy.idxmax()
        cif_table = cif_table.ix[idx]
        return cif_table.set_index(['auth_seq_id'])

    # in case of insertion code, add it to the index
    elif len(cif_table.pdbx_PDB_ins_code.unique()) > 1:
        cif_table.pdbx_PDB_ins_code.replace("?", "", inplace=True)
        cif_table.auth_seq_id = (cif_table.auth_seq_id.astype(str) +
                                 cif_table.pdbx_PDB_ins_code)
        return cif_table.set_index(['auth_seq_id'])

    # otherwise try using the pdbe_label_seq_id
    elif not cif_table['pdbe_label_seq_id'].duplicated().any():
        return cif_table.set_index(['pdbe_label_seq_id'])
    # TODO verify for alternative id even with multiple chains/models
    log.error('Failed to find unique index for {}'.format(cif_path))
    return cif_table.set_index(['auth_seq_id'])


def residues_as_cetroid(table):
    columns_to_agg = {col: "first" if table[col].dtype == 'object' else 'mean'
                      for col in table.columns
                      if col not in UNIFIED_COL}
    columns_to_agg['auth_atom_id'] = 'unique'
    return table.groupby(by=UNIFIED_COL, as_index=False).agg(columns_to_agg)


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
        dssp_path = fetch_files(pdb_id, sources='dssp',
                                directory=defaults.db_dssp)[0]
        dssp_table = _dssp(dssp_path)
    except StopIteration:
        raise IOError('{} is unreadable.'.format(dssp_path))

    if chains:
        dssp_table = table_selector(dssp_table, 'chain_id', chains)

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
        # stardatise column types
        for col in sift_table:
            #  bool columns
            if col.startswith('is'):
                sift_table[col].fillna(False)
    if chains is None:
        return sift_table
    else:
        return table_selector(sift_table, 'PDB_dbChainId', chains)


def select_validation(pdb_id, chains=None):
    """
    Produces table from PDB validation XML file.

    :param pdb_id: PDB identifier
    :param chains: PDB protein chain
    :return: pandas dataframe
    """
    val_path = path.join(defaults.db_pdb,
                         pdb_id + defaults.validation_extension)
    try:
        val_table = _pdb_validation_to_table(val_path)
    except IOError:
        val_path = fetch_files(pdb_id, sources='validation',
                               directory=defaults.db_pdb)[0]
        val_table = _pdb_validation_to_table(val_path)

    if chains:
        val_table = table_selector(val_table, 'chain', chains)
    if not val_table.empty:
        val_table.columns = ["val_" + name for name in val_table.columns]
        return val_table
    raise ValueError('Parsing {} resulted in a empty table'.format(val_path))


def sifts_best(identifier, first=False):
    """
    Retrieves the best structures from the SIFTS endpoint in the PDBe api.

    :param identifier: Uniprot ID
    :param first: gets the first entry
    :return: url content or url content in json data structure.
    """

    sifts_endpoint = "mappings/best_structures/"
    url = defaults.api_pdbe + sifts_endpoint + str(identifier)
    response = get_url_or_retry(url, json=True)
    return response if not first else response[identifier][0]


def _rcsb_description(pdb_id, tag, key):
    api = defaults.api_rcsb
    endpoint = 'describeMol'
    query = '?structureId=' + pdb_id

    url = api + endpoint + query

    tree = etree.fromstring(get_url_or_retry(url))
    values = []
    for i in tree.iter(tag):
        values.append(i.attrib[key])

    return values


if __name__ == '__main__':
    X = select_cif('2pah', atoms='backbone_centroid')
    X = select_cif('2pah', chains='A', atoms='CA')
    C = _get_contacts_from_table(X)


