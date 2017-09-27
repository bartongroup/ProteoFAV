# -*- coding: utf-8 -*-


import os
import logging
import pandas as pd
from lxml import etree
from io import StringIO
from requests import HTTPError
from scipy.spatial import cKDTree
from string import ascii_uppercase

from proteofav.downloaders import Downloader
from proteofav.fetchers import get_preferred_assembly_id
from proteofav.library import aa_default_atoms
from proteofav.library import pdbx_types
from proteofav.parsers import parse_mmcif_atoms_from_file
from proteofav.parsers import parse_pdb_atoms_from_file
from proteofav.utils import InputFileHandler
from proteofav.utils import OutputFileHandler
from proteofav.utils import constrain_column_types
from proteofav.utils import exclude_columns
from proteofav.utils import fetch_files
from proteofav.utils import get_url_or_retry
from proteofav.utils import get_preferred_assembly_id
from proteofav.utils import row_selector
from proteofav.utils import GenericInput

from proteofav.config import defaults as config

log = logging.getLogger('proteofav')
__all__ = ['_pdb_validation_to_table',
           '_rcsb_description', '_get_contacts_from_table',
           # '_residues_as_centroid', '_import_dssp_chains_ids',
           'select_cif', 'select_dssp', 'select_sifts', 'select_validation', 'sifts_best']

UNIFIED_COL = ['pdbx_PDB_model_num', 'auth_asym_id', 'auth_seq_id']

_PDB_FORMAT = "%s%5i %-4s%c%3s %c%4s%c   %8.3f%8.3f%8.3f%s%6.2f      %4s%2s%2s\n"


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


class Structure(GenericInput):
    def __init__(self, identifier=None, filename=None, table=None):
        self.identifier = identifier
        self.filename = filename
        self.table = table
        self._file_format = None

    def read(self, filename=None, excluded_cols=None,
             models='first', chains=None, res=None, res_full=None,
             comps=None, atoms=None, lines=None, category='label',
             residue_agg=False, agg_method='centroid',
             add_res_full=True, add_atom_altloc=False, reset_atom_id=True,
             remove_altloc=False, remove_hydrogens=True, remove_partial_res=False,
             pdb_fix_ins_code=True, pdb_fix_label_alt_id=True, pdb_fix_type_symbol=True):
        """
        Parse ATOM and HETATM lines from mmCIF or PDB files.

        :param filename: path to the mmCIF/PDB file
        :param excluded_cols: option to exclude mmCIF columns
        :param models: (tuple) model IDs or 'first' or None
        :param chains: (tuple) chain IDs or None
        :param res: (tuple) res IDs or None
        :param res_full: (tuple) full res IDs ('*_seq_id' + 'pdbx_PDB_ins_code') or None
        :param comps: (tuple) comp IDs or None
        :param atoms: (tuple) atom IDs or None
        :param lines: (tuple) 'ATOM' or 'HETATM' or None (both)
        :param category: data category to be used as precedence in _atom_site.*_*
            asym_id, seq_id and atom_id
        :param residue_agg: (boolean) whether to do residue aggregation
        :param agg_method: current values: 'centroid', 'first', 'mean' and 'unique'
        :param add_res_full: option to extend the table with 'res_full'
            i.e. res_number + insertion_code (e.g. '12A')
        :param add_atom_altloc: boolean (new string join)
        :param reset_atom_id: boolean
        :param remove_altloc: boolean
        :param remove_hydrogens: boolean
        :param remove_partial_res: (boolean) removes amino acids with missing atoms
        :param pdb_fix_ins_code: boolean
        :param pdb_fix_label_alt_id: boolean
        :param pdb_fix_type_symbol: boolean
        :return: returns a pandas DataFrame
        """

        filename = self._get_filename(filename)
        InputFileHandler(filename)

        if excluded_cols is None:
            excluded_cols = ('Cartn_x_esd', 'Cartn_y_esd', 'Cartn_z_esd',
                             'occupancy_esd', 'B_iso_or_equiv_esd',
                             'pdbx_formal_charge')

        if self._file_format == 'mmcif':
            table = parse_mmcif_atoms_from_file(filename)
        elif self._file_format == 'pdb':
            table = parse_pdb_atoms_from_file(filename)

            # fixes the 'pdbx_PDB_ins_code'
            if pdb_fix_ins_code:
                table = _fix_pdb_ins_code(table)
            # fixes the 'label_alt_id
            if pdb_fix_label_alt_id:
                table = _fix_label_alt_id(table)
            # fixes 'type_symbol' if missing
            if pdb_fix_type_symbol:
                table = _fix_type_symbol(table)

        # excluding columns
        table = exclude_columns(table, excluded=excluded_cols)

        # enforce some specific column types
        table = constrain_column_types(table, col_type_dict=pdbx_types)

        # selections / filtering
        # excluding rows
        # if only first model (>1 in NMR structures)
        if models is not None:
            table = row_selector(table, 'pdbx_PDB_model_num', models)
            log.info("mmCIF/PDB table filtered by pdbx_PDB_model_num...")

        if chains is not None:
            table = row_selector(table, '{}_asym_id'.format(category), chains)
            log.info("mmCIF/PDB table filtered by %s_asym_id...", category)

        if res is not None:
            table = row_selector(table, '{}_seq_id'.format(category), res)
            log.info("mmCIF/PDB table filtered by %s_seq_id...", category)

        if res_full is not None:
            table = row_selector(table, '{}_seq_id_full'.format(category), res_full)
            log.info("mmCIF/PDB table filtered by %s_seq_id_full...", category)

        if comps is not None:
            table = row_selector(table, '{}_comp_id'.format(category), comps)
            log.info("mmCIF/PDB table filtered by %s_comp_id...", category)

        if atoms is not None:
            table = row_selector(table, '{}_atom_id'.format(category), atoms)
            log.info("mmCIF/PDB table filtered by %s_atom_id...", category)

        if lines is not None:
            table = row_selector(table, 'group_PDB', lines)
            log.info("mmCIF/PDB table filtered by group_PDB...")

        # table modular extensions or selections
        if add_res_full:
            table = _add_mmcif_res_full(table)
            log.info("mmCIF/PDB added full res (res + ins_code)...")

        if add_atom_altloc:
            table = _add_mmcif_atom_altloc(table)
            log.info("mmCIF/PDB added full atom (atom + altloc)...")

        if remove_altloc:
            table = _remove_multiple_altlocs(table)
            reset_atom_id = True
            log.info("mmCIF/PDB removed altlocs...")

        if remove_hydrogens:
            table = row_selector(table, key='type_symbol', value='H', reverse=True)
            log.info("mmCIF/PDB removed existing hydrogens...")

        if remove_partial_res:
            table = _remove_partial_residues(table)
            log.info("mmCIF/PDB removed incomplete residues...")

        if reset_atom_id:
            table.reset_index(inplace=True)
            table = table.drop(['index'], axis=1)
            table['id'] = table.index + 1
            log.info("mmCIF/PDB reset atom numbers...")

        if residue_agg:
            table = residues_aggregation(table, agg_method=agg_method,
                                         category=category)

        if table.empty:
            raise ValueError('{} resulted in an empty DataFrame...'.format(filename))

        self.table = table
        return self.table

    def write(self, table=None, filename=None, overwrite=False, category='label'):
        """
        Writes 'ATOM' lines in PDB or mmCIF format.

        :param table: pandas DataFrame object
        :param filename: path to the mmCIF/PDB file
        :param overwrite: boolean
        :param category: data category to be used as precedence in _atom_site.*_*
            asym_id, seq_id and atom_id (generally either 'label' or 'auth')
        :return: (side effects) writes to a file
        """

        table = self._get_table(table)
        filename = self._get_filename(filename)

        if table is None:
            table = self.table

        OutputFileHandler(filename, overwrite=overwrite)

        if self._file_format == 'mmcif':
            write_mmcif_from_table(filename=filename, table=table,
                                   override=overwrite)
        elif self._file_format == 'pdb':
            write_pdb_from_table(filename=filename, table=table,
                                 override=overwrite, category=category)

    def download(self, identifier=None, filename=None,
                 bio_unit=False, bio_unit_preferred=False, bio_unit_id="1",
                 overwrite=False, decompress=True):
        """
        Downloads a structure from the PDBe to the filesystem.
        BioUnits only work for mmCIF. If `pdb=True` the bio_unit arguments
        are ignored.

        :param identifier: (str) PDB accession ID
        :param filename: path to the mmCIF/PDB file
        :param bio_unit: (boolean) whether to get BioUnit instead of AsymUnit
        :param bio_unit_preferred: (boolean) if True fetches the 'preferred'
            biological assembly instead
        :param bio_unit_id: biological unit identifier
        :param overwrite: boolean
        :param decompress: boolean
        :return: (side effects) writes to a file
        """

        identifier = self._get_identifier(identifier)
        filename = self._get_filename(filename)
        OutputFileHandler(filename, overwrite=overwrite)

        if self._file_format == 'mmcif':
            if bio_unit:
                # atom lines only?
                # url_endpoint = ("static/entry/download/"
                #                "{}-assembly-{}_atom_site.cif.gz".format(identifier, pref))
                if bio_unit_preferred:
                    pref = get_preferred_assembly_id(identifier=identifier)
                else:
                    pref = bio_unit_id
                url_endpoint = ("static/entry/download/"
                                "{}-assembly-{}.cif.gz".format(identifier, pref))
            else:
                # original mmCIF?
                # url_endpoint = "entry-files/download/{}.cif".format(pdbid)
                url_endpoint = "entry-files/download/{}_updated.cif".format(identifier)
        elif self._file_format == 'pdb':
            url_endpoint = "entry-files/download/pdb{}.ent".format(identifier)

        url_root = config.http_pdbe
        url = url_root + url_endpoint

        Downloader(url=url, filename=filename,
                   decompress=decompress, overwrite=overwrite)


class mmCIF(Structure):
    def __init__(self, identifier=None, filename=None, table=None):
        self.identifier = identifier
        self.filename = filename
        self.table = table
        self._file_format = 'mmcif'


class PDB(Structure):
    def __init__(self, identifier=None, filename=None, table=None):
        self.identifier = identifier
        self.filename = filename
        self.table = table
        self._file_format = 'pdb'


# TODO fix for speed
def _fix_pdb_ins_code(table):
    """
    Utility that fixes the 'pdbx_PDB_ins_code' column to match is expected
    in the mmCIF format.

    :param table: pandas DataFrame object
    :return: returns a modified pandas DataFrame
    """

    ins_codes = []
    for i in table.index:
        value = table.loc[i, "pdbx_PDB_ins_code"]
        if value == '' or value == ' ' or value == '?':
            value = '?'
        ins_codes.append(value)
    table["pdbx_PDB_ins_code"] = ins_codes
    table['pdbx_PDB_ins_code'] = table['pdbx_PDB_ins_code'].fillna('?').astype(str)
    return table


# TODO fix for speed
def _fix_label_alt_id(table):
    """
    Utility that fixes the 'label_alt_id' column to match what is
    expected in the mmCIF format.

    :param table: pandas DataFrame object
    :return: returns a modified pandas DataFrame
    """

    alt_locs = []
    for i in table.index:
        value = table.loc[i, "label_alt_id"]
        if value == '' or value == ' ' or value == '?':
            value = '.'
        alt_locs.append(value)
    table["label_alt_id"] = alt_locs
    table['label_alt_id'] = table['label_alt_id'].fillna('.').astype(str)
    return table


# TODO fix for speed
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


# TODO fix for speed
def _add_mmcif_res_full(table):
    """
    Utility that adds a new column to the table.
    Adds a new column with the 'full res' (i.e. seq_id + ins_code).

    :param table: pandas DataFrame object
    :return: returns a modified pandas DataFrame
    """

    # adds both 'label' and 'auth' entries
    if 'label_seq_id' in table:
        seqs_full = []
        for ix in table.index:
            seq = "{}{}".format(table.loc[ix, 'label_seq_id'],
                                table.loc[ix, 'pdbx_PDB_ins_code']).replace('?', '')
            seqs_full.append(seq)
        assert len(seqs_full) == len(table)
        table['label_seq_id_full'] = seqs_full
    if 'auth_seq_id' in table:
        seqs_full = []
        for ix in table.index:
            seq = "{}{}".format(table.loc[ix, 'auth_seq_id'],
                                table.loc[ix, 'pdbx_PDB_ins_code']).replace('?', '')
            seqs_full.append(seq)
        assert len(seqs_full) == len(table)
        table['auth_seq_id_full'] = seqs_full

    return table


# TODO fix for speed
def _add_mmcif_atom_altloc(table):
    """
    Utility that adds new columns to the table.
    adds: 'label_atom_altloc_id', and 'auth_atom_altloc_id' which is a
        string join between '*_atom_id' + 'label_alt_id'

    :param table: pandas DataFrame object
    :return: returns a modified pandas DataFrame
    """

    def join_atom_altloc(table, category='label'):
        atom = table['{}_atom_id'.format(category)]
        altloc = table['label_alt_id']
        if altloc == "." or altloc == '' or altloc == ' ':
            return atom
        else:
            return atom + '.' + altloc

    table.is_copy = False
    table['label_atom_altloc_id'] = table.apply(join_atom_altloc,
                                                axis=1, args=('label',))
    table['auth_atom_altloc_id'] = table.apply(join_atom_altloc,
                                               axis=1, args=('auth',))
    return table


# TODO fix for speed
def _remove_multiple_altlocs(table):
    """
    Removes alternative locations (i.e. 'rows') leaving only the first.
    Needs to find rows with alt_id != '.' and then find following rows until
    '.' appears again (expects atoms to be consequent).

    :param table: pandas DataFrame object
    :return: returns a modified pandas DataFrame
    """

    drop_ixs = []
    for ix in table.index:
        altloc = table.loc[ix, 'label_alt_id']
        if altloc != '.':
            # table.loc[ix, 'label_alt_id'] = '.'
            table.set_value(ix, 'label_alt_id', '.')
            atomid = table.loc[ix, 'label_atom_id']
            try:
                for nx in range(1, 100, 1):
                    altnx = table.loc[ix + nx, 'label_alt_id']
                    atomnx = table.loc[ix + nx, 'label_atom_id']
                    if altnx != '.' and atomnx == atomid:
                        # store indexes of the rows to be dropped
                        drop_ixs.append(ix + nx)
                    else:
                        break
            except KeyError:
                break
    return table.drop(table.index[drop_ixs])


# TODO fix for speed
def _remove_partial_residues(table, category='label'):
    """
    Removes residues that contain missing atoms. Needs to check
    which atoms are available for each residue. Also removes residues with
    same '*_seq_id' as the previous residue.

    :param table: pandas DataFrame object
    :param category: data category to be used as precedence in _atom_site.*_*
        asym_id, seq_id and atom_id
    :return: returns a modified pandas DataFrame
    """

    drop_ixs = []
    curr_ixs = []
    curr_atoms = []
    prev_res = ''
    prev_seq = ''
    next_res_for_rm = False
    table.reset_index(inplace=True)
    table = table.drop(['index'], axis=1)
    for ix in table.index:
        group = table.loc[ix, 'group_PDB']
        if group == 'ATOM':
            curr_res = table.loc[ix, '{}_comp_id'.format(category)]
            curr_seq = table.loc[ix, '{}_seq_id'.format(category)]
            if curr_res in aa_default_atoms:
                curr_atom = table.loc[ix, '{}_atom_id'.format(category)]
                if prev_res == curr_res and prev_seq == curr_seq:
                    curr_ixs.append(ix)
                    curr_atoms.append(curr_atom)
                else:
                    if curr_ixs:
                        # check available atoms
                        default_atoms = aa_default_atoms[prev_res]
                        intersection = list(set(default_atoms) - set(curr_atoms))
                        if intersection != [] or next_res_for_rm:
                            # missing atoms: means that there are atoms in the 'default' list
                            # that are not observed in the structure
                            drop_ixs += curr_ixs
                            next_res_for_rm = False
                        elif prev_seq == curr_seq:
                            # duplicated *_seq_id: means that the next residue has got the same
                            # *_seq_id as the previous res (could come from removing altlocs)
                            next_res_for_rm = True
                    # resetting variables
                    prev_res = curr_res
                    prev_seq = curr_seq
                    curr_ixs = [ix]
                    curr_atoms = [curr_atom]

    return table.drop(table.index[drop_ixs])


# TODO add backbone centroid
#     elif atoms == 'backbone_centroid':
#         cif_table = row_selector(
#             cif_table, 'label_atom_id', ('CA', 'N', 'C', 'O'))
#         cif_table = _residues_as_centroid(cif_table)
def residues_aggregation(table, agg_method='centroid', category='label'):
    """
    Gets the residues' atoms and their centroids (mean).

    :param table: pandas DataFrame object
    :param agg_method: current values: 'centroid', 'first', 'mean' and 'unique'
    :param category: data category to be used as precedence in _atom_site.*_*
        asym_id, seq_id and atom_id
    :return: returns a modified pandas DataFrame
    """

    agg_generic = agg_method
    agg_cols = ['pdbx_PDB_model_num', '{}_asym_id'.format(category),
                '{}_seq_id'.format(category)]
    if agg_method not in ['centroid', 'first', 'unique', 'mean']:
        raise ValueError('Method {} is not currently implemented...'
                         ''.format(agg_method))
    if agg_method == 'centroid' or agg_method == 'mean':
        agg_generic = 'first'
        agg_method = 'mean'
    columns_to_agg = {col: agg_generic if table[col].dtype == 'object' else agg_method
                      for col in table.columns if col not in agg_cols}
    columns_to_agg['id'] = 'first'
    table = table.groupby(by=agg_cols, as_index=False).agg(columns_to_agg)
    table = table.sort_values(by='id').reset_index()
    return table


def write_mmcif_from_table(table, filename, override=False):
    """
    Generic method that writes 'atom' lines in mmCIF format.

    :param table: pandas DataFrame object
    :param filename: path to the mmCIF file
    :param override: boolean
    :return: (side effects) writes to file
    """

    atom_lines = ['data_mmCIF_generated_by_ProteoFAV', 'loop_']
    atom_lines += ["_atom_site.{}".format(v) for v in list(table)]
    for i in table.index:
        line = ' '.join([str(v) for v in list(table.loc[i, :])])
        atom_lines.append(line)

    # write the final output
    if not os.path.exists(filename) or override:
        with open(filename, 'w') as outlines:
            outlines.write("\n".join(atom_lines))
    else:
        log.info("mmCIF for %s already available...", filename)
    return


def write_pdb_from_table(table, filename, override=False, category='auth'):
    """
    Generic method that writes 'atom' lines in PDB format.

    :param table: pandas DataFrame object
    :param filename: path to the PDB file
    :param override: boolean
    :param category: data category to be used as precedence in _atom_site.*_*
        asym_id, seq_id and atom_id
    :return: (side effects) writes to file
    """

    atom_lines = ['REMARK 100 PDB generated by ProteoFAV\n']
    atom_number = 0
    for i in table.index:
        atom_number += 1
        atom_lines.append(_get_atom_line(table=table, index=i,
                                         atom_number=atom_number,
                                         category=category))

    # write the final output
    if not os.path.exists(filename) or override:
        with open(filename, 'w') as outlines:
            outlines.write("".join(atom_lines))
    else:
        log.info("PDB for %s already available...", filename)
    return


def _get_atom_line(table, index, atom_number, category='auth'):
    """
    Returns an ATOM PDB-formatted string.
    (Based on code from the PDB module in Biopython.)

    :param table: pandas DataFrame object
    :param index: atom index
    :param atom_number: incremental number
    :param category: data category to be used as precedence in _atom_site.*_*
        asym_id, seq_id and atom_id
    :return: returns a PDB-formatted ATOM/HETATM line
    """

    ix = index
    record_type = table.loc[ix, 'group_PDB']
    if record_type == "ATOM":
        record_type = "ATOM  "

    """
    ATOM     16  CB  ASN A   2      22.780  31.612   8.039  1.00 97.98           C
    ATOM     17  CG  ASN A   2      23.735  31.870   9.167  1.00100.56           C
    ATOM     18  OD1 ASN A   2      23.345  32.366  10.218  1.00 97.84           O
    ATOM     19  ND2 ASN A   2      25.014  31.606   8.922  1.00106.62           N
    ATOM     20  H   ASN A   2      24.256  34.106   6.858  1.00  0.00           H
    ATOM     21 HD21 ASN A   2      25.654  31.751   9.644  1.00  0.00           H
    ATOM     22 HD22 ASN A   2      25.276  31.270   8.035  1.00  0.00           H
    """

    name = table.loc[ix, '{}_atom_id'.format(category)]
    if len(name) == 1:
        name = " {}  ".format(name.strip())
    elif len(name) == 2:
        name = " {} ".format(name.strip())
    elif len(name) == 3:
        name = " {}".format(name.strip())
    elif len(name) == 4:
        name = name.strip()

    altloc = table.loc[ix, 'label_alt_id']
    if altloc == ".":
        altloc = " "

    resname = table.loc[ix, '{}_comp_id'.format(category)]

    chain_id = table.loc[ix, '{}_asym_id'.format(category)]
    chain_id = chain_id[0]

    resseq = str(table.loc[ix, '{}_seq_id'.format(category)])

    icode = table.loc[ix, 'pdbx_PDB_ins_code']
    if icode == "?":
        icode = " "

    x = float(table.loc[ix, 'Cartn_x'])
    y = float(table.loc[ix, 'Cartn_y'])
    z = float(table.loc[ix, 'Cartn_z'])

    occupancy_str = "%6.2f" % float(table.loc[ix, 'occupancy'])

    bfactor = float(table.loc[ix, 'B_iso_or_equiv'])

    segid = ""

    element = table.loc[ix, 'type_symbol']
    element = element.strip().upper()

    charge = "  "

    values = (record_type, atom_number, name, altloc, resname, chain_id,
              resseq, icode, x, y, z, occupancy_str, bfactor, segid,
              element, charge)

    return _PDB_FORMAT % values


if __name__ == '__main__':
    pass
