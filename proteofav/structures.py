# -*- coding: utf-8 -*-

"""
Functions that handle the reading data files and extracting their information
as a pandas.DataFrame. Also include wrapper functions that select and index
the information. Prefers the use o the wrapper instead the private functions
for better error handling. Both levels are covered by test cases.
"""

import os
import logging
import pandas as pd
from string import ascii_uppercase

try:
    # python 2.7
    from StringIO import StringIO
except ImportError:
    from io import StringIO

from proteofav.config import defaults
from proteofav.utils import (fetch_from_url_or_retry, row_selector,
                             InputFileHandler, OutputFileHandler, GenericInputs,
                             constrain_column_types, exclude_columns, Downloader,
                             check_sequence)
from proteofav.library import pdbx_types, aa_default_atoms, scop_3to1

log = logging.getLogger('proteofav.config')

__all__ = ['parse_mmcif_atoms', 'residues_aggregation',
           'fetch_summary_properties_pdbe', 'get_preferred_assembly_id',
           'filter_structures', 'select_structures', 'write_mmcif_from_table', 'write_pdb_from_table',
           'read_structures', 'download_structures', 'write_structures', 'PDB', 'mmCIF']

UNIFIED_COL = ['pdbx_PDB_model_num', 'auth_asym_id', 'auth_seq_id']

PDB_FORMAT = "%s%5i %-4s%c%3s %c%4s%c   %8.3f%8.3f%8.3f%s%6.2f      %4s%2s%2s\n"


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

    log.debug("Parsing mmCIF atoms from lines...")

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
    log.debug("Removed columns from mmCIF: {}..."
              "".format(', '.join(list(excluded_cols))))

    # enforce some specific column types
    table = constrain_column_types(table, col_type_dict=pdbx_types)

    if table.empty:
        raise ValueError('mmCIF file {} resulted in a empty Dataframe'
                         ''.format(filename))
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

    log.debug("Parsing PDB atoms from lines...")

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
    log.debug("Removed columns from PDB: {}..."
              "".format(', '.join(list(excluded_cols))))

    # enforce some specific column types
    table = constrain_column_types(table, col_type_dict=pdbx_types)

    if table.empty:
        raise ValueError('PDB file {} resulted in a empty Dataframe'
                         ''.format(filename))
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


def _add_mmcif_res_full(table):
    """
    Utility that adds a new column to the table.
    Adds a new column with the 'full res' (i.e. seq_id + ins_code).

    :param table: pandas DataFrame object
    :return: returns a modified pandas DataFrame
    """

    # adds both 'label' and 'auth' entries
    if 'label_seq_id' in table:
        table['label_seq_id_full'] = (table['label_seq_id'] +
                                      table['pdbx_PDB_ins_code'].str.replace('?', ''))
    if 'auth_seq_id' in table:
        table['auth_seq_id_full'] = (table['auth_seq_id'] +
                                     table['pdbx_PDB_ins_code'].str.replace('?', ''))
    return table


def _add_mmcif_atom_altloc(table):
    """
    Utility that adds new columns to the table.
    adds: 'label_atom_altloc_id', and 'auth_atom_altloc_id' which is a
        string join between '*_atom_id' + 'label_alt_id'

    :param table: pandas DataFrame object
    :return: returns a modified pandas DataFrame
    """

    def join_atom_altloc(table, category='auth'):
        atom = table['{}_atom_id'.format(category)]
        altloc = table['label_alt_id']

        new_column = (atom + '.' + altloc)
        has_no_alt = altloc.isin(['.', '', ' '])
        new_column[has_no_alt] = atom.loc[has_no_alt]

        return new_column

    # NB. Use of assign prevents inplace modification
    table = table.assign(label_atom_altloc_id=join_atom_altloc(table, 'label'))
    table = table.assign(auth_atom_altloc_id=join_atom_altloc(table, 'auth'))
    return table


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


def _remove_partial_residues(table, category='auth'):
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


def residues_aggregation(table, agg_method='centroid', category='auth'):
    """
    Gets the residues' atoms and their centroids (mean).

    :param table: pandas DataFrame object
    :param agg_method: current values: 'centroid', 'backbone_centroid'',
        first', 'mean' and 'unique'
    :param category: data category to be used as precedence in _atom_site.*_*
        asym_id, seq_id and atom_id
    :return: returns a modified pandas DataFrame
    """

    agg_generic = agg_method
    agg_cols = ['pdbx_PDB_model_num', '{}_asym_id'.format(category),
                '{}_seq_id'.format(category)]
    if agg_method not in ['centroid', 'first', 'unique', 'mean', 'backbone_centroid']:
        raise ValueError('Method {} is not currently implemented...'.format(agg_method))
    if agg_method == 'backbone_centroid':
        table = row_selector(table, '{}_atom_id'.format(category), ('CA', 'N', 'C', 'O'))
        agg_method = 'centroid'
    if agg_method == 'centroid' or agg_method == 'mean':
        agg_generic = 'first'
        agg_method = 'mean'
    columns_to_agg = {col: agg_generic if table[col].dtype == 'object' else agg_method
                      for col in table.columns if col not in agg_cols}
    columns_to_agg['id'] = 'first'
    table = table.groupby(by=agg_cols, as_index=False).agg(columns_to_agg)
    table = table.sort_values(by='id').reset_index()
    log.debug("Aggregated residues based on agg_method '{}'".format(agg_method))
    return table


def write_mmcif_from_table(table, filename, overwrite=False):
    """
    Generic method that writes 'atom' lines in mmCIF format.

    :param table: pandas DataFrame object
    :param filename: path to the mmCIF file
    :param overwrite: boolean
    :return: (side effects) writes to file
    """

    atom_lines = ['data_mmCIF_generated_by_ProteoFAV', 'loop_']
    atom_lines += ["_atom_site.{}".format(v) for v in list(table)]
    for i in table.index:
        line = ' '.join([str(v) for v in list(table.loc[i, :])])
        atom_lines.append(line)

    # write the final output
    if not os.path.exists(filename) or overwrite:
        with open(filename, 'w') as outlines:
            outlines.write("\n".join(atom_lines))
        log.debug("Wrote mmCIF-formatted file {}...".format(filename))
    else:
        log.info("mmCIF for %s already available...", filename)
    return


def write_pdb_from_table(table, filename, overwrite=False, category='auth'):
    """
    Generic method that writes 'atom' lines in PDB format.

    :param table: pandas DataFrame object
    :param filename: path to the PDB file
    :param overwrite: boolean
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
    if not os.path.exists(filename) or overwrite:
        with open(filename, 'w') as outlines:
            outlines.write("".join(atom_lines))
        log.debug("Wrote PDB-formatted file {}...".format(filename))
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

    return PDB_FORMAT % values


def fetch_summary_properties_pdbe(identifier, retry_in=(429,)):
    """
    Queries the PDBe api to get summary validation report.

    :param identifier: PDB ID
    :param retry_in: http code for retrying connections
    :return: Requests object
    """
    pdbe_endpoint = "pdb/entry/summary/{}".format(identifier)
    url = defaults.api_pdbe + pdbe_endpoint
    response = fetch_from_url_or_retry(url, json=True, retry_in=retry_in)
    return response


def get_preferred_assembly_id(identifier):
    """
    Gets the preferred assembly id for the given PDB ID, from the PDBe API.

    :param identifier: PDB ID
    :return: str
    """

    # getting the preferred biological assembly from the PDBe API
    try:
        data = fetch_summary_properties_pdbe(identifier).json()
    except Exception as e:
        log.error("Something went wrong for {}... {}".format(identifier, e))
    try:
        nassemblies = data[identifier][0]["assemblies"]
        if len(nassemblies) > 1:
            for entry in nassemblies:
                if entry["preferred"]:
                    pref_assembly = entry["assembly_id"]
                    break
        else:
            pref_assembly = data[identifier][0]["assemblies"][0]["assembly_id"]
    except Exception:
        pref_assembly = "1"

    bio_best = str(pref_assembly)
    return bio_best


def get_sequence(table, category='auth', ambiguous='X'):
    """
    Get the sequence for the PDBx table.

    :param table: pandas DataFrame from PDB/mmCIF.read()
    :param category: data category to be used as precedence in _atom_site.*_*
        asym_id, seq_id and atom_id
    :param ambiguous: (str) 1-letter symbol for ambiguous residues
    :returns: (str) 1-letter sequence
    """

    # assumes it's a valid PDBx table
    if isinstance(table, pd.DataFrame):
        sequence = ""
        for ix in table.index:
            aa3 = table.loc[ix, "{}_comp_id".format(category)]
            aa1 = scop_3to1[aa3]
            if len(aa1) == 1:
                sequence += aa1
            else:
                sequence += 'X'
        sequence = check_sequence(sequence=sequence, ambiguous=ambiguous)
    else:
        return ValueError("Pandas DataFrame is not valid...")
    return sequence


##############################################################################
# Public methods
##############################################################################
def select_structures(identifier=None, excluded_cols=None,
                      bio_unit=False, bio_unit_preferred=False, bio_unit_id="1",
                      overwrite=False, **kwargs):
    """
    Produce table read from mmCIF file.

    :param identifier: PDB/mmCIF accession ID
    :param excluded_cols: option to exclude mmCIF columns
    :param bio_unit: (boolean) whether to get BioUnit instead of AsymUnit
    :param bio_unit_preferred: (boolean) if True fetches the 'preferred'
        biological assembly instead
    :param bio_unit_id: biological unit identifier
    :param overwrite: boolean
    :return: returns a pandas DataFrame
    """

    if bio_unit:
        filename = os.path.join(defaults.db_mmcif, "{}_bio.cif".format(identifier))
    else:
        filename = os.path.join(defaults.db_mmcif, "{}.cif".format(identifier))

    download_structures(identifier=identifier, filename=filename,
                        bio_unit=bio_unit, bio_unit_preferred=bio_unit_preferred,
                        bio_unit_id=bio_unit_id, overwrite=overwrite)

    table = read_structures(filename=filename, excluded_cols=excluded_cols)

    table = filter_structures(table=table, excluded_cols=excluded_cols, **kwargs)

    table = constrain_column_types(table, col_type_dict=pdbx_types)

    # id is the atom identifier and it is need for all atoms tables.
    if table[UNIFIED_COL + ['id']].duplicated().any():
        log.error('Failed to find unique index for {}'.format(filename))

    return table


def filter_structures(table, excluded_cols=None,
                      models='first', chains=None, res=None, res_full=None,
                      comps=None, atoms=None, lines=None, category='auth',
                      residue_agg=False, agg_method='centroid',
                      add_res_full=True, add_atom_altloc=False, reset_atom_id=True,
                      remove_altloc=False, remove_hydrogens=True, remove_partial_res=False):
    """
    Filter and residue aggregation for mmCIF and PDB Pandas Dataframes.
    Improves consistency  ATOM and HETATM lines from mmCIF or PDB files.

    :param table: pandas DataFrame object
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
    :return: returns a pandas DataFrame
    """

    # selections / filtering
    # excluding columns
    table = exclude_columns(table, excluded=excluded_cols)

    # excluding rows
    # select the models
    # if only first model (>1 in NMR structures)
    if models:
        table = row_selector(table, 'pdbx_PDB_model_num', models)
        log.debug("mmCIF/PDB table filtered by pdbx_PDB_model_num...")

    # select chains
    if chains:
        table = row_selector(table, '{}_asym_id'.format(category), chains)
        log.debug("mmCIF/PDB table filtered by %s_asym_id...", category)

    # select lines
    if lines:
        table = row_selector(table, 'group_PDB', lines)
        log.debug("mmCIF/PDB table filtered by group_PDB...")

    # table modular extensions or selections
    if add_res_full:
        table = _add_mmcif_res_full(table)
        log.debug("mmCIF/PDB added full res (res + ins_code)...")

    if add_atom_altloc:
        table = _add_mmcif_atom_altloc(table)
        log.debug("mmCIF/PDB added full atom (atom + altloc)...")

    if remove_hydrogens:
        table = row_selector(table, key='type_symbol', value='H', reverse=True)
        log.debug("mmCIF/PDB removed existing hydrogens...")

    if remove_altloc:
        table = _remove_multiple_altlocs(table)
        reset_atom_id = True
        log.debug("mmCIF/PDB removed altlocs...")

    if remove_partial_res:
        table = _remove_partial_residues(table)
        log.debug("mmCIF/PDB removed incomplete residues...")

    if reset_atom_id:
        table.reset_index(inplace=True)
        table = table.drop(['index'], axis=1)
        table['id'] = table.index + 1
        log.debug("mmCIF/PDB reset atom numbers...")

    # excluding rows
    # select residues
    if res:
        table = row_selector(table, '{}_seq_id'.format(category), res)
        log.debug("mmCIF/PDB table filtered by %s_seq_id...", category)

    if res_full:
        table = row_selector(table, '{}_seq_id_full'.format(category), res_full)
        log.debug("mmCIF/PDB table filtered by %s_seq_id_full...", category)

    # select amino acids/molecules
    if comps:
        table = row_selector(table, '{}_comp_id'.format(category), comps)
        log.debug("mmCIF/PDB table filtered by %s_comp_id...", category)

    # select atoms
    if atoms == 'centroid' or atoms == 'backbone_centroid':
        table = residues_aggregation(table, agg_method=atoms)
    elif atoms:
        table = row_selector(table, '{}_atom_id'.format(category), atoms)
        log.debug("mmCIF/PDB table filtered by %s_atom_id...", category)

    if residue_agg:
        table = residues_aggregation(table, agg_method=agg_method,
                                     category=category)

    if table.empty:
        raise ValueError("The filters resulted in an empty DataFrame...")
    return table


def read_structures(filename=None, input_format=None, excluded_cols=None,
                    pdb_fix_ins_code=True, pdb_fix_label_alt_id=True, pdb_fix_type_symbol=True):
    """
    Parse ATOM and HETATM lines from mmCIF or PDB files.

    :param filename: path to the mmCIF/PDB file
    :param input_format: either 'mmcif', 'cif', 'pdb' or 'ent'
    :param excluded_cols: option to exclude mmCIF columns
    :param pdb_fix_ins_code: boolean
    :param pdb_fix_label_alt_id: boolean
    :param pdb_fix_type_symbol: boolean
    :return: returns a pandas DataFrame
    """

    InputFileHandler(filename)

    # try to guess the correct format
    if input_format is None:
        if filename.endswith('.pdb') or filename.endswith('.ent'):
            input_format = "pdb"
        elif filename.endswith('.cif') or filename.endswith('.mmcif'):
            input_format = "mmcif"
        else:
            raise ValueError("Could not guess the format of the input file... "
                             "Please define it by passing 'input_format'='<name>'")
        log.debug("Input format seems to be {}...".format(input_format))

    if input_format == 'mmcif':
        table = parse_mmcif_atoms(filename, excluded_cols=excluded_cols)
    elif input_format == 'pdb':
        table = parse_pdb_atoms(filename, excluded_cols=excluded_cols,
                                fix_ins_code=pdb_fix_ins_code,
                                fix_label_alt_id=pdb_fix_label_alt_id,
                                fix_type_symbol=pdb_fix_type_symbol)
    return table


def write_structures(table=None, filename=None, overwrite=False,
                     output_format=None, category='auth'):
    """
    Writes 'ATOM' lines in PDB or mmCIF format.

    :param table: pandas DataFrame object
    :param filename: path to the mmCIF/PDB file
    :param overwrite: boolean
    :param category: data category to be used as precedence in _atom_site.*_*
        asym_id, seq_id and atom_id (generally either 'label' or 'auth')
    :param output_format: either 'mmcif', 'cif', 'pdb' or 'ent'
    :return: (side effects) writes to a file
    """

    if not os.path.exists(filename) or overwrite:
        OutputFileHandler(filename, overwrite=overwrite)

        if output_format is None:
            if filename.endswith('.pdb') or filename.endswith('.ent'):
                output_format = "pdb"
            elif filename.endswith('.cif') or filename.endswith('.mmcif'):
                output_format = "mmcif"
            else:
                raise ValueError("Could not guess the format of the input file... "
                                 "Please define it by passing 'input_format'='<name>'")

        if output_format == 'mmcif' or output_format == 'cif':
            write_mmcif_from_table(filename=filename, table=table,
                                   overwrite=overwrite)
        elif output_format == 'pdb' or output_format == 'ent':
            write_pdb_from_table(filename=filename, table=table,
                                 overwrite=overwrite, category=category)


def download_structures(identifier, filename, output_format='mmcif',
                        bio_unit=False, bio_unit_preferred=False, bio_unit_id="1",
                        mmcif_version="_updated", overwrite=False):
    """
    Downloads a structure from the PDBe to the filesystem.
    BioUnits only work for mmCIF. If `pdb=True` the bio_unit arguments
    are ignored.

    :param identifier: (str) PDB accession ID
    :param filename: path to the mmCIF/PDB file
    :param output_format: either 'mmcif', 'cif', 'pdb' or 'ent'
    :param bio_unit: (boolean) whether to get BioUnit instead of AsymUnit
    :param bio_unit_preferred: (boolean) if True fetches the 'preferred'
        biological assembly instead
    :param bio_unit_id: biological unit identifier
    :param overwrite: boolean
    :param mmcif_version: "_updated" or simply ""
    :return: (side effects) writes to a file
    """

    if output_format is None:
        if filename.endswith('.pdb') or filename.endswith('.ent'):
            output_format = "pdb"
        elif filename.endswith('.cif') or filename.endswith('.mmcif'):
            output_format = "mmcif"
        else:
            raise ValueError("Could not guess the format of the input file... "
                             "Please define it by passing 'input_format'='<name>'")

    decompress = False
    if output_format == 'mmcif' or output_format == 'cif':
        if bio_unit:
            # atom lines only?
            # url_endpoint = ("{}-assembly-{}_atom_site.cif.gz".format(identifier, pref))
            if bio_unit_preferred:
                pref = get_preferred_assembly_id(identifier=identifier)
            else:
                pref = bio_unit_id
            url_endpoint = ("{}-assembly-{}.cif.gz".format(identifier, pref))
            decompress = True
            url_root = defaults.bio_fetch
            url = url_root + url_endpoint
        else:
            # original mmCIF?
            # url_endpoint = "download/{}.cif".format(pdbid)
            url_endpoint = "download/{}{}.cif".format(identifier, mmcif_version)
            url_root = defaults.pdbe_fetch
            url = url_root + url_endpoint

    elif output_format == 'pdb' or output_format == 'ent':
        url_endpoint = "download/pdb{}.ent".format(identifier)
        url_root = defaults.pdbe_fetch
        url = url_root + url_endpoint

    else:
        raise ValueError("Could not guess the format of the output file... "
                         "Please define it by passing 'output_format'='<name>'")

    Downloader(url=url, filename=filename,
               decompress=decompress, overwrite=overwrite)


class Structure(GenericInputs):
    def read(self, filename=None, **kwargs):
        filename = self._get_filename(filename)
        self.table = read_structures(filename=filename, **kwargs)
        return self.table

    def write(self, table=None, filename=None, **kwargs):
        table = self._get_table(table)
        filename = self._get_filename(filename)
        return write_structures(table=table, filename=filename, **kwargs)

    def download(self, identifier=None, filename=None, **kwargs):
        identifier = self._get_identifier(identifier)
        filename = self._get_filename(filename)
        return download_structures(identifier=identifier, filename=filename, **kwargs)

    def select(self, identifier=None, **kwargs):
        identifier = self._get_identifier(identifier)
        self.table = select_structures(identifier=identifier, **kwargs)
        return self.table


PDB = Structure()
mmCIF = Structure()
