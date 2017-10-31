# -*- coding: utf-8 -*-

import os
import logging
import pandas as pd
from string import digits
from string import ascii_uppercase

try:
    # python 2.7
    from StringIO import StringIO
except ImportError:
    from io import StringIO

from proteofav.structures import select_structures
from proteofav.utils import (row_selector, InputFileHandler,
                             constrain_column_types, exclude_columns,
                             GenericInputs, Downloader)
from proteofav.library import (scop_3to1, dssp_types, aa_codes_1to3_extended)
from proteofav.library import (ASA_Miller, ASA_Wilke, ASA_Sander)

from proteofav.config import defaults

log = logging.getLogger('proteofav.config')

__all__ = ['parse_dssp_residues', '_import_dssp_chains_ids', 'select_dssp',
           'filter_dssp', 'get_rsa', 'get_rsa_class', 'download_dssp', 'DSSP']


def parse_dssp_residues(filename, excluded_cols=None):
    """
    Parse lines of the DSSP file to get entries for every Residue
    in each CHAIN. The hierarchy is maintained. CHAIN->RESIDUE->[...].

    :param filename: path to the DSSP file
    :param excluded_cols: list of columns to be excluded
    :return: returns a pandas DataFrame
    """

    log.debug("Parsing DSSP from lines...")

    # example lines with some problems
    """
      #  RESIDUE AA STRUCTURE BP1 BP2  ACC     N-H-->O    O-->H-N    N-H-->O    O-->H-N    TCO  KAPPA ALPHA  PHI   PSI    X-CA   Y-CA   Z-CA
        1    1 A M              0   0  127      0, 0.0   345,-0.1     0, 0.0     3,-0.1   0.000 360.0 360.0 360.0 162.0  -18.7   21.6  -55.4
        2    2 A R        +     0   0  117      1,-0.1    28,-0.4   343,-0.1     2,-0.3   0.455 360.0  81.5-136.8 -28.7  -17.0   22.3  -52.1

      381  394 A K              0   0  125     -2,-0.4   -21,-0.1   -21,-0.2    -2,-0.0  -0.421 360.0 360.0  64.1 360.0  -22.5   44.2  -25.4
      382        !*             0   0    0      0, 0.0     0, 0.0     0, 0.0     0, 0.0   0.000 360.0 360.0 360.0 360.0    0.0    0.0    0.0
      383    1 A M              0   0  127      0, 0.0   345,-0.1     0, 0.0     3,-0.1   0.000 360.0 360.0 360.0 162.0  -10.0   71.4  -55.4

    10278  103 H H  E     -XZ1023010269W  69     -9,-2.3    -9,-2.2    -2,-0.3     2,-1.0  -0.884  22.6-128.4-108.1 141.6  -97.0   28.7  112.2
    10279  104 H I  E     +XZ1022910268W   0    -50,-2.2   -50,-0.6    -2,-0.4   -11,-0.3  -0.801  30.6 175.4 -90.4  95.6  -98.5   32.0  111.3
    10280  105 H L  E     +     0   0   21    -13,-1.7   -55,-2.5    -2,-1.0     2,-0.3   0.812  62.6   4.9 -70.5 -35.5  -96.3   34.5  113.1

    # missing segment break
      145        !              0   0    0      0, 0.0     0, 0.0     0, 0.0     0, 0.0   0.000 360.0 360.0 360.0 360.0    0.0    0.0    0.0

    # chain break
      382        !*             0   0    0      0, 0.0     0, 0.0     0, 0.0     0, 0.0   0.000 360.0 360.0 360.0 360.0    0.0    0.0    0.0
    """

    InputFileHandler(filename)

    lines = []
    parse = False
    with open(filename) as inlines:
        for line in inlines:
            line = line.rstrip()
            if parse:
                lines.append(line)
            if line.startswith("  #"):
                parse = True
    lines = "\n".join(lines)

    # column width descriptors
    header = ("LINE", "RES", "RES_FULL", "INSCODE", "CHAIN", "AA", "SS", "STRUCTURE",
              "BP1", "BP2", "BP2_CHAIN", "ACC",
              "NH_O_1", "NH_O_1_nrg", "O_HN_1", "O_HN_1_nrg",
              "NH_O_2", "NH_O_2_nrg", "O_HN_2", "O_HN_2_nrg",
              "TCO", "KAPPA", "ALPHA", "PHI", "PSI",
              "X-CA", "Y-CA", "Z-CA")

    widths = ((0, 5), (5, 10), (5, 11), (10, 11), (11, 12), (12, 15), (16, 17), (17, 25),
              (25, 29), (29, 33), (33, 34), (34, 38),
              (38, 45), (46, 50), (50, 56), (57, 61),
              (61, 67), (68, 72), (72, 78), (79, 84),
              (85, 91), (91, 97), (97, 103), (103, 109), (109, 115),
              (115, 123), (123, 130), (130, 137))

    all_str = {key: str for key in header}
    table = pd.read_fwf(StringIO(lines), names=header, colspecs=widths,
                        compression=None, converters=all_str, keep_default_na=False)

    # excluding columns
    if excluded_cols is None:
        excluded_cols = ("LINE", "STRUCTURE", "BP1", "BP2", "BP2_CHAIN",
                         "NH_O_1", "NH_O_1_nrg", "O_HN_1", "O_HN_1_nrg",
                         "NH_O_2", "NH_O_2_nrg", "O_HN_2", "O_HN_2_nrg",
                         "X-CA", "Y-CA", "Z-CA")

    table = exclude_columns(table, excluded=excluded_cols)

    # enforce some specific column types
    table = constrain_column_types(table, col_type_dict=dssp_types)

    if table.empty:
        raise ValueError('DSSP file {} resulted in a empty Dataframe'
                         ''.format(filename))
    return table


def _import_dssp_chains_ids(pdb_id):
    """Imports mmCIF chain identifier to DSSP.

    :param pdb_id:
    :return: DSSP table with corrected chain ids.
    """
    dssp_table = select_dssp(pdb_id)
    cif_table = select_structures(pdb_id)
    cif_seq = cif_table.auth_comp_id.apply(scop_3to1.get)
    dssp_has_seq = dssp_table.aa.isin(scop_3to1.values())
    dssp_seq = dssp_table.aa[dssp_has_seq]
    # Import only if the sequences are identical
    if not (cif_seq == dssp_seq).all():
        err = ('Inconsitent DSSP / mmCIF sequence for {} protein structure cannot be fixed'
               'by import_dssp_chains_ids')
        raise ValueError(err.format(pdb_id))
    dssp_table.loc[dssp_has_seq, 'CHAIN'] = cif_table.auth_asym_id
    return dssp_table


def _add_dssp_rsa(data, method="Sander"):
    """
    Utility that adds a new column to the table.
    Adds a new column with Relative Solvent Accessibility (RSA).

    :param data: pandas DataFrame object
    :param method: name of the method
    :return: returns a modified pandas DataFrame
    """

    table = data
    rsas = []
    for i in table.index:
        rsas.append(get_rsa(table.loc[i, "ACC"], table.loc[i, "AA"],
                            method=method))
    table["RSA"] = rsas
    return table


def _add_dssp_full_chain(data):
    """
    Utility that adds a new column to the table.
    Specific to DSSP outputs that are generated from mmCIF files containing
    multiple char chain IDs (e.g. 'AA' and 'BA'). These are found in the
    Biological Unit mmCIF structures.

    :param data: pandas DataFrame object
    :return: returns a modified pandas DataFrame
    """

    # BioUnits chain naming seems to follow the pattern:
    # chain A becomes AA then
    # AA->AZ then A0->A9 [A-Z then 0-9] and then AAA->AAZ and AA0->AA9
    # then ABA->ABZ and AB0->AB9
    alpha1 = [k for k in ascii_uppercase + digits]
    alpha2 = ['A' + k for k in alpha1]
    alpha3 = ['B' + k for k in alpha1]
    new_alphabet = alpha1 + alpha2 + alpha3

    table = data
    chains_full = []
    c = -1
    for ix in table.index:
        chain_id = table.loc[ix, "CHAIN"]
        aa_id = table.loc[ix, "AA"]
        if aa_id == "!*":
            if table.loc[ix - 1, "CHAIN"] == table.loc[ix + 1, "CHAIN"]:
                c += 1
            else:
                c = -1
        if c != -1 and aa_id != "!*" and aa_id != "!":
            if c >= len(new_alphabet):
                raise IndexError('Alphabet needs update to accommodate '
                                 'such high number of chains...')
            chain_id += new_alphabet[c]
        chains_full.append(chain_id)
    if not chains_full:
        table["CHAIN_FULL"] = table["CHAIN"]
    else:
        table["CHAIN_FULL"] = chains_full
    return table


def _add_dssp_rsa_class(data, rsa_col='RSA'):
    """
    Utility that adds a new column to the table.
    Adds a new column with Relative Solvent Accessibility (RSA) classes.

    :param data: pandas DataFrame object
    :param rsa_col: column name
    :return: returns a modified pandas DataFrame
    """

    table = data
    rsas_class = []
    for i in table.index:
        rsas_class.append(get_rsa_class(table.loc[i, "{}".format(rsa_col)]))
    table["{}_CLASS".format(rsa_col)] = rsas_class
    return table


def _add_dssp_ss_reduced(data):
    """
    Utility that adds a new column to the table.
    Adds a reduced-stated Secondary Structure (SS).

    :param data: pandas DataFrame object
    :return: returns a modified pandas DataFrame
    """

    table = data
    alphas = ['H']
    betas = ['E']
    coils = ['G', 'I', 'B', 'C', 'T', 'S', '', ' ']

    # replace some NaN with custom strings
    # table['SS'] = table.SS.fillnan('-')
    sss = []
    for ix in table.index:
        ss = table.loc[ix, "SS"]
        if ss in alphas:
            ss = 'H'
        elif ss in betas:
            ss = 'E'
        elif ss in coils:
            ss = 'C'
        else:
            ss = '-'
        sss.append(ss)

    table["SS_CLASS"] = sss

    return table


def get_rsa(acc, resname, method="Sander"):
    """
    Computes Relative Solvent Accessibility (RSA) from an input
    DSSP ACC value, and according to ASA standard values.

    :param acc: (int) DSSP ACC
    :param resname: single letter residue name
    :param method: name of the method
    :return: (float) RSA value
    """

    rsa = ""
    try:
        acc = float(acc)
    except ValueError:
        return rsa

    try:
        assert len(resname) == 1
    except AssertionError:
        return rsa

    if method == "Miller":
        sasa = ASA_Miller
    elif method == "Wilke":
        sasa = ASA_Wilke
    elif method == "Sander":
        sasa = ASA_Sander
    else:
        raise ValueError("Method {} is not implemented...".format(method))

    try:
        rsa = round((acc / sasa[aa_codes_1to3_extended[resname]] * 100), 3)
    except KeyError:
        return rsa

    log.debug("Computed RSA from DSSP ACC with method '{}'".format(method))
    return rsa


def get_rsa_class(rsa, lower_threshold=5.0, upper_threshold=25.0):
    """
    Gets a class based on the RSA value

    :param rsa: (float) RSA score or string
    :param lower_threshold: (float) threshold
    :param upper_threshold: (float) threshold
    :return: RSA class
    """
    rsa_class = '-'
    try:
        rsa = float(rsa)
        # surface is rsa >= 25% (value = 2)
        # exposed is rsa >= 5% and rsa < 25% (value = 1)
        # core is rsa < 5% (value = 0)
        if rsa >= upper_threshold:
            rsa_class = 'Surface'
        elif lower_threshold <= rsa < upper_threshold:
            rsa_class = 'Part. Exposed'
        else:
            rsa_class = 'Core'
    except ValueError:
        # returns a string
        pass
    log.debug("RSA Classes are based on RSA cut-offs: {} and {}"
              "".format(lower_threshold, upper_threshold))
    return rsa_class


def select_dssp(identifier, excluded_cols=None, overwrite=False, **kwargs):
    """
    Produce table from DSSP file output.

    :param identifier: PDB/mmCIF accession ID
    :param excluded_cols: option to exclude DSSP columns
    :param overwrite: boolean
    :return: returns a pandas DataFrame
    """

    filename = os.path.join(defaults.db_dssp, "{}.dssp".format(identifier))

    download_dssp(identifier=identifier, filename=filename, overwrite=overwrite)

    table = parse_dssp_residues(filename=filename, excluded_cols=excluded_cols)

    table = filter_dssp(table=table, excluded_cols=excluded_cols, **kwargs)
    table = constrain_column_types(table, col_type_dict=dssp_types)

    if table.duplicated(['RES_FULL', 'CHAIN']).any():
        log.warning('DSSP file for {} has not unique index'.format(identifier))
    return table


def filter_dssp(table, excluded_cols=None,
                chains=None, chains_full=None, res=None,
                add_full_chain=True, add_ss_reduced=False,
                add_rsa=True, rsa_method="Sander", add_rsa_class=False,
                reset_res_id=False):
    """
    Filter for DSSP Pandas Dataframes.

    :param table: pandas DataFrame object
    :param excluded_cols: option to exclude DSSP columns
    :param chains: (tuple) chain IDs or None
    :param chains_full: (tuple) alternative chain IDs or None
    :param res: (tuple) res IDs or None
    :param add_full_chain: boolean
    :param add_ss_reduced: boolean
    :param add_rsa: boolean
    :param rsa_method: "Sander",  "Miller" or "Wilke"
    :param add_rsa_class: boolean
    :param reset_res_id: boolean
    :return: returns a pandas DataFrame
    """

    # selections / filtering
    # excluding columns
    table = exclude_columns(table, excluded=excluded_cols)

    # table modular extensions
    if add_full_chain:
        table = _add_dssp_full_chain(table)
        log.debug("DSSP added full chain...")

    table['SS'] = table.SS.fillna('-')
    if add_ss_reduced:
        table = _add_dssp_ss_reduced(table)
        log.debug("DSSP added reduced SS...")

    if add_rsa:
        table = _add_dssp_rsa(table, method=rsa_method)
        log.debug("DSSP added RSA...")

    if add_rsa_class:
        table = _add_dssp_rsa_class(table)
        log.debug("DSSP added RSA class...")

    # drop missing residues ("!")  and chain breaks ("!*")
    table = table[table['AA'] != '!']
    table = table[table['AA'] != '!*']

    # excluding rows
    if chains is not None:
        table = row_selector(table, 'CHAIN', chains)
        log.debug("DSSP table filtered by CHAIN...")

    if chains_full is not None:
        table = row_selector(table, 'CHAIN_FULL', chains_full)
        log.debug("DSSP table filtered by CHAIN_FULL...")

    if res is not None:
        table = row_selector(table, 'RES', res)
        log.debug("DSSP table filtered by RES...")

    if reset_res_id:
        table.reset_index(inplace=True)
        table = table.drop(['index'], axis=1)
        table['LINE'] = table.index + 1
        log.debug("DSSP reset residue number...")

    if table.empty:
        raise ValueError("The filters resulted in an empty DataFrame...")
    return table


def download_dssp(identifier, filename, overwrite=False):
    """
    Downloads a pre-computed DSSP from the CMBI Netherlands FTP
    to the filesystem.

    :param identifier: (str) PDB accession ID
    :param filename: path to the DSSP file
    :param overwrite: (boolean)
    :return: (side effects) output file path
    """

    url_root = defaults.dssp_fetch
    url_endpoint = "{}.dssp".format(identifier)
    url = url_root + url_endpoint
    Downloader(url=url, filename=filename,
               decompress=False, overwrite=overwrite)


class DSSP(GenericInputs):
    def read(self, filename=None, **kwargs):
        filename = self._get_filename(filename)
        self.table = parse_dssp_residues(filename=filename, **kwargs)
        return self.table

    def download(self, identifier=None, filename=None, **kwargs):
        identifier = self._get_identifier(identifier)
        filename = self._get_filename(filename)
        return download_dssp(identifier=identifier, filename=filename, **kwargs)

    def select(self, identifier=None, **kwargs):
        identifier = self._get_identifier(identifier)
        self.table = select_dssp(identifier=identifier, **kwargs)
        return self.table


DSSP = DSSP()
