# -*- coding: utf-8

import logging
from string import digits
from string import ascii_uppercase

from proteofav.parsers import parse_dssp_from_file
from proteofav.downloaders import Downloader

from proteofav.utils import InputFileHandler
from proteofav.utils import OutputFileHandler
from proteofav.utils import GenericInput
from proteofav.utils import exclude_columns
from proteofav.utils import constrain_column_types
from proteofav.utils import row_selector
from proteofav.library import dssp_types
from proteofav.library import aa_codes_1to3_extended
from proteofav.library import (ASA_Miller, ASA_Wilke, ASA_Sander)

from proteofav.config import defaults as config

log = logging.getLogger('proteofav')


class DSSP(GenericInput):
    def read(self, filename=None, excluded_cols=None,
             chains=None, chains_full=None, res=None,
             add_full_chain=True, add_ss_reduced=False,
             add_rsa=True, rsa_method="Sander", add_rsa_class=False,
             reset_res_id=False):
        """
        Parse lines of the DSSP file to get entries for every Residue
        in each CHAIN. The hierarchy is maintained. CHAIN->RESIDUE->[...].

        :param filename: path to the DSSP file
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

        filename = self._get_filename(filename)
        InputFileHandler(filename)

        if excluded_cols is None:
            excluded_cols = ("LINE", "STRUCTURE", "BP1", "BP2", "BP2_CHAIN",
                             "NH_O_1", "NH_O_1_nrg", "O_HN_1", "O_HN_1_nrg",
                             "NH_O_2", "NH_O_2_nrg", "O_HN_2", "O_HN_2_nrg",
                             "X-CA", "Y-CA", "Z-CA")

        table = parse_dssp_from_file(filename)

        # excluding columns
        table = exclude_columns(table, excluded=excluded_cols)

        # table modular extensions
        if add_full_chain:
            table = _add_dssp_full_chain(table)
            log.info("DSSP added full chain...")

        table['SS'] = table.SS.fillna('-')
        if add_ss_reduced:
            table = _add_dssp_ss_reduced(table)
            log.info("DSSP added reduced SS...")

        if add_rsa:
            table = _add_dssp_rsa(table, method=rsa_method)
            log.info("DSSP added RSA...")

        if add_rsa_class:
            table = _add_dssp_rsa_class(table)
            log.info("DSSP added RSA class...")

        # drop missing residues ("!")  and chain breaks ("!*")
        table = table[table['AA'] != '!']
        table = table[table['AA'] != '!*']

        if reset_res_id:
            table.reset_index(inplace=True)
            table = table.drop(['index'], axis=1)
            table['LINE'] = table.index + 1
            log.info("DSSP reset residue number...")

        # excluding rows
        if chains is not None:
            table = row_selector(table, 'CHAIN', chains)
            log.info("DSSP table filtered by CHAIN...")

        if chains_full is not None:
            table = row_selector(table, 'CHAIN_FULL', chains_full)
            log.info("DSSP table filtered by CHAIN_FULL...")

        if res is not None:
            table = row_selector(table, 'RES', res)
            log.info("DSSP table filtered by RES...")

        # enforce some specific column types
        table = constrain_column_types(table, dssp_types)

        if table.empty:
            raise ValueError('{} resulted in an empty DataFrame...'.format(filename))

        self.table = table
        return self.table

    def download(self, identifier=None, filename=None,
                 overwrite=False, decompress=True):
        """
        Downloads a pre-computed DSSP from the CMBI Netherlands FTP
        to the filesystem.

        :param identifier: (str) PDB accession ID
        :param filename: path to the DSSP file
        :param overwrite: (boolean)
        :param decompress: (boolean) Decompresses the file
        :return: (side effects) output file path
        """

        identifier = self._get_identifier(identifier)
        filename = self._get_filename(filename)
        OutputFileHandler(filename, overwrite=overwrite)

        url_root = config.ftp_dssp
        url_endpoint = "{}.dssp".format(identifier)
        url = url_root + url_endpoint
        Downloader(url=url, filename=filename,
                   decompress=decompress, overwrite=overwrite)


# TODO fix for speed
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


# TODO fix for speed
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


# TODO fix for speed
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

    return rsa


def get_rsa_class(rsa):
    """
    Gets a class based on the RSA value

    :param rsa: (float) RSA score or string
    :return: RSA class
    """
    rsa_class = '-'
    try:
        rsa = float(rsa)
        # surface is rsa >= 25% (value = 2)
        # exposed is rsa >= 5% and rsa < 25% (value = 1)
        # core is rsa < 5% (value = 0)
        if rsa >= 25.0:
            rsa_class = 'Surface'
        elif 5.0 <= rsa < 25.0:
            rsa_class = 'Part. Exposed'
        else:
            rsa_class = 'Core'
    except ValueError:
        # returns a string
        pass
    return rsa_class
