#!/usr/bin/env python
# -*- coding: utf-8 -*-

import gzip
import json
import logging
import os
import shutil
import socket
import time
import urllib

import numpy as np
import pandas as pd
import requests

from .config import defaults
from .library import valid_ensembl_species_variation

__all__ = ["get_url_or_retry", "check_local_or_fetch", "fetch_files", "get_preferred_assembly_id",
           "IDNotValidError", "raise_if_not_ok", "_pdb_uniprot_sifts_mapping",
           "_uniprot_pdb_sifts_mapping", "icgc_missense_variant", "is_valid",
           "is_valid_ensembl_id", "confirm_column_types"]

socket.setdefaulttimeout(15)
log = logging.getLogger('proteofav.config')


def get_url_or_retry(url, retry_in=None, wait=1, json=False, header=None, **params):
    """
    Fetch an url using Requests or retry fetching it if the server is
        complaining with retry_in error.

    :param retry_in: list or array of http status codes
    :param json: boolean
    :param header: dictionary with head params
    :param url: url to be fetched as a string
    :param wait: sleeping between tries in seconds
    :param params: request.get kwargs.
    :return: url content or url content in json data structure.
    """
    if not header:
        header = {}
    if not retry_in:
        retry_in = []
    if json:
        header.update({"Content-Type": "application/json"})
    response = requests.get(url, headers=header, params=params)

    if response.ok and json:
        return response.json()
    elif response.ok:
        return response.content
    elif response.status_code in retry_in:
        time.sleep(wait)
        return get_url_or_retry(url, retry_in, wait, json, header, **params)
    else:
        log.error(response.status_code)
        response.raise_for_status()


def check_local_or_fetch(identifier, update=False):
    """

    :param update:
    :param identifier:
    :return:
    """
    # TODO fetch what?
    pass


def fetch_files(identifier, directory=None, sources=("cif", "dssp", "sifts")):
    """
    Small routine to fetch data files from their respective repositories.

    Use defaults to fetch files from servers. Defaults server are defined in
        the default config.txt. Returns None since it has side effects.

    Created on 11:41 24/07/15 2015
    :param identifier: protein identifier as PDB identifier
    :param directory: path to download. The default downloads all files to
        default.temp folder. If is a string and valid path downloads all files
        to that folder. If its a iterable and all valid path downloads to the
        respective folders
    :param sources: where to fetch the data. Must be in the config file.
    :return: list of path
    :rtype: list
    :raise: TypeError
    """
    if isinstance(sources, str):
        sources = [sources]
    if not directory:
        directory = defaults.tmp
    elif isinstance(directory, str):
        if not os.path.isdir(directory):
            raise IOError(directory + " is not a valid path.")
        directory = [directory] * len(sources)
    elif hasattr(directory, "__iter__"):
        if len(directory) != len(sources):
            raise IOError(directory + " you need one directory for each source,"
                                      "or a path for all.")
        for d in directory:
            if not os.path.isdir(d):
                raise IOError(d + " is not a valid path.")
    else:
        raise TypeError('Unidentified source|directory combination.')

    result = []
    for source, destination in zip(sources, directory):
        filename = identifier + getattr(defaults, source + '_extension')
        url = getattr(defaults, source + '_fetch') + filename
        try:
            try:
                urllib.urlretrieve(url, destination + filename)
            except AttributeError:
                # python 3.5
                urllib.request.urlretrieve(url, destination + filename)
        except IOError as e:
            log.error('Unable to retrieve {} for {}'.format(url, str(e)))
            raise
        if filename.endswith('.gz'):
            with gzip.open(destination + filename, 'rb') as input_f, \
                    open(destination + filename.replace('.gz', ''),
                         'wb') as output_f:
                shutil.copyfileobj(input_f, output_f)
                filename = filename.replace('.gz', '')
        result.append(destination + filename)
    return result


def _fetch_summary_properties_pdbe(pdbid, retry_in=(429,)):
    """
    Queries the PDBe api to get summary validation report.

    :param pdbid: PDB ID
    :param retry_in: http code for retrying connections
    :return: dictionary
    """
    pdbe_endpoint = "pdb/entry/summary/"
    url = defaults.api_pdbe + pdbe_endpoint + pdbid
    rows = get_url_or_retry(url, retry_in=retry_in, json=True)
    return rows


def get_preferred_assembly_id(pdbid, verbose=False):
    """
    Gets the preferred assembly id for the given PDB ID, from the PDBe API.

    :param pdbid: PDB ID
    :param verbose: boolean
    :return: str
    """

    # getting the preferred biological assembly from the PDBe API
    try:
        data = _fetch_summary_properties_pdbe(pdbid)
    except Exception as e:
        message = "Something went wrong for {}... {}".format(pdbid, e)
        if verbose:
            log.error(message)
    try:
        nassemblies = data[pdbid][0]["assemblies"]
        if len(nassemblies) > 1:
            for entry in nassemblies:
                if entry["preferred"]:
                    pref_assembly = entry["assembly_id"]
                    break
        else:
            pref_assembly = data[pdbid][0]["assemblies"][0]["assembly_id"]
    except Exception as e:
        pref_assembly = "1"

    bio_best = str(pref_assembly)
    return bio_best


##############################################################################
# Some other utils not currently in use
##############################################################################
class IDNotValidError(Exception):
    """
    Base class for database related exceptions.
    Databases: UniProt, PDB, Ensembl, etc.
    """
    pass


def raise_if_not_ok(response):
    """

    :param response:
    """
    if not response.ok:
        response.raise_for_status()


def _pdb_uniprot_sifts_mapping(identifier):
    """
    Queries the PDBe API for SIFTS mapping between PDB - UniProt.
    One to many relationship expected.

    :param identifier: PDB id
    :return: pandas table dataframe
    """

    sifts_endpoint = 'mappings/uniprot/'
    url = defaults.api_pdbe + sifts_endpoint + identifier
    information = get_url_or_retry(url, json=True)

    rows = []
    for uniprot in information[identifier]['UniProt']:
        uniprots = {'uniprot_id': uniprot}
        rows.append(uniprots)
    return pd.DataFrame(rows)


def _uniprot_pdb_sifts_mapping(identifier):
    """
    Queries the PDBe API for SIFTS mapping between UniProt - PDB entries.
    One to many relationship expected.

    :param identifier: UniProt ID
    :return: pandas table dataframe
    """
    sifts_endpoint = 'mappings/best_structures/'
    url = defaults.api_pdbe + sifts_endpoint + str(identifier)
    information = get_url_or_retry(url, json=True)

    rows = []
    for entry in information[identifier]:
        rows.append(entry)
    return pd.DataFrame(rows)


def icgc_missense_variant(ensembl_gene_id):
    """Fetch a gene missense variants from ICGC data portal.

    :param ensembl_gene_id: ensembl gene accession
    :type ensembl_gene_id: str
    :return: DataFrame with one mutation per row.
    :rtype : pandas.DataFrame

    :Example:


        >>> table = icgc_missense_variant('ENSG00000012048')
        >>> table.loc[0, ['id', 'mutation', 'type', 'start']]
        id                          MU601299
        mutation                         G>A
        type        single base substitution
        start                       41245683
        Name: 0, dtype: object


    .. note:: ICGC doc https://dcc.icgc.org/docs/#!/genes/findMutations_get_8

    """
    base_url = defaults.api_icgc + ensembl_gene_id + "/mutations/"
    headers = {'content-type': 'application/json'}
    filt = json.dumps({"mutation": {"consequenceType": {"is": ['missense_variant']}}})
    params = {"filters": filt}

    # First counts the number of mutation entries
    counts_resp = requests.get(base_url + "counts/", headers=headers, params=params)
    raise_if_not_ok(counts_resp)
    total = counts_resp.json()['Total']

    # then iterate the pages of entries, max 100 entries per page
    hits = []
    params['size'] = 100
    for i in range(total // 100 + 1):
        params['from'] = i * 100 + 1
        mutation_resp = requests.get(base_url, headers=headers, params=params)
        raise_if_not_ok(mutation_resp)
        hits.extend(mutation_resp.json()['hits'])

    return pd.DataFrame(hits)


def is_valid(identifier, database=None, url=None):
    """
    Generic method to check if a given id is valid.

    :param url: if given use this instead
    :param identifier: accession id
    :param database: database to test against
    :return: simply a true or false
    :rtype: boolean
    """

    try:
        identifier = str(identifier)
    except ValueError:
        # raise IDNotValidError
        return False

    if len(str(identifier)) < 1:
        # raise IDNotValidError
        return False
    if not url:
        url = getattr(defaults, 'api_' + database) + identifier
    r = requests.get(url)
    if not r.ok:
        # not the best approach
        try:
            raise IDNotValidError('{} not found at {}: (url check:{}'.format(
                identifier, database, r.url))
        except IDNotValidError:
            return False
    else:
        return True


def is_valid_ensembl_id(identifier, species='human', variant=False):
    """
    Checks if an Ensembl id is valid.

    :param identifier: testing ID
    :param species: Ensembl species
    :param variant: boolean if True uses the variant endpoint
    :return: simply a true or false
    :rtype: boolean
    """

    try:
        identifier = str(identifier)
    except ValueError:
        # raise IDNotValidError
        return False

    if len(str(identifier)) < 1:
        # raise IDNotValidError
        return False

    if variant:
        if species not in valid_ensembl_species_variation:
            raise ValueError('Provided species {} is not valid'.format(species))
        ensembl_endpoint = "variation/{}/".format(species)
    else:
        ensembl_endpoint = "lookup/id/"
    try:
        if identifier != '':
            try:
                url = defaults.api_ensembl + ensembl_endpoint + urllib.quote(identifier, safe='')
            except AttributeError:
                # python 3.5
                url = defaults.api_ensembl + ensembl_endpoint + urllib.parse.quote(identifier,
                                                                                   safe='')
            data = get_url_or_retry(url, json=True)
            if 'error' in data:
                return False
            return True
        else:
            raise IDNotValidError
    except IDNotValidError:
        return False
    except requests.HTTPError:
        return False


def confirm_column_types(table):
    """
    Check a table's column types against a defined column name/type dictionary
    and correct them if necessary.

    :param table: A pandas data frame produced by a to_* function
    :return: A pandas data frame of the same data with correct column types

    .. note::
    There are fewer pandas `dtypes` than the corresponding numpy type classes, but all numpy types
    can be accommodated.

    NaNs and Upcasting:
    The upcasting of dtypes on columns that contain NaNs has been an issue and can lead to
    inconsistencies between the column types of different tables that contain equivalent data.
    E.g., a `merged_table` from a PDB entry may have NaNs in UniProt_dbResNum and so is cast to
    float64 whilst a UniProt variants table will have no NaNs in this field and so remains int64
    by default. The main issue here is that it creates additional work for merging these tables as
    the keys must be of the same type. It's also a problem if you want to test elements say for
    example you expect to be testing integer equality but the values have been coerced to floats.

    Ideally then we need to ensure that all equivalent column types are by default the least
    generic dtype that can contain all possible values that can be seen in the column, i.e. if
    NaNs are possible then even when not present, a column of integers that can contain NaNs
    should always be at least float.

    dtypes and element types:
    It seems that coercion to certain dtypes alters the element types too. For instance, int64 ->
    float64 will make an equivalent change to individual value types. However, the original element
    types can be preserved with the generic `object` dtype, so this will be used in preference for
    integer columns that can contain NaNs.
    """
    # TODO: update this to accomodate the new sifts table headers
    column_types_long = {
        # SIFTs mappings
        'CATH_dbAccessionId': 'string',
        'CATH_dbChainId': 'string',
        'CATH_dbCoordSys': 'string',
        'CATH_dbResName': 'string',
        'CATH_dbResNum': 'string',
        'InterPro_dbAccessionId': 'string',
        'InterPro_dbCoordSys': 'string',
        'InterPro_dbEvidence': 'string',
        'InterPro_dbResName': 'string',
        'InterPro_dbResNum': 'string',
        'NCBI_dbAccessionId': 'string',
        'NCBI_dbCoordSys': 'string',
        'NCBI_dbResName': 'string',
        'NCBI_dbResNum': 'string',
        'PDB_dbAccessionId': 'string',
        'PDB_dbChainId': 'string',
        'PDB_dbCoordSys': 'string',
        'PDB_dbResName': 'string',
        'PDB_dbResNum': 'string',
        'Pfam_dbAccessionId': 'string',
        'Pfam_dbCoordSys': 'string',
        'Pfam_dbResName': 'string',
        'Pfam_dbResNum': 'string',
        'REF_codeSecondaryStructure': 'string',
        'REF_dbCoordSys': 'string',
        'REF_dbResName': 'string',
        'REF_dbResNum': 'string',
        'REF_nameSecondaryStructure': 'string',
        'UniProt_dbAccessionId': 'string',
        'UniProt_dbCoordSys': 'string',
        'UniProt_dbResName': 'string',
        'UniProt_dbResNum': 'string',
        # mmCIF fields
        'auth_asym_id': 'string',
        'auth_atom_id': 'string',
        'auth_comp_id': 'string',
        'auth_seq_id': 'string',
        'B_iso_or_equiv': 'float',
        'B_iso_or_equiv_esd': 'float',
        'Cartn_x': 'float',
        'Cartn_x_esd': 'float',
        'Cartn_y': 'float',
        'Cartn_y_esd': 'float',
        'Cartn_z': 'float',
        'Cartn_z_esd': 'float',
        'label_alt_id': 'string',
        'label_asym_id': 'string',
        'label_atom_id': 'string',
        'label_comp_id': 'string',
        'label_entity_id': 'string',
        'label_seq_id': 'integer',
        'occupancy': 'float',
        'occupancy_esd': 'float',
        'pdbe_label_seq_id': 'integer',
        'pdbx_formal_charge': 'integer',
        'pdbx_PDB_ins_code': 'string',
        'pdbx_PDB_model_num': 'integer',
        'type_symbol': 'string',
        'group_PDB': 'string',
        'id': 'string',
        # DSSP
        'chain_id': 'string',
        'aa': 'string',
        'ss': 'string',
        'acc': 'float',
        'phi': 'float',
        'psi': 'float',
        # Merged table
        'dssp_aa': 'string',
        'cif_aa': 'string',
        'sifts_aa': 'string',
        # UniProt variant table
        'resn': 'string',
        'mut': 'string',
        'disease': 'string',
        # UniProt variants table 2
        'translation': 'string',
        'id': 'string',
        'start': 'string',
        'residues': 'string',
        # Derived boolean columns
        'is_expression_tag': 'bool',
        'is_not_observed': 'bool'
    }

    type_to_dtype = {
        'string': 'object',
        'float': 'float64',
        'integer': 'int64',
        'bool': 'bool',
        'O': 'object'
    }

    type_to_dtype_if_contains_nan = {
        'string': 'object',
        'float': 'float64',
        'integer': 'object',
        'bool': 'object',
        'O': 'object'
    }

    # Columns in this dictionary will undergo `.replace(value[0], value[1])`
    column_replacements = {
        'Cartn_x_esd': ['?', np.nan],
        'Cartn_y_esd': ['?', np.nan],
        'Cartn_z_esd': ['?', np.nan],
        'occupancy_esd': ['?', np.nan],
        'B_iso_or_equiv_esd': ['?', np.nan]
    }

    for column in table:
        type_should_be = column_types_long.get(column)

        # Get dtype from column type depending on whether can contain NaN
        can_be_nan = True  # TODO: Either this should be a test or just get rid of the if block
        if can_be_nan:
            dtype_should_be = type_to_dtype_if_contains_nan.get(type_should_be)
        else:
            dtype_should_be = type_to_dtype.get(type_should_be)

        if dtype_should_be is None:
            logging.warning('Column `{}` not recognised'.format(column))
            continue

        # Element replacements as required
        if column in column_replacements:
            to_replace, replacement = column_replacements.get(column)
            logging.debug(
                'Replacing {} with {} in column {}'.format(to_replace, replacement, column))
            table[column] = table[column].replace(to_replace, replacement)

        # Coerce column if neccessary
        current_dtype = table[column].dtype
        if current_dtype != dtype_should_be:
            logging.debug('Coercing `{}` to `{}`'.format(column, dtype_should_be))
            table[column] = table[column].astype(dtype_should_be)

    # Now check that the index is the correct type
    column = table.index.name
    type_should_be = column_types_long.get(column)
    dtype_should_be = type_to_dtype_if_contains_nan.get(type_should_be)
    if dtype_should_be is None:
        logging.warning('Index column `{}` not recognised'.format(column))
    current_dtype = table.index.dtype
    if current_dtype != dtype_should_be:
        logging.debug('Coercing index `{}` to `{}`'.format(column, dtype_should_be))
        table.index = table.index.astype(dtype_should_be)
    return table


if __name__ == '__main__':
    pass
