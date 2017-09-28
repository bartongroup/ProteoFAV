# -*- coding: utf-8 -*-

import os
import gzip
import json
import copy
import time
import shutil
import socket
import logging

try:
    from urllib import quote as urllib_quote
    from urllib import urlretrieve
except ImportError:  # python 3.5
    from urllib.parse import quote as urllib_quote
    from urllib.request import urlretrieve

import numpy as np
import pandas as pd
import requests

from proteofav.config import defaults
from proteofav.library import valid_ensembl_species_variation

__all__ = ["get_url_or_retry", "check_local_or_fetch", "fetch_files", "get_preferred_assembly_id",
           "IDNotValidError", "raise_if_not_ok", "_pdb_uniprot_sifts_mapping",
           "_uniprot_pdb_sifts_mapping", "icgc_missense_variant", "is_valid",
           "is_valid_ensembl_id", "confirm_column_types"]

log = logging.getLogger('proteofav')
socket.setdefaulttimeout(15)  # avoid infinite hanging


def get_url_or_retry(url, retry_in=None, wait=1, n_retries=10, json=False, header=None, **params):
    """
    Fetch an url using Requests or retry fetching it if the server is
        complaining with retry_in error.

    :param retry_in: list or array of http status codes
    :param json: boolean
    :param header: dictionary with head params
    :param url: url to be fetched as a string
    :param wait: sleeping between tries in seconds
    :param n_retries: number of retry attempts
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
    elif response.status_code in retry_in and n_retries >= 0:
        time.sleep(wait)
        return get_url_or_retry(url, retry_in, (n_retries - 1), wait, json, header, **params)
    else:
        log.error(response.status_code)
        response.raise_for_status()


def fetch_files(identifier, directory=None, sources=("cif", "dssp", "sifts")):
    """
    Small routine to fetch data files from their respective repositories.
        Use defaults to fetch files from servers. Defaults server are defined
        in the default config.txt. Returns None since it has side effects.

    :param identifier: protein identifier as PDB identifier
    :param directory: path to download. The default downloads all files to
        default.temp folder. If is a string and valid path downloads all files
        to that folder. If its a iterable and all valid path downloads to the
        respective folders

    :param sources: where to fetch the data. Must be in the config file.
    :return list: list of paths
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
            urlretrieve(url, destination + filename)
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

    :param str ensembl_gene_id: ensembl gene accession
    :return pd.DataFrame: DataFrame with one mutation per row.
    :example:

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
            url = defaults.api_ensembl + ensembl_endpoint + urllib_quote(identifier, safe='')
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

    :param pd.DataFrame table: A table produced by ProteoFAV
    :return pd.DataFrame: the table with the same data, but normalised column types

    .. note:: There are fewer pandas `dtypes` than the corresponding numpy
        type, but all numpy types can be accommodated.

    NaNs and Upcasting:
    The upcasting of dtypes on columns that contain NaNs has been an issue and
        can lead to inconsistencies between the column types of different
        tables that contain equivalent data. E.g., a `merged_table` from a PDB
        entry may have NaNs in UniProt_dbResNum and so is cast to float64
        whilst a UniProt variants table will have no NaNs in this field and
        so remains int64 by default. The main issue here is that it creates
        additional work for merging these tables as the keys must be of the
        same type. It's also a problem if you want to test elements say for
        example you expect to be testing integer equality but the values have
        been coerced to floats.

    Ideally then we need to ensure that all equivalent column types are by
        default the least generic dtype that can contain all possible values
        that can be seen in the column, i.e. if NaNs are possible then even
        when not present, a column of integers that can contain NaNs should
        always be at least float.

    dtypes and element types:
    It seems that coercion to certain dtypes alters the element types too.
        For instance, int64 -> float64 will make an equivalent change to
        individual value types. However, the original element types can be
        preserved with the generic `object` dtype, so this will be used in
        preference for integer columns that can contain NaNs.

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


def row_selector(data, key=None, value=None, reverse=False):
    """
    Generic method to filter columns
    :param data: pandas DataFrame
    :param key: pandas DataFrame column name
    :param value: value(s) to be looked for
    :param reverse: opposite behavior (e.g. 'isin' becomes 'isnotin'
        and 'equals' becomes 'differs')
    :return: returns a modified pandas DataFrame
    """

    table = data
    assert type(table) is pd.core.frame.DataFrame
    if key is not None and value is not None:
        assert type(key) is str
        if key in table:
            if value == 'first':
                value = table[key].iloc[0]
                table = table.loc[table[key] == value]
            elif (hasattr(value, '__iter__') and
                      (type(value) is tuple or type(value) is list)):
                if not reverse:
                    table = table.loc[table[key].isin(value)]
                else:
                    table = table.loc[~table[key].isin(value)]
            else:
                if not reverse:
                    table = table.loc[table[key] == value]
                else:
                    table = table.loc[table[key] != value]
        else:
            log.debug("%s not in the DataFrame...", key)

    if table.empty:
        message = 'Column {} does not contain {} value(s)...'.format(key, value)
        log.debug(message)
        raise ValueError(message)

    return table


def constrain_column_types(data, col_type_dict=None, nan_value_dict=None,
                           replace_value_dict=None):
    """
    Helper method that helps in constrain data types for the
    various DataFrame columns.

    :param data: pandas DataFrame
    :param col_type_dict: (dict) optional defines common types
    :param nan_value_dict: (dict) optional new value passed to replace NaNs
    :param replace_value_dict: (dict) optional new value passed to replace
        specific values
    :return: modified pandas DataFrame
    """

    table = data
    for col in table:
        if col_type_dict is not None and col in col_type_dict:
            try:
                table[col] = table[col].astype(col_type_dict[col])
            except (ValueError, KeyError):
                # probably there are some NaNs in there
                pass
        if nan_value_dict is not None and col in nan_value_dict:
            if table[col].isnull().any().any():
                table[col] = table[col].fillna(nan_value_dict[col])
        if replace_value_dict is not None and col in replace_value_dict:
            table[col] = table[col].replace(replace_value_dict[col][0],
                                            replace_value_dict[col][1])

    return table


def exclude_columns(data, excluded=()):
    """
    Helper method that helps in filtering out columns based
    on the column name.

    :param data: pandas DataFrame
    :param excluded: (tuple) optional columns to be excluded
    :return: modified pandas DataFrame
    """

    table = data
    if excluded is not None:
        assert hasattr(excluded, '__iter__')
        try:
            table = table.drop(list(excluded), axis=1)
        except ValueError:
            # most likely theses are not in there
            pass
    return table


def merging_down_by_key(table, key="xrefs_id"):
    """
    Helper method that mergers the rows  containing data
    that have the same value (e.g. ID), according
    to the 'key' passed to function.
    This works as a collapse down (from many-to-one rows).
    Aggregation is possible since multi-value cells are stored as
    tuples which are hashable.

    :param table: pandas DataFrame
    :param key: key to base the 'merge down' upon
    :return: modified pandas DataFrame
    """

    new_table = table.copy()
    duplicated = {}
    drop_indexes = []
    for ix in table.index:
        pid = table.loc[ix, key]
        dup = table[table[key] == pid].index.tolist()
        duplicated[pid] = dup
        if len(dup) > 1:
            drop_indexes.append(ix)

    new_table = new_table.drop(new_table.index[drop_indexes])
    for key, val in duplicated.items():
        if type(val) is list and len(val) > 1:
            rows = []
            d = {}
            for i in range(len(list(table))):
                values = []
                for j in range(len(val)):
                    v = table.loc[val[j], list(table)[i]]
                    if v not in values:
                        if type(v) is tuple or type(v) is list:
                            for g in v:
                                values.append(g)
                        else:
                            values.append(v)
                if not values:
                    values = np.nan
                elif len(values) == 1:
                    values = values[0]
                else:
                    values = [v for v in values if not pd.isnull(v)]
                    if not values:
                        values = np.nan
                    elif len(values) == 1:
                        values = values[0]
                    else:
                        values = tuple(set(values))

                d[list(table)[i]] = values
            rows.append(d)
            combined = pd.DataFrame(rows)
            new_table = new_table.append(combined)

    return new_table.reset_index(drop=True)


def splitting_up_by_key(table, key="xrefs_id"):
    """
    Helper method that splits the rows containing
    multi-value entries (e.g. [ID_1, ID_2]), according
    to the 'key' passed to the function.

    :param table: pandas DataFrame
    :param key: key to base the 'split up' upon
    :return: modified pandas DataFrame
    """
    rows = []
    for ix in table.index:
        entries = {k: table.loc[ix, k] for k in list(table)}
        val = table.loc[ix, key]
        if type(val) is tuple or type(val) is list:
            for v in val:
                nentries = copy.deepcopy(entries)
                nentries[key] = v
                rows.append(nentries)
        else:
            rows.append(entries)

    new_table = pd.DataFrame(rows)

    return new_table.reset_index(drop=True)


def flatten_nested_structure(data, dictionary, keys=None, values=None):
    """
    Flattens a deeply nested json structure to unique columns (keys),
    where columns with multiple values are aggregated to the same column,
    and the items grouped by order in a list. Also copes with

    :param data: json.load() content
    :param dictionary: mutable dictionary
    :param keys: dict keys (used during recursion)
    :param values: dict values (used during recursion)
    :return: (side-effects) updates an input dictionary
    """
    if type(data) is tuple or type(data) is list:
        for e in data:
            flatten_nested_structure(e, dictionary, keys, values)
    elif type(data) is dict:
        for k, v in data.items():
            if keys is not None:
                k = '_'.join([keys, k])
            flatten_nested_structure(v, dictionary, k, v)
    else:
        if keys is not None and values is not None:
            if keys not in dictionary:
                dictionary[keys] = [values]
            else:
                if values not in dictionary[keys]:
                    dictionary[keys].append(values)


def refactor_key_val_singletons(dictionary):
    """
    Simply updates a dictionary values that are singletons
    i.e. len(val) == 1:

    :param dictionary: mutable dictionary
    :return: updated dictionary
    """
    new_dictionary = {}
    for k, v in dictionary.items():
        if (type(v) is list or type(v) is tuple) and len(v) == 1:
            new_dictionary[k] = v[0]
        else:
            new_dictionary[k] = v
    return new_dictionary


class InputFileHandler(object):
    """Validates input file paths."""

    def __init__(self, filename):
        self.__filename = filename
        self.__validate()

    def __validate(self):
        if not os.path.isfile(self.__filename):
            raise IOError("File '%s' not available..." % self.__filename)


class OutputFileHandler(object):
    """Validates output file paths."""

    def __init__(self, filename, overwrite=False):
        self.__filename = filename
        self.__overwrite = overwrite
        self.__validate()

    def __validate(self):
        if os.path.exists(self.__filename) and not self.__overwrite:
            raise OSError("File '%s' already exists..." % self.__filename)
        elif not os.access(os.path.dirname(self.__filename), os.W_OK):
            raise OSError("File '%s' cannot be written..." % self.__filename)


class GenericInput(object):
    def __init__(self, identifier=None, filename=None, table=None):
        self.identifier = identifier
        self.filename = filename
        self.table = table

    def _get_identifier(self, identifier):
        if identifier is None and self.identifier is None:
            raise ValueError("Input identifier needed!")
        elif identifier is None:
            identifier = self.identifier
        return identifier

    def _get_filename(self, filename):
        if filename is None and self.filename is None:
            raise ValueError("Input/Output filename needed!")
        elif filename is None:
            filename = self.filename
        return filename

    def _get_table(self, table):
        if table is None and self.table is None:
            raise ValueError("Pandas DataFrame needed!")
        elif table is None:
            table = self.table
        return table


class TableMergerError(Exception):
    pass


def fetch_from_url_or_retry(url, json=True, headers=None, post=False, data=None,
                            retry_in=None, wait=1, n_retries=10, stream=False, **params):
    """
    Fetch an url using Requests or retry fetching it if the server is
    complaining with retry_in error. There is a limit to the number of retries.

    Retry code examples: 429, 500 and 503

    :param url: url to be fetched as a string
    :param json: json output
    :param headers: dictionary
    :param post: boolean
    :param data: dictionary: only if post is True
    :param retry_in: http codes for retrying
    :param wait: sleeping between tries in seconds
    :param n_retries: number of retry attempts
    :param stream: boolean
    :param params: request.get kwargs.
    :return: url content
    """

    if retry_in is None:
        retry_in = ()
    else:
        assert type(retry_in) is tuple or type(retry_in) is list

    if headers is None:
        headers = {}
    else:
        assert type(headers) is dict

    if json:
        headers.update({"Content-Type": "application/json"})
    else:
        if "Content-Type" not in headers:
            headers.update({"Content-Type": "text/plain"})

    if post:
        if data is not None:
            assert type(data) is dict or type(data) is str
            response = requests.post(url, headers=headers, data=data, params=params)
        else:
            return None
    else:
        response = requests.get(url, headers=headers, params=params, stream=stream)

    if response.ok:
        return response
    elif response.status_code in retry_in and n_retries >= 0:
        time.sleep(wait)
        return fetch_from_url_or_retry(url, json, headers, post, data, retry_in, wait,
                                       (n_retries - 1), stream, **params)
    else:
        try:
            response.raise_for_status()
        except requests.exceptions.HTTPError as e:
            log.debug("%s: Unable to retrieve %s for %s",
                      response.status_code, url, e)


if __name__ == '__main__':
    pass
