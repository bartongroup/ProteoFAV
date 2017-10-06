# -*- coding: utf-8 -*-

import os
import gzip
import json
import time
import shutil
import socket
import logging
import tempfile
import requests
import pandas as pd
import requests_cache

try:
    from urllib import quote as urllib_quote
    from urllib import urlretrieve
except ImportError:  # python 3.5
    from urllib.parse import quote as urllib_quote
    from urllib.request import urlretrieve

from proteofav.library import valid_ensembl_species_variation

from proteofav.config import defaults

log = logging.getLogger('proteofav.config')
socket.setdefaulttimeout(15)  # avoid infinite hanging

requests_cache.install_cache('proteofav')

__all__ = ["fetch_from_url_or_retry", "check_local_or_fetch", "fetch_files", "get_preferred_assembly_id",
           "IDNotValidError", "raise_if_not_ok", "_pdb_uniprot_sifts_mapping",
           "_uniprot_pdb_sifts_mapping", "icgc_missense_variant", "is_valid",
           "is_valid_ensembl_id", "confirm_column_types"]


def fetch_from_url_or_retry(url, json=True, header=None, post=False, data=None,
                            retry_in=None, wait=1, n_retries=10, stream=False, **params):
    """
    Fetch an url using Requests or retry fetching it if the server is
    complaining with retry_in error. There is a limit to the number of retries.

    Retry code examples: 429, 500 and 503

    :param url: url to be fetched as a string
    :param json: json output
    :param header: dictionary
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

    if header is None:
        header = {}
    else:
        assert type(header) is dict

    if json:
        header.update({"Content-Type": "application/json"})
    else:
        if "Content-Type" not in header:
            header.update({"Content-Type": "text/plain"})

    if post:
        if data is not None:
            assert type(data) is dict or type(data) is str
            response = requests.post(url, headers=header, data=data)
        else:
            return None
    else:
        response = requests.get(url, headers=header, params=params, stream=stream)

    if response.ok:
        return response
    elif response.status_code in retry_in and n_retries >= 0:
        time.sleep(wait)
        return fetch_from_url_or_retry(url, json, header, post, data, retry_in, wait,
                                       (n_retries - 1), stream, **params)
    else:
        try:
            response.raise_for_status()
        except requests.exceptions.HTTPError as e:
            log.debug("%s: Unable to retrieve %s for %s",
                      response.status_code, url, e)


def fetch_files(identifier, directory=None, sources=("cif", "dssp", "sifts"),
                overwrite=False):
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
    :param overwrite: (boolean) Overrides any existing file, if available
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
        if filename.endswith('.gz'):
            filename = destination + filename.replace('.gz', '')
            Downloader(url=url, filename=filename, overwrite=overwrite,
                       decompress=True)
        else:
            filename = destination + filename
            Downloader(url=url, filename=filename, overwrite=overwrite,
                       decompress=False)
        result.append(filename)
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
    rows = fetch_from_url_or_retry(url, retry_in=retry_in, json=True).json()
    return rows


def get_preferred_assembly_id(identifier, verbose=False):
    """
    Gets the preferred assembly id for the given PDB ID, from the PDBe API.

    :param identifier: PDB ID
    :param verbose: boolean
    :return: str
    """

    # getting the preferred biological assembly from the PDBe API
    try:
        data = _fetch_summary_properties_pdbe(identifier)
    except Exception as e:
        message = "Something went wrong for {}... {}".format(identifier, e)
        if verbose:
            log.error(message)
    try:
        nassemblies = data[identifier][0]["assemblies"]
        if len(nassemblies) > 1:
            for entry in nassemblies:
                if entry["preferred"]:
                    pref_assembly = entry["assembly_id"]
                    break
        else:
            pref_assembly = data[identifier][0]["assemblies"][0]["assembly_id"]
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
    information = fetch_from_url_or_retry(url, json=True).json()

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
    information = fetch_from_url_or_retry(url, json=True).json()

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
            data = fetch_from_url_or_retry(url, json=True).json()
            if 'error' in data:
                return False
            return True
        else:
            raise IDNotValidError
    except IDNotValidError:
        return False
    except requests.HTTPError:
        return False


def row_selector(table, key=None, value=None, reverse=False):
    """
    Generic method to filter columns
    :param table: pandas DataFrame
    :param key: pandas DataFrame column name
    :param value: value(s) to be looked for
    :param reverse: opposite behavior (e.g. 'isin' becomes 'isnotin'
        and 'equals' becomes 'differs')
    :return: returns a modified pandas DataFrame
    """

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


def constrain_column_types(table, col_type_dict=None, nan_value_dict=None,
                           replace_value_dict=None):
    """
    Helper method that helps in constraining data types for the
    various DataFrame columns.

    This function checks a table's column types against a defined
    column name/type dictionary and correct them if necessary.

    .. notes::

    There are fewer pandas `dtypes` than the corresponding numpy
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

    :param table: pandas DataFrame
    :param col_type_dict: (dict) optional defines common types
    :param nan_value_dict: (dict) optional new value passed to replace NaNs
    :param replace_value_dict: (dict) optional new value passed to replace
        specific values
    :return: modified pandas DataFrame
        (the table with the same data, but normalised column types)
    """

    for col in table:
        if col_type_dict is not None and col in col_type_dict:
            try:
                table[col] = table[col].astype(col_type_dict[col])
            except (ValueError, KeyError, TypeError):
                # probably there are some NaNs in there
                # and it is impossible to coerce the column to the
                # pre-defined dtype
                pass
        if nan_value_dict is not None and col in nan_value_dict:
            if table[col].isnull().any().any():
                table[col] = table[col].fillna(nan_value_dict[col])
        if replace_value_dict is not None and col in replace_value_dict:
            table[col] = table[col].replace(replace_value_dict[col][0],
                                            replace_value_dict[col][1])

    return table


def exclude_columns(table, excluded=()):
    """
    Helper method that helps in filtering out columns based
    on the column name.

    :param table: pandas DataFrame
    :param excluded: (tuple) optional columns to be excluded
    :return: modified pandas DataFrame
    """

    if excluded is not None:
        assert hasattr(excluded, '__iter__')
        try:
            table = table.drop(list(excluded), axis=1)
        except ValueError:
            # most likely theses are not in there
            pass
    return table


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


class Downloader(object):
    def __init__(self, url, filename, decompress=False, overwrite=False):
        """
        :param filename: (str) Output filename
        :param decompress: (boolean) Decompresses the file
        :param overwrite: (boolean) Overrides any existing file, if available
        """

        self._url = url
        self._filename = filename
        self._decompress = decompress
        self._overwrite = overwrite

        if not os.path.exists(self._filename) or self._overwrite:
            OutputFileHandler(self._filename, overwrite=self._overwrite)
            self._tempfile = tempfile.NamedTemporaryFile().name
            OutputFileHandler(self._tempfile, overwrite=self._overwrite)

            self._download()
            if self._decompress:
                self._uncompress()

    def _download(self):
        try:
            try:
                import urllib.request
                from urllib.error import URLError, HTTPError
                with urllib.request.urlopen(self._url) as response, \
                        open(self._tempfile, 'wb') as outfile:
                    shutil.copyfileobj(response, outfile)
            except (AttributeError, ImportError):
                import urllib
                urllib.urlretrieve(self._url, self._tempfile)
            InputFileHandler(self._tempfile)
            with open(self._tempfile, 'rb') as infile, \
                    open(self._filename, 'wb') as outfile:
                shutil.copyfileobj(infile, outfile)
                os.remove(self._tempfile)
        except (URLError, HTTPError, IOError, Exception) as e:
            log.debug("Unable to retrieve %s for %s", self._url, e)

    def _uncompress(self):
        InputFileHandler(self._filename)
        with open(self._filename, 'rb') as infile, \
                open(self._tempfile, 'wb') as outfile:
            shutil.copyfileobj(infile, outfile)
        InputFileHandler(self._tempfile)
        with gzip.open(self._tempfile, 'rb') as infile, \
                open(self._filename, 'wb') as outfile:
            shutil.copyfileobj(infile, outfile)
            os.remove(self._tempfile)
            log.info("Decompressed %s", self._filename)


class GenericInputs(object):
    def __init__(self, identifier=None, filename=None, table=None):
        self._identifier = identifier
        self._filename = filename
        self._table = table

    def _get_identifier(self, identifier=None):
        if identifier is None and self._identifier is None:
            raise ValueError("An Identifier is needed!")
        elif identifier is not None:
            self._identifier = identifier
        return self._identifier

    def _get_filename(self, filename=None):
        if filename is None and self._filename is None:
            raise ValueError("A filename is needed!")
        elif filename is not None:
            self._filename = filename
        return self._filename

    def _get_table(self, table=None):
        if table is None and self._table is None:
            raise ValueError("A Pandas DataFrame is needed!")
        elif table is not None:
            self._table = table
        return self._table


if __name__ == '__main__':
    pass
