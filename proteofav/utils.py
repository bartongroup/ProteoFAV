# -*- coding: utf-8 -*-

import os
import gzip
import time
import copy
import shutil
import socket
import logging
import tempfile
import requests
import numpy as np
import pandas as pd
import requests_cache

try:
    from urllib import quote as urllib_quote
    from urllib import urlretrieve
except ImportError:  # python 3.5
    from urllib.parse import quote as urllib_quote
    from urllib.request import urlretrieve

from proteofav.library import aa_codes_1to3_extended

log = logging.getLogger('proteofav.config')
socket.setdefaulttimeout(15)  # avoid infinite hanging

requests_cache.install_cache('proteofav')

__all__ = ["fetch_from_url_or_retry", "row_selector", "constrain_column_types", "exclude_columns",
           "splitting_up_by_key", "merging_down_by_key",
           "flatten_nested_structure", "refactor_key_val_singletons",
           "InputFileHandler", "OutputFileHandler", "Downloader", "GenericInputs"]


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
        log.debug("Fetched data from {}...".format(url))
        return response
    elif response.status_code in retry_in and n_retries >= 0:
        time.sleep(wait)
        return fetch_from_url_or_retry(url, json, header, post, data, retry_in, wait,
                                       (n_retries - 1), stream, **params)
    else:
        try:
            response.raise_for_status()
        except requests.exceptions.HTTPError as e:
            log.warning("%s: Unable to retrieve %s for %s",
                      response.status_code, url, e)


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
        raise ValueError('Column {} does not contain {} value(s)...'
                         ''.format(key, value))

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


def check_sequence(sequence, gap_symbol='-', new_gap_symbol='-', ambiguous='X'):
    """
    Checks an input sequence for uncommon residue symbols.

    :param sequence: (str) protein sequence
    :param gap_symbol: (str) 1-letter symbol for gaps
    :param new_gap_symbol: (str) 1-letter symbol for gaps
    :param ambiguous: (str) 1-letter symbol for ambiguous residues
    :return: returns modified sequence
    """

    new_sequence = "".join([ambiguous if aa not in list(aa_codes_1to3_extended.keys()) else aa
                            for aa in sequence])
    if gap_symbol != new_gap_symbol:
        new_sequence = new_sequence.replace(gap_symbol, new_gap_symbol)

    return new_sequence


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
            log.debug("Downloaded data from {}".format(self._url))
        except (URLError, HTTPError, IOError, Exception) as e:
            log.error("Unable to retrieve %s for %s", self._url, e)

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
            log.debug("Decompressed %s", self._filename)


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
