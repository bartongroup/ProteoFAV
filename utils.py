#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function

import logging
import os
import re
import sys
import time
import urllib
import gzip
import shutil
import socket
import colorsys
import requests
from os import path
from datetime import datetime
import numpy as np
from Bio import pairwise2
from library import valid_ensembl_species_variation
from config import defaults

socket.setdefaulttimeout(15)
log = logging.getLogger(__name__)


def is_valid_file(parser, arg):
    """
    Check if arg is a valid file and throw a parsing error if not.

    :param parser: argparse.ArgumentParser
    :param arg: argument
    :return: Open file handle
    """
    try:
        return open(arg, 'r')
    except:
        parser.error("Not a valid file: %s" % arg)


def delete_file(filename):
    """

    :param filename: File to delete
    :return: None
    """
    os.remove(filename)


class IDNotValidError(Exception):
    """
    Base class for database related exceptions.
    Databases: UniProt, PDB, Ensembl, etc.
    """
    pass


def current_date():
    """
    Gets the current date.

    :return: outputs formatted date as 'Day/Month/Year'
    :rtype: str
    """
    date = datetime.now()
    month = date.strftime("%b")
    day = date.strftime("%d")
    year = date.strftime("%Y")
    return "{}/{}/{}".format(day, month, year)


def current_time():
    """
    Gets current date and time.

    :return: outputs formatted time as 'Day/Month/Year H:M:S'
    :rtype: str
    """
    date = datetime.now()
    year = date.strftime("%Y")
    month = date.strftime("%m")
    day = date.strftime("%d")
    hour = date.strftime("%H")
    minute = date.strftime("%M")
    second = date.strftime("%S")
    return "{}/{}/{} {}:{}:{}".format(day, month, year, hour, minute, second)


def create_directory(directory):
    """
    Creates a directory structure if it does not exist.

    :param directory: directory name (expects full path)
    :return: creates a directory if it does not exist yet
    """
    return os.makedirs(directory)


def flash(message):
    """
    Flashes a message out.

    :param message: input message str()
    """
    print(str(message))
    sys.stdout.flush()
    return


def string_split(s):
    """
    Splits a string by its numeric values:

    :param s: input string
    :return: a list of strings
    :rtype: list

    :Example:

        >>> a = "foo234bar"
        >>> b = string_split(a)
        >>> print(b)
        ['foo', '234', 'bar']

    """
    return filter(None, re.split(r'(\d+)', s))


def get_url_or_retry(url, retry_in=None, wait=1, json=False, header=None,
                     **params):
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

    if response.ok:
        if json:
            return response.json()
        else:
            return response.content
    elif response.status_code in retry_in:
        time.sleep(wait)
        return get_url_or_retry(
                url, retry_in, wait, json, header, **params)
    else:
        print(response.url)
        log.error(response.status_code)
        response.raise_for_status()


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
        url = getattr(defaults, 'http_' + database) + identifier
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
            url = defaults.api_ensembl + ensembl_endpoint + urllib.quote(identifier, safe='')
            data = requests.get(url)
            if data.status_code is not 200:
                return False
            elif 'error' not in data.text:
                return True
            else:
                return False
        else:
            raise IDNotValidError
    except IDNotValidError:
        return False
    except requests.HTTPError:
        return False


def fetch_sifts_best(identifier, first=False):
    """
    Gets the best structures from the SIFTS endpoint in the
    PDBe api.

    :param identifier: UniProt ID
    :param first: gets the first entry
    :return: url content or url content in json data structure.
    """

    sifts_endpoint = "mappings/best_structures/"
    url = defaults.api_pdbe + sifts_endpoint + str(identifier)
    response = get_url_or_retry(url, json=True)
    return response if not first else response[identifier][0]


def compare_uniprot_ensembl_sequence(sequence1, sequence2,
                                     permissive=True):
    """
    Compares two given sequences in terms of length and
    sequence content.

    :param sequence1: sequence 1
    :param sequence2: sequence 2
    :param permissive: if True it allows for same size
                      sequences with mismatches
    :return: simply a true or false
    :rtype: boolean
    """

    if permissive:
        if len(sequence1) == len(sequence2):
            return True
        else:
            return False
    else:
        if not len(sequence1) == len(sequence2):
            return False
        else:
            for i, j in zip(sequence1, sequence2):
                if i != j:
                    return False
            return True


def count_mismatches(sequence1, sequence2):
    """
    Counts the number of mismatches between two sequences
    of the same length.

    :param sequence1: sequence 1
    :param sequence2: sequence 2
    :return: The number of mismatches between sequences 1 and 2.
    """
    if len(sequence1) == len(sequence2):
        mismatches = []
        for i, j in zip(sequence1, sequence2):
                    if i != j:
                        mismatches.append((i, j))
    else:
        raise ValueError('Sequences are different lengths.')
    return(len(mismatches))


def map_sequence_indexes(from_seq, to_seq):
    """
    Gets a map between sequences.

    :param from_seq: input sequence
    :param to_seq: input sequence
    :return: a map between sequences.
    :rtype: dict
    """

    def aligned_seq_indexes(seq):
        seq_indexes = []
        i = 0
        for res in seq:
            if res != '-':
                seq_indexes.append(i)
                i += 1
            else:
                seq_indexes.append('-')
        return seq_indexes

    # build the local alignment
    alignments = pairwise2.align.localxx(from_seq, to_seq)
    scores = zip(*alignments)[2]
    message = "Alignment score(s): {}".format(scores)
    logging.info(message)
    message = "First alignment:\n" + pairwise2.format_alignment(*alignments[0])
    logging.debug(message)
    if len(scores) > 1:
        message = "Found multiple alignments, arbitrarily proceeding with the first."
        logging.warning(message)

    # create the index mapping
    seq_one = aligned_seq_indexes(alignments[0][0])
    seq_two = aligned_seq_indexes(alignments[0][1])
    outmap = dict(zip(seq_one, seq_two))

    return outmap


# TODO: documentation
def apply_sequence_index_map(indexes, imap):
    """

    :param indexes:
    :param map:
    :return:
    """

    # perform the raw translation
    translation = []
    for i in indexes:
        equivalent = imap.get(i)
        translation.append(equivalent)

    return translation


def get_colors(num_colors):
    """
    See
    http://stackoverflow.com/questions/470690/how-to-automatically-generate-n-distinct-colors

    :param num_colors: number of color
    :return: a list of colors
    :rtype: list
    """
    colors = []
    for i in np.arange(0., 360., 360. / num_colors):
        hue = i / 360.
        lightness = (50 + np.random.rand() * 10) / 100.
        saturation = (90 + np.random.rand() * 10) / 100.
        colors.append(colorsys.hls_to_rgb(hue, lightness, saturation))
    return colors


def _mmcif_unit_cell(pdb_id):
    """
    Loader of mmCIF unit cell parameters.

    :param pdb_id: PDB id
    :return: pandas table dataframe
    :rtype: dict
    """

    cif_path = path.join(defaults.db_mmcif, pdb_id + '.cif')

    lines = []
    with open(cif_path) as inlines:
        for line in inlines:
            if line.startswith("_cell."):
                lines.append(line.split('.')[1].rstrip())

    l = [i.split() for i in lines]
    d = {k: v for k, v in l}

    return d


def fractional_to_cartesian(coords, pdb_id, matrix_only=False):
    """
    Converts fractional unit cell coords to cartesian.

    :param coords: Atom coordinates
    :param pdb_id: PDB id
    :param matrix_only: boolean
    :return: cartesian coordinates
    :rtype: np.array
    """

    # retrieve and parse unit cell parameters
    d = _mmcif_unit_cell(pdb_id)
    a2r = np.pi / 180.
    alpha = a2r * float(d['angle_alpha'])
    beta = a2r * float(d['angle_beta'])
    gamma = a2r * float(d['angle_gamma'])
    a = float(d['length_a'])
    b = float(d['length_b'])
    c = float(d['length_c'])

    # unit cell volume
    v = np.sqrt(1 - np.cos(alpha) * np.cos(alpha) - np.cos(beta) * np.cos(beta) -
                np.cos(gamma) * np.cos(gamma) + 2 * np.cos(alpha) * np.cos(beta) * np.cos(gamma))

    # build the transformation matrix
    tr = np.matrix([
        [a, b * np.cos(gamma), c * np.cos(beta)],
        [0, b * np.sin(gamma),
         c * ((np.cos(alpha) - np.cos(beta) * np.cos(gamma)) / np.sin(gamma))],
        [0, 0, c * (v / np.sin(gamma))]
    ])

    if matrix_only:
        return tr

    # now apply the transformation to the coordinates
    # TODO: type check this?
    coords = np.matrix(coords)
    coords = coords * tr

    # return the Nx3 results
    return np.array(coords)


def autoscale_axes(xyz, margin=5):
    """
    Auto scales the xyz axes by a margin.

    :param xyz: Atom coordinates
    :param margin: spacial margin size/length
    :return: array
    """

    x = xyz[:, 0]
    y = xyz[:, 1]
    z = xyz[:, 2]

    rx = max(x) - min(x)
    ry = max(y) - min(y)
    rz = max(z) - min(z)

    max_range = max([rx, ry, rz]) + margin

    pad = []
    for i in [rx, ry, rz]:
        pad.append((max_range - i) / 2)

    return [[max(x) + pad[0], min(x) - pad[0]],
            [max(y) + pad[1], min(y) - pad[1]],
            [max(z) + pad[2], min(z) - pad[2]]]


def fetch_files(identifier, directory=None, sources=("cif", "dssp", "sifts")):
    """
    Small routine to fetch data files from their respective repositories.

    Use defaults to fetch files from servers. Defaults server are defined in
     the default config.txt. Returns None since it has side effects.

    Created on 11:41 24/07/15 2015
    :param identifier: protein identifier as PDB identifier
    :param directory: path to download. The default downloads all files to
     default.temp folder. If is a string and valid path downloads all files to
     that folder. If its a iterable and all valid path downloads to the
     respective folders
    :param sources: where to fetch the data. Must be in the config file.
    :return: None
    """
    if isinstance(sources, str):
        sources = [sources]
    if not directory:
        directory = defaults.temp
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
            urllib.urlretrieve(url, destination + filename)
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
    # testing routines
    pass
