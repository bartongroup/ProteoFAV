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
    """
    column_types_long = {'CATH_dbAccessionId': 'string',
                         'InterPro_dbAccessionId': 'string',
                         'NCBI_dbAccessionId': 'string',
                         'PDB_dbAccessionId': 'string',
                         'PFAM__dbAccessionId': 'string',
                         'UniProt_dbAccessionId': 'string'
                         }

    column_types_short = {'*_dbAccessionId': 'string',
                          '*_dbChainId': 'string'}


if __name__ == '__main__':
    # testing routines
    pass
