#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function

import logging
import os
import re
import sys
import random
import time
from datetime import datetime

import requests

from library import valid_ensembl_species_variation

from config import defaults

from Bio import pairwise2

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
    return filter(None, re.split(r'(\d+)', s))


def write_log(message, output_file):
    """
    Appends a message to a log file.

    @param message: message to be printed
    @param output_file: path to the output file
    """
    with open(output_file, "a") as outlog:
        outlog.write(message + "\n")
    return


def get_random_string_with_n_digits(n):
    """
    gets a random number with n digits.

    @param n: number of digits
    @return: returns a random integer (int)
    """

    range_start = 10 ** (n - 1)
    range_end = (10 ** n) - 1
    return random.randint(range_start, range_end)


def get_url_or_retry(url, retry_in=None, wait=1, json=False, header={},
                     **params):
    """
    Fetch an url using Requests or retry fetching it if the server is
    complaining with retry_in error.

    :param url: url to be fetched as a string
    :param wait: sleeping between tries in seconds
    :param params: request.get kwargs.
    :return: url content or url content in json data structure.
    """
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

def isvalid_uniprot_id(identifier):
    """
    'Quickly' checks if a UniProt id is valid.

    :param identifier: testing ID
    :return: boolean
    """
    try:
        if identifier != '':
            url = defaults.http_uniprot + str(identifier)
            get_url_or_retry(url)
            return True
        else:
            raise IDNotValidError
    except IDNotValidError:
        return False
    except requests.HTTPError:
        return False


def isvalid_pdb_id(identifier):
    """
    'Quickly' checks if a PDB id is valid.

    :param identifier: testing ID
    :return: boolean
    """
    try:
        if identifier != '':
            url = defaults.http_pdbe + str(identifier)
            get_url_or_retry(url)
            return True
        else:
            raise IDNotValidError
    except IDNotValidError:
        return False
    except requests.HTTPError:
        return False


def isvalid_ensembl_id(identifier, species='human', variant=False):
    """
    'Quickly' checks if an Ensembl id is valid.

    :param identifier: testing ID
    :param species: Ensembl species
    :param variant: boolean if True uses the variant endpoint
    :return: boolean
    """

    if variant:
        if species not in valid_ensembl_species_variation:
            raise ValueError('Provided species {} is not valid'.format(species))
        ensembl_endpoint = "variation/{}/".format(species)
    else:
        ensembl_endpoint = "lookup/id/"
    try:
        if identifier != '':
            url = defaults.api_ensembl + ensembl_endpoint + str(identifier)
            data = requests.get(url)
            if 'error' not in data:
                return True
            else:
                return False
        else:
            raise IDNotValidError
    except IDNotValidError:
        return False
    except requests.HTTPError:
        return False


def _fetch_sifts_best(identifier, first=False):
    """
    Gets the best structures from the SIFTS endpoint in the
    PDBe api.

    :param identifier: Uniprot ID
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
    :return: boolean
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

    :param from_seq:
    :param to_seq:
    :return:
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


    # Build the global alignment
    alignments = pairwise2.align.globalxx(from_seq, to_seq)
    scores = zip(*alignments)[2]
    message = "Alignment score(s): {}".format(scores)
    logging.info(message)
    message = "First alignment:\n" + pairwise2.format_alignment(*alignments[0])
    logging.debug(message)
    if len(scores) > 1:
        message = "Found multiple alignments, arbitrarily proceeding with the first."
        logging.warning(message)

    # Create the index mapping
    seq_one = aligned_seq_indexes(alignments[0][0])
    seq_two = aligned_seq_indexes(alignments[0][1])
    map = dict(zip(seq_one, seq_two))

    return map


def apply_sequence_index_map(indexes, map):
    """

    :param indexes:
    :return:
    """

    # Perform the raw translation
    translation = []
    for i in indexes:
        equivalent = map.get(i)
        translation.append(equivalent)

    return translation


if __name__ == '__main__':
    # testing routines
    pass
