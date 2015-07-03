#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import print_function
import logging
import os
from os import path
import re
import sys
import random
import time
import json
from datetime import datetime
import numpy as np

import requests

from config import defaults, Defaults
import to_table

log = logging.getLogger(__name__)

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


def request_info_url(url, params=None, verbose=False):
    """
    Gets content from the provided url.

    If it fails, tries to access it a few more times to make sure it
    was not a temporary server hang.

    @param url: input URL
    @param params: parameters to be encoded and requested
    @param verbose: boolean
    @return: returns a data object from *requests*
    """

    if params is None:
        params = {}

    # limiting ourselves to up to 10 requests per second
    time.sleep(0.1)

    try:
        request = requests.get(url, params=params)
    except Exception as e:
        message = "{}\t{}".format(url, e)
        flash(message)
        raise Exception(message)

    if verbose:
        message = "{}\t{}".format(url, request.status_code)
        flash(message)

    if request.ok:
        return request
    # Ensembl status for exceeding number of connections per second
    elif request.status_code == 429:
        time.sleep(0.5)
        request_info_url(url, params=params, verbose=verbose)
    else:
        message = "{}\t{}".format(url, request.status_code)
        flash(message)
        raise Exception(message)


def isvalid_uniprot(identifier):
    """
    'Quickly' checks if a UniProt id is valid.

    :param identifier: testing ID
    :return: boolean
    """

    # TODO: Improve this!
    try:
        if identifier != '':
            request_info_url(defaults.http_uniprot + identifier)
            return True
        else:
            raise Exception
    except Exception:
        return False


def isvalid_pdb(identifier):
    """
    'Quickly' checks if a PDB id is valid.

    :param identifier: testing ID
    :return: boolean
    """

    # TODO: Improve this!
    try:
        if identifier != '':
            request_info_url("{}{}".format(defaults.http_pdbe,
                                           str(identifier)))
            return True
        else:
            raise Exception
    except Exception:
        return False


def isvalid_ensembl(identifier, variant=False):
    """
    'Quickly' checks if an Ensembl id is valid.

    :param identifier: testing ID
    :param variant: boolean if True uses the variant endpoint
    :return: boolean
    """

    # TODO: Improve this!

    if variant:
        ensembl_endpoint = 'variation/human/'
    else:
        ensembl_endpoint = 'lookup/id/'
    params = {'content-type': 'application/json'}
    try:
        if identifier != '':
            request = request_info_url("{}{}{}".format(defaults.api_ensembl,
                                                       ensembl_endpoint,
                                                       str(identifier)),
                                       params=params)
            print(request.url)
            data = json.loads(request.text)
            if 'error' not in data:
                return True
            else:
                return False
        else:
            raise Exception
    except Exception:
        return False


def compare_uniprot_ensembl_sequence():
    pass


def get_url_or_retry(url, retry_in=(), wait=1, json=False, header={}, **params):
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
        log.error(response.status_code)
        response.raise_for_status()


def _fetch_sifts_best(uniprot_id, first=False):
    url = 'http://wwwdev.ebi.ac.uk/pdbe/api/mappings/best_structures/'
    url = url + uniprot_id
    response = get_url_or_retry(url, json=True)
    return response if not first else response[uniprot_id][0]


def merge_tables(uniprot_id, pdb_id=None, chain=None, groupby='CA',
                 defaults=None):
    """
    Merges the output from multiple to table method.
    if no pdb_id uses sifts_best_structure
    if no chain uses first
    or if use_all_chains use all chains of the pdb_id
    """

    def to_list(series):
        return series.tolist()

    def to_unique(series):
        return series.unique()

    if not defaults:
        from config import defaults

    groupby_opts = {'list': {'label_comp_id': 'unique',
                             'label_atom_id': to_list,
                             'label_asym_id': 'unique',
                             'Cartn_x': to_list,
                             'Cartn_y': to_list,
                             'Cartn_z': to_list,
                             'occupancy': 'unique',
                             'B_iso_or_equiv': np.mean},
                    'CA': {'label_comp_id': to_unique,
                           'label_atom_id': to_unique,
                           'label_asym_id': to_unique,
                           'Cartn_x': to_list,
                           'Cartn_y': to_list,
                           'Cartn_z': to_list,
                           'occupancy': to_unique,
                           'B_iso_or_equiv': to_unique}}  # TODO centroid?

    try:
        groupby_opts[groupby]
    except KeyError:
        raise TypeError('Groupby paramenter should be one of {}'.format(
            ', '.join(groupby_opts.keys())))

    if not pdb_id:
        best_pdb = _fetch_sifts_best(uniprot_id, first=True)
        pdb_id = best_pdb['pdb_id']
        chain = best_pdb['chain_id']
    if not chain:
        best_pdb = _fetch_sifts_best(uniprot_id)[pdb_id]
        chain = best_pdb['chain_id']

    cif_path = path.join(defaults.db_mmcif, pdb_id + ".cif")
    dssp_path = path.join(defaults.db_dssp, pdb_id + ".dssp")
    sifts_path = path.join(defaults.db_sifts, pdb_id + ".xml")

    cif_table = to_table._mmcif_atom_to_table(cif_path)
    cif_table = cif_table.query("label_asym_id == @chain & group_PDB == 'ATOM'")
    # next line raises a SettingWithCopyWarning but seems to be correct
    cif_table.loc[:, ("label_seq_id")] = cif_table.loc[:, "label_seq_id"].astype(np.int)
    # TODO keep heteroatoms?

    dssp_table = to_table._dssp_to_table(dssp_path)
    dssp_table = dssp_table.query('chain_id == @chain')
    dssp_table.set_index(['icode'], inplace=True)

    sifts_table = to_table._sifts_residues_to_table(sifts_path)
    sifts_table = sifts_table.query('PDB_dbChainId == @chain')
    sifts_table.set_index(['PDB_dbResNum'], inplace=True)

    if groupby == 'CA':
        cif_table = cif_table.query("label_atom_id == 'CA'")
    cif_table = cif_table.groupby('label_seq_id').agg(groupby_opts[groupby])

    # This will remove all atoms annotated in sifts but not in mmcif seqres
    # To keep those, sift_table should be the left in join
    cif_sifts = cif_table.join(sifts_table)
    cif_sifts_dssp = cif_sifts.join(dssp_table)
    return cif_sifts_dssp

if __name__ == '__main__':
    # testing routines

    pass