#!/usr/bin/env python
# -*- coding: utf-8
import logging
import requests
import time

# import pandas as pd
#
from pdbs.to_table import _mmcif_atom_to_table as mmcif_to_table
from pdbs.to_table import _dssp_to_table as dssp_to_table
from sifts.to_table import _sifts_residues_to_table as sifts_to_table

# import utils

logger = logging.getLogger(__name__)

def _fetch_sifts_best(uniprot_id, first=False):
    url = "http://wwwdev.ebi.ac.uk/pdbe/api/mappings/best_structures/"
    url = url + uniprot_id
    response = get_url_or_retry(url, json=True)
    return response if not first else response[uniprot_id][0]

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
        logger.error(response.status_code)
        logger.error(response.url)
        response.raise_for_status()

def merge_tables(uniprot_id, pdb_id=None, chain=None):
    """
    Merges the output from multiple to table method.
    if no pdb_id uses sifts_best_structure
    if no chain uses first
    or if use_all_chains use all chains of the pdb_id
    """
    if not pdb_id:
        best_pdb = _fetch_sifts_best(uniprot_id, first=True)
        pdb_id = best_pdb['pdb_id']
        chain = best_pdb['chain_id']
        # start stop
    if not chain:
        raise NotImplementedError
        best_pdb = _fetch_sifts_best(uniprot_id)[pdb_id]
        chain = best_pdb['chain_id']







class Structural(object):
    pass

if __name__ == '__main__':
    merge_tables("Q16566")