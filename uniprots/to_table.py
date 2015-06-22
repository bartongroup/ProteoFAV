#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Created on 11/06/2015

"""
from StringIO import StringIO

import logging

import pandas as pd

from utils.config import defaults
from utils.utils import request_info_url
from utils.utils import get_url_or_retry

logger = logging.getLogger(__name__)

# def _uniprot_info_to_table(identifier, verbose=False):
#     """
#     Fetches some information including the sequence of a particular
#     UniProt entry.
#
#     :param identifier: UniProt accession identifier
#     :param verbose: boolean
#     :return: pandas table dataframe
#     """
#
#     information = {}
#     rows = []
#
#     params = {'query': 'accession:' + identifier,
#               'columns': 'entry name,reviewed,protein names,genes,organism,sequence,length',
#               'format': 'tab',
#               'contact': defaults.contact_email}
#     request = request_info_url(defaults.http_uniprot, params, verbose=verbose)
#
#     data = request.text.split('\n')
#     for i, line in enumerate(data):
#         if i == 1 and line != '':
#             line = line.split('\t')
#             information['Name'] = line[0]
#             information['Status'] = line[1]
#             information['Protein'] = line[2]
#             information['Genes'] = line[3]
#             information['Organism'] = line[4]
#             information['Sequence'] = line[5]
#             information['Length'] = int(line[6])
#
#     rows.append(information)
#     return pd.DataFrame(rows)

def _uniprot_info_to_table(identifier, retry_in=(503, 500), cols=None):
    """
    Retrive uniprot information from the database.

    :param identifier: UniProt accession identifier
    :param verbose: boolean
    :return: pandas table dataframe
    """

    if not cols:
        cols = ('entry name', 'reviewed', 'protein names', 'genes', 'organism',
                'sequence', 'length')
    elif isinstance(cols, str):
        cols = ('entry name', cols)

    params = {'query': 'accession:' + identifier,
              'columns': ",".join(cols),
              'format': 'tab',
              'contact': ""}
    url = "http://www.uniprot.org/uniprot/"
    response = get_url_or_retry(url=url, retry_in=retry_in, **params)
    try:
        data = pd.read_table(StringIO(response))
    except ValueError as e:
        #logger
        data = response
    return data


def _uniprot_ensembl_mapping_to_table(identifier, verbose=False):
    """
    Uses the UniProt mapping service to try and get Ensembl IDs for
    the UniProt accession identifier provided.

    :param identifier: UniProt accession identifier
    :param verbose: boolean
    :return: pandas table dataframe
    """

    information = {}
    rows = []

    # TODO: keeps failing due to server issues - perhaps use Ensembl endpoints
    # for this mapping
    ensembl_mappings = ["ENSEMBL_ID", "ENSEMBL_PRO_ID", "ENSEMBL_TRS_ID"]
    for ensembl in ensembl_mappings:
        params = {'from': 'ACC',
                  'to': ensembl,
                  'format': 'tab',
                  'query': identifier,
                  'contact': defaults.contact_email}

        request = request_info_url(defaults.http_uniprot_mapping, params,
                                   verbose=verbose)

        data = request.text.split('\n')
        for i, line in enumerate(data):
            if i >= 1 and line != '':
                line = line.split('\t')
                try:
                    if line[1] in information[ensembl]:
                        continue
                    information[ensembl].append(line[1])
                except KeyError:
                    information[ensembl] = line[1]
                except AttributeError:
                    information[ensembl] = [information[ensembl]]
                    information[ensembl].append(line[1])

    rows.append(information)
    return pd.DataFrame(rows)



if __name__ == '__main__':
    # testing routines
    pass
