#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Created on 11/06/2015

"""

import sys
sys.path.insert(0, '../')
import logging

import pandas as pd

from utils.config import get_config
from utils.utils import request_info_url
from utils.utils import isvalid_uniprot

logger = logging.getLogger(__name__)


def _uniprot_info_to_table(identifier, verbose=False):
    """
    Fetches some information including the sequence of a particular
    UniProt entry.

    :param identifier: UniProt accession identifier
    :param verbose: boolean
    :return: pandas table dataframe
    """

    if not isvalid_uniprot(identifier):
        raise ValueError("{} is not a valid UniProt Accession.".format(identifier))

    information = {}
    rows = []

    params = {'query': 'accession:%s' % identifier,
              'columns': 'entry name,reviewed,protein names,genes,organism,sequence,length',
              'format': 'tab',
              'contact': 'fmmarquesmadeira@dundee.ac.uk'}
    config = get_config('http_uniprot')
    request = request_info_url(config.http_uniprot, params, verbose=verbose)

    data = request.text.split('\n')
    for i, line in enumerate(data):
        if i == 1 and line != '':
            line = line.split('\t')
            information['Name'] = line[0]
            information['Status'] = line[1]
            information['Protein'] = line[2]
            information['Genes'] = line[3]
            information['Organism'] = line[4]
            information['Sequence'] = line[5]
            information['Length'] = int(line[6])

    rows.append(information)
    return pd.DataFrame(rows)


def _uniprot_ensembl_mapping_to_table(identifier, verbose=False):
    """
    Uses the UniProt mapping service to try and get Ensembl IDs for
    the UniProt accession identifier provided.

    :param identifier: UniProt accession identifier
    :param verbose: boolean
    :return: pandas table dataframe
    """

    if not isvalid_uniprot(identifier):
        raise ValueError("{} is not a valid UniProt Accession.".format(identifier))

    information = {}
    rows = []

    ensembl_mappings = ["ENSEMBL_ID", "ENSEMBL_PRO_ID", "ENSEMBL_TRS_ID"]
    for ensembl in ensembl_mappings:
        params = {'from': 'ACC',
                  'to': ensembl,
                  'format': 'tab',
                  'query': identifier,
                  'contact': 'fmmarquesmadeira@dundee.ac.uk'}

        config = get_config('http_uniprot_mapping')
        request = request_info_url(config.http_uniprot_mapping, params,
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
