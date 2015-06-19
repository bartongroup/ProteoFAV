#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Created on 11/06/2015

"""

import sys
sys.path.insert(0, '../')
import logging
import json
import pandas as pd

from utils.config import defaults
from utils.utils import request_info_url
from utils.utils import isvalid_ensembl

logger = logging.getLogger(__name__)


def _transcript_variants_ensembl_to_table(identifier, verbose=False):
    """
    Queries the Ensembl API for transcript variants (mostly dbSNP)
    based on Ensembl Protein identifiers (e.g. ENSP00000326864).

    :param identifier: Ensembl Protein ID
    :param verbose: boolean
    :return: pandas table dataframe
    """

    if not isvalid_ensembl(identifier):
        raise ValueError("{} is not a valid Ensembl Accession.".format(identifier))

    ensembl_endpoint = 'overlap/translation/'
    params = {'feature': 'transcript_variation',
              'content-type': 'application/json'}
    request = request_info_url("{}{}{}".format(defaults.api_ensembl,
                                               ensembl_endpoint,
                                               str(identifier)),
                               params=params,
                               verbose=verbose)
    rows = json.loads(request.text)
    return pd.DataFrame(rows)


def _somatic_variants_ensembl_to_table(identifier, verbose=False):
    """
    Queries the Ensembl API for somatic transcript variants (COSMIC)
    based on Ensembl Protein identifiers (e.g. ENSP00000326864).

    :param identifier: Ensembl Protein ID
    :param verbose: boolean
    :return: pandas table dataframe
    """

    if not isvalid_ensembl(identifier):
        raise ValueError("{} is not a valid Ensembl Accession.".format(identifier))

    ensembl_endpoint = 'overlap/translation/'
    params = {'feature': 'somatic_transcript_variation',
              'content-type': 'application/json'}
    request = request_info_url("{}{}{}".format(defaults.api_ensembl,
                                               ensembl_endpoint,
                                               str(identifier)),
                               params=params,
                               verbose=verbose)
    rows = json.loads(request.text)
    return pd.DataFrame(rows)


def _ensembl_variant_to_table(identifier, verbose=False):
    """
    Queries the Ensembl API for variant IDs (e.g rs376845802 or COSM302853).

    :param identifier: variant ID
    :param verbose: boolean
    :return: pandas table dataframe
    """

    if not isvalid_ensembl(identifier, variant=True):
        raise ValueError("{} is not a valid Variation Accession.".format(identifier))


    # TODO: fix this assuming human variation
    ensembl_endpoint = 'variation/human/'
    params = {'content-type': 'application/json'}
    # other params are {'pops': '1', 'phenotypes': '1', 'genotypes': '1'}
    request = request_info_url("{}{}{}".format(defaults.api_ensembl,
                                               ensembl_endpoint,
                                               str(identifier)),
                               params=params,
                               verbose=verbose)
    data = json.loads(request.text)

    rows = []
    information = {}
    for parent in data:
        if parent == "mappings":
            for entry in data[parent]:
                for key in entry:
                    try:
                        if entry[key] in information[key]:
                            continue
                        information[key].append(entry[key])
                    except KeyError:
                        information[key] = entry[key]
                    except AttributeError:
                        information[key] = [information[key]]
                        information[key].append(entry[key])
        else:
            information[parent] = data[parent]

    rows.append(information)
    return pd.DataFrame(rows)

if __name__ == '__main__':
    # testing routines
    pass
