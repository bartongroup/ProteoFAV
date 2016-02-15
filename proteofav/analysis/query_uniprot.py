#!/usr/bin/env python
# -*- coding: utf-8 -*-

import requests
import logging


logger = logging.getLogger(__name__)


# TODO: move this to variants.py?
def query_uniprot(search_terms=('keyword:Disease', 'reviewed:yes', 'organism:human', 'database:(type:pdb)')):
    """
    Query the UniProt API for proteins that have particualar characteristics.

    :param search_terms: A tuple of UniProt Query search terms.
    :return: A list of UniProt IDs
    """
    url = 'http://www.uniprot.org/uniprot'
    params = {'query': ' AND '.join(search_terms),
              'format': 'tab', 'columns': 'id'}
    logger.info('Querying UniProt DB...')
    r = requests.get(url, params=params)
    uniprots = r.content.split('\n')[1:]
    logger.info('Retreived {} UniProt IDs matching query.'.format(len(uniprots)))

    return uniprots


if __name__ == '__main__':
    uniprot_ids = query_uniprot()
    for i in uniprot_ids:
        print i