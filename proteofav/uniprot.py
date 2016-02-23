#!/usr/bin/env python
# -*- coding: utf-8
"""
Created on 17:26 19/02/2016 2016 

"""
from StringIO import StringIO
from urlparse import parse_qs

import pandas as pd

from proteofav import defaults
from proteofav.utils import get_url_or_retry
from proteofav.variants import log


def get_uniprot_sequence(uniprot_id):
    """

    :param uniprot_id:
    :return:
    """
    return _uniprot_info(uniprot_id, cols='sequence').iloc[0, 1]


def get_uniprot_formal_specie(uniprot_id):
    """

    :param uniprot_id:
    """
    full_specie = _uniprot_info(uniprot_id, cols='organism').iloc[0, 1]
    return full_specie.rsplit(None, 1)[0]


def _uniprot_info(uniprot_id, retry_in=(503, 500), cols=None, check_id=False):
    """
    Retrieve uniprot information from the database.

    :param uniprot_id: UniProt accession identifier
    :return: pandas table dataframe
    """

    if not cols:
        cols = ('entry name', 'reviewed', 'protein names', 'genes', 'organism',
                'sequence', 'length')
    elif isinstance(cols, str):
        cols = ('entry name', cols)

    params = {'query': 'accession:' + str(uniprot_id),
              'columns': ",".join(cols),
              'format': 'tab',
              'contact': ""}
    url = defaults.api_uniprot
    response = get_url_or_retry(url=url, retry_in=retry_in, **params)
    try:
        data = pd.read_table(StringIO(response))
    except ValueError as e:
        log.error(e)
        data = response
    return data


def _fetch_uniprot_gff(uniprot_id):
    """
    Retrieve UniProt data from the GFF file

    :param uniprot_id: UniProt accession identifier
    :type uniprot_id: str
    :return: table
    :return type: pandas.DataFrame
    """
    url = defaults.api_uniprot + uniprot_id + ".gff"
    cols = "NAME SOURCE TYPE START END SCORE STRAND FRAME GROUP empty".split()

    data = pd.read_table(url, skiprows=2, names=cols)
    groups = data.GROUP.apply(parse_qs)
    groups = pd.DataFrame.from_records(groups)
    data = data.merge(groups, left_index=True, right_index=True)

    return data


def map_gff_features_to_sequence(uniprot_id, query_type='', group_residues=True, drop_types=(
        'Helix', 'Beta strand', 'Turn', 'Chain')):
    """Remaps features in the uniprot gff file to the sequence.

    :param group_residues:
    :param query_type:
    :param uniprot_id: UniProt-SP accession
    :param drop_types: Annotation type to be dropped
    :return: table read to be joined to main table
    """

    def annotation_writer(gff_row):
        """
        Establish a set of rules to annotate uniprot GFF.

        :param gff_row: each line in the GFF file
        :return: template filled with type-specific fields.
        """
        if not gff_row.ID and not gff_row.Note:
            return gff_row.TYPE
        elif not gff_row.ID:
            return '{0.TYPE}: {0.Note}'.format(gff_row)
        elif not gff_row.Note:
            return '{0.TYPE} ({0.ID})'.format(gff_row)
        else:
            return '{0.TYPE}: {0.Note} ({0.ID})'.format(gff_row)

    table = _fetch_uniprot_gff(uniprot_id)
    if query_type:
        table = table[table.TYPE == query_type]
    elif drop_types:
        table = table[~table.TYPE.isin(drop_types)]

    lines = []
    for i, row in table.iterrows():
        lines.extend({'idx': i, 'annotation': annotation_writer(row)}
                     for i in range(row.START, row.END + 1))
    table = pd.DataFrame(lines)
    if group_residues:
        return table.groupby('idx').agg({'annotation': lambda x: ', '.join(x)})
    else:
        return table


def _uniprot_to_ensembl_xref(uniprot_id, species='homo_sapiens'):
    """
        Return Gene, transcripts and translational ids from Ensembl to Uniprot.
        Ensembl -> Uniprot reference if better than otherwise due Ensembl
        move quicker than Uniprot.

        :param uniprot_id: Uniprot ID
        :param species: organism name
        :return:
        :rtype: pandas.DataFrame
        """
    url = "{}xrefs/symbol/{}/{}?content-type=application/json".format(
            defaults.api_ensembl, species, uniprot_id)
    return pd.read_json(url)


if __name__ == '__main__':
    # testing routines
    pass
