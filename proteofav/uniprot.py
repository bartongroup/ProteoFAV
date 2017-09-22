#!/usr/bin/env python
# -*- coding: utf-8
"""
Created on 17:26 19/02/2016 2016 
Define auxiliary functions for interacting with Uniprot.
"""
from __future__ import absolute_import
import logging

try:
    # python 2.7
    from StringIO import StringIO
except ImportError:
    from io import StringIO, BytesIO
try:
    # python 2.7
    from urlparse import parse_qs
except ImportError:
    from urllib.parse import parse_qs

import pandas as pd

from proteofav.config import defaults
from proteofav.utils import get_url_or_retry
from proteofav.utils import row_selector

__all__ = ["fetch_uniprot_sequence",
           "fetch_uniprot_formal_specie",
           "_uniprot_info",
           "_fetch_uniprot_gff",
           "map_gff_features_to_sequence",
           "_uniprot_to_ensembl_xref"]
log = logging.getLogger('proteofav.config')


def fetch_uniprot_sequence(uniprot_id):
    """
    Gets current sequence of a Uniprot entry.

    :param str uniprot_id: Uniprot accession
    :return str: the sequence

    >>> print(fetch_uniprot_formal_specie('P17612'))
    Homo sapiens

    """

    return _uniprot_info(uniprot_id, cols='sequence').iloc[0, 1]


def fetch_uniprot_formal_specie(uniprot_id, remove_isoform=True):
    """
    Gets the species name of an organism expressing a protein.

    :param Bool remove_isoform: whether to remove the isoform identifier.
    :param str uniprot_id: Uniprot accession
    :return: the species name (two words)
    :rtype: str or None

    >>> print(fetch_uniprot_sequence('P17612'))[:20]
    MGNAAAAKKGSEQESVKEFL

    """
    if remove_isoform:
        uniprot_id = uniprot_id.split('-')[0]

    full_specie = _uniprot_info(uniprot_id, cols='organism').iloc[0, 1]

    try:
        return " ".join(full_specie.split()[0:2])
    except AttributeError:
        log.error('Could not retrieve {} information. Maybe obsolete UniProt accession?'.format(
            uniprot_id))
        return None


def _uniprot_info(uniprot_id, retry_in=(503, 500), cols=None):
    """
    Retrieve Uniprot information from the database.

    :param str uniprot_id: Uniprot accession identifier
    :param retry_in: iterable of status_code to be retried upon error.
    :type retry_in: list of [int]
    :return pandas.DataFrame: table from Uniprot.
    Default table columns:
    :raises requests.HTTPError: when hits a bad status_code
    """

    if not cols:
        cols = ('id', 'entry name', 'reviewed', 'protein names', 'genes', 'organism',
                'sequence', 'length')
    elif isinstance(cols, str):
        cols = ('id', cols)

    params = {'query': 'accession:' + str(uniprot_id),
              'columns': ",".join(cols),
              'format': 'tab',
              'contact': ""}
    url = defaults.api_uniprot
    response = get_url_or_retry(url=url, retry_in=retry_in, **params)
    try:
        data = pd.read_table(StringIO(response))
    except TypeError:
        # python 3.5
        data = pd.read_table(BytesIO(response))
    except ValueError as e:
        log.error(e)
        return None
    # id column is called Entry in the table
    return row_selector(data, 'Entry', uniprot_id)


def _fetch_uniprot_gff(uniprot_id):
    """
    Retrieve Uniprot data from the GFF file.

    :param str uniprot_id: Uniprot accession
    :return pandas.DataFrame: table
    :raises requests.HTTPError: when hits a bad status_code
    """
    url = defaults.api_uniprot + uniprot_id + ".gff"
    cols = "NAME SOURCE TYPE START END SCORE STRAND FRAME GROUP empty".split()

    data = pd.read_table(url, skiprows=2, names=cols)
    groups = data.GROUP.apply(parse_qs)
    groups = pd.DataFrame.from_records(groups)

    return data.merge(groups, left_index=True, right_index=True)


def map_gff_features_to_sequence(uniprot_id,
                                 query_type='',
                                 group_residues=True,
                                 drop_types=('Helix', 'Beta strand', 'Turn', 'Chain')):
    """
    Map Uniprot GFF features to the protein sequence.

    :param str uniprot_id: Uniprot accession
    :param str query_type: Select type of feature
    :param bool group_residues: by default each row in the resulting table,
        maps to a residue. When set to False, each row represent a feature
        per residue.

    :param tuple drop_types: Filter out some of the features, important to
        remove fetures that spam.

    :return pd.DataFrame: table. Columns will depend on paramenters.
    """

    def annotation_writer(gff_row):
        """
        Establish a set of rules to annotate Uniprot GFF.

        :param pd.Series gff_row: each line in the GFF file.
        :return str: template filled with type-specific fields.
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

    if table.empty:
        return table

    if group_residues:

        return table.groupby('idx').agg({'annotation': lambda x: ', '.join(x)})
    else:

        return table


def _uniprot_to_ensembl_xref(uniprot_id, species='homo_sapiens'):
    """
    Return Gene, transcripts and translational ids from Ensembl to Uniprot.
    Ensembl -> Uniprot reference is ideal because Ensembl database change more
    often the Uniprot'smove quicker than Uniprot.

    :param str uniprot_id: Uniprot accession
    :param str species: species name
    :return pandas.DataFrame: table with columns
    """

    url = "{}xrefs/symbol/{}/{}?content-type=application/json".format(
        defaults.api_ensembl, species, uniprot_id)

    return pd.read_json(url)


if __name__ == '__main__':
    pass
