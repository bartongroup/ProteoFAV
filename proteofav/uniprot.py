#!/usr/bin/env python
# -*- coding: utf-8
"""
Created on 17:26 19/02/2016 2016 

"""
from urlparse import parse_qs

import pandas as pd

from proteofav.config import defaults

__author__ = 'tbrittoborges'


def map_gff_features_to_sequence(identifier, query_type='', drop_types=(
        'Helix', 'Beta strand', 'Turn', 'Chain')):
    """Remaps features in the uniprot gff file to the sequence.

    :param query_type:
    :param identifier: UniProt-SP accession
    :param drop_types: Annotation type to be dropped
    :return: table read to be joined to main table
    """

    def annotation_writer(gff_row):
        """
        Establish a set of rules to annotate uniprot GFF.

        :param gff_row: each line in the GFF file
        :return:
        """
        if not gff_row.ID and not gff_row.Note:
            return gff_row.TYPE
        elif not gff_row.ID:
            return '{0.TYPE}: {0.Note}'.format(gff_row)
        elif not gff_row.Note:
            return '{0.TYPE} ({0.ID})'.format(gff_row)
        else:
            return '{0.TYPE}: {0.Note} ({0.ID})'.format(gff_row)

    table = _fetch_uniprot_gff(identifier)
    if query_type:
        table = table[table.TYPE == query_type]
    elif drop_types:
        table = table[~table.TYPE.isin(drop_types)]

    lines = []
    for i, row in table.iterrows():
        lines.extend({'idx': i, 'annotation': annotation_writer(row)}
                     for i in range(row.START, row.END + 1))
    table = pd.DataFrame(lines)
    return table.groupby('idx').agg({'annotation': lambda x: ', '.join(x)})


def _fetch_uniprot_gff(identifier):
    """
    Retrieve UniProt data from the GFF file

    :param identifier: UniProt accession identifier
    :type identifier: str
    :return: table
    :return type: pandas.DataFrame
    """
    url = defaults.api_uniprot + identifier + ".gff"
    cols = "NAME SOURCE TYPE START END SCORE STRAND FRAME GROUP empty".split()

    data = pd.read_table(url, skiprows=2, names=cols)
    groups = data.GROUP.apply(parse_qs)
    groups = pd.DataFrame.from_records(groups)
    data = data.merge(groups, left_index=True, right_index=True)

    return data
