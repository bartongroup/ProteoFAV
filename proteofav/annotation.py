# -*- coding: utf-8

"""
Created on 17:26 19/02/2016 2016 
Define auxiliary functions for interacting with Uniprot.
"""

import logging

import pandas as pd

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

from proteofav.config import defaults

log = logging.getLogger('proteofav.config')

__all__ = ["_fetch_uniprot_gff",
           "map_gff_features_to_sequence"]


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
