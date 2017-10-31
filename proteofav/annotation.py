# -*- coding: utf-8

"""
Created on 17:26 19/02/2016 2016 
Define auxiliary functions for interacting with Uniprot.
"""

import os
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
from proteofav.utils import (exclude_columns, constrain_column_types,
                             Downloader, GenericInputs)
from proteofav.library import annotation_types

log = logging.getLogger('proteofav.config')

__all__ = ["parse_gff_features", "filter_annotation",
           "select_annotation", "download_annotation", "Annotation"]


def parse_gff_features(filename, excluded_cols=None):
    """
    Map Uniprot GFF features to the protein sequence.

    :param filename: path to the Annotation file
    :param excluded_cols: option to exclude Validation columns
    :return: returns a pandas DataFrame
    """

    log.debug("Parsing GFF annotations from file...")

    header = ('NAME', 'SOURCE', 'TYPE', 'START', 'END',
              'SCORE', 'STRAND', 'FRAME', 'GROUP', 'empty')
    data = pd.read_table(filename, skiprows=2, names=header)
    groups = data.GROUP.apply(parse_qs)
    groups = pd.DataFrame.from_records(groups)

    table = data.merge(groups, left_index=True, right_index=True)

    # excluding columns
    if excluded_cols is None:
        excluded_cols = ('empty',)
    table = exclude_columns(table, excluded=excluded_cols)

    # enforce some specific column types
    table = constrain_column_types(table, col_type_dict=annotation_types)

    if table.empty:
        raise ValueError("The filters resulted in an empty DataFrame...")
    return table


def annotation_aggregation(table, identifier=None, query_type='', group_residues=True,
                           drop_types=('Helix', 'Beta strand', 'Turn', 'Chain')):
    """
    Filtering and making the final Annotation Table.

    :param table: pandas DataFrame object
    :param identifier: to add a new column with 'accession' ID
    :param str query_type: Select type of feature
    :param bool group_residues: by default each row in the resulting table,
        maps to a residue. When set to False, each row represent a feature
        per residue.
    :param tuple drop_types: Filter out some of the features, important to
        remove fetures that spam.
    :return pd.DataFrame: table. Columns will depend on parameters.
    """
    if query_type:
        table = table[table.TYPE == query_type]
    elif drop_types:
        table = table[~table.TYPE.isin(drop_types)]

    lines = []
    for i, row in table.iterrows():
        lines.extend({'idx': i, 'annotation': _annotation_writer(row)}
                     for i in range(row.START, row.END + 1))
    table = pd.DataFrame(lines)

    if group_residues:
        table = table.groupby('idx').agg({'annotation': lambda x: ', '.join(x)})

    table['site'] = table.index.astype(str)
    table['accession'] = [identifier] * len(table)

    log.debug("Aggregated annotations...")

    if table.empty:
        raise ValueError("The filters resulted in an empty DataFrame...")
    return table


def _annotation_writer(gff_row):
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


def filter_annotation(table, identifier=None, annotation_agg=False, **kwargs):
    """
    Filtering and making the final Annotation Table.

    :param table: pandas DataFrame object
    :param identifier: to add a new column with 'accession' ID
    :param annotation_agg: boolean
    :return: returns a pandas DataFrame
    """

    if annotation_agg:
        table = annotation_aggregation(table, identifier=identifier, **kwargs)

    if table.empty:
        raise ValueError("The filters resulted in an empty DataFrame...")
    return table


def select_annotation(identifier, excluded_cols=None, overwrite=False, **kwargs):
    """
    Produces table from PDB validation XML file.

    :param identifier: UniProt accession ID
    :param excluded_cols: option to exclude columns
    :param overwrite: boolean
    :return: returns a pandas DataFrame
    """
    filename = os.path.join(defaults.db_annotation, "{}.gff".format(identifier))

    download_annotation(identifier=identifier, filename=filename,
                        overwrite=overwrite)

    table = parse_gff_features(filename=filename, excluded_cols=excluded_cols)
    table = filter_annotation(table, identifier, **kwargs)
    table = constrain_column_types(table, col_type_dict=annotation_types)
    return table


def download_annotation(identifier, filename, overwrite=False):
    """
    Downloads a Annotation GFF from the UniProt.

    :param identifier: (str) UniProt accession ID
    :param filename: path to the Validation file
    :param overwrite: (boolean)
    :return: (side effects) output file path
    """

    url_root = defaults.api_uniprot
    url_endpoint = "{}.gff".format(identifier)
    url = url_root + url_endpoint
    Downloader(url=url, filename=filename,
               decompress=False, overwrite=overwrite)


class Annotation(GenericInputs):
    def read(self, filename=None, **kwargs):
        filename = self._get_filename(filename)
        self.table = parse_gff_features(filename=filename, **kwargs)
        return self.table

    def download(self, identifier=None, filename=None, **kwargs):
        identifier = self._get_identifier(identifier)
        filename = self._get_filename(filename)
        return download_annotation(identifier=identifier, filename=filename, **kwargs)

    def select(self, identifier=None, **kwargs):
        identifier = self._get_identifier(identifier)
        self.table = select_annotation(identifier=identifier, **kwargs)
        return self.table


Annotation = Annotation()
