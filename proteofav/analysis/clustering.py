#!/usr/bin/env python
# -*- coding: utf-8
"""
Pipeline for spatial clustering of genetic variants using ProteoFAV
"""
import logging
import subprocess
from tempfile import NamedTemporaryFile

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform

import proteofav

log = logging.getLogger(__name__)
__author__ = 'smacgowan, tbrittoborges'


def generate_chimera_attrFile(column):
    """Write a chimera attribute fil
    :param pd.Series column: column with attributes

    ..example:
    generate_chimera_attrFile(variants.applymap(counts).tcga)
    """
    with open('{}.txt'.format(column.name), 'w') as open_f:
        open_f.write(
            'attribute: {}\nmatch mode: 1-to-1\nrecipient: residues\n'.format(column.name))
        open_f.write(column.reset_index()
                     .to_string(index=False, header=False,
                                formatters=["\t:{}".format, "\t{}".format]).
                     replace(" ", ""))


def counts(x):
    """
    Count elements inside column of iterable

    :param object x:
    :return int: number of elements per
    """
    try:
        return len(x)
    except TypeError:
        return 0


def cluster_significance(uniprot_id, column='somatic_variants'):
    # cluster significance experiment
    # variants_cols=None cluster which variants
    # n_samples=10 number of times to bootstrap
    # project_name where to save the files

    # create a directory

    # load the structure table and variants

    # TODO use atoms=None and centroid position
    """

    :param uniprot_id:
    :param column:
    :return:
    """
    table = proteofav.merge_tables(uniprot_id=uniprot_id)
    variants = proteofav.variants.select_variants(uniprot_id, column)
    table.UniProt_dbResNum = table.UniProt_dbResNum.astype(float)
    table = table.merge(variants, left_on='UniProt_dbResNum', right_index=True, how='left')

    # fall back to PDB reference
    table.set_index('PDB_dbResNum', inplace=True)

    # multiple variants in the same residue are represented as multiple nodes in the same location
    # TODO past version also supported clustering residues with/out variants
    temp = pd.DataFrame(table[column].tolist(),
                        index=table.index)
    temp = temp.stack()
    temp.index = temp.index.droplevel(1)
    temp.name = column
    table = table.drop([column], axis=1).join(temp)

    return clustering(table, column)
    #     for i in range(n_samples):
    #         temp = variants.copy().sample(frac=1, replace=True)
    #         temp.index = variants.index


def distance_to_similarity(distance, method='max_minus_d', threshold=10., gamma=1. / 3):
    """
    Transform distance matrix (or array) to similarity by normalisation.

    :param np.array distance: A reduced distance matrix, such as produced by `pdist` or the
     expanded, squarematrix form.
    :param str method: The distance-to-similarity conversion to employ before MCL analysis.
      Choose from: 'max_minus_d' or 'reciprocal'
    :param float threshold: The threshold to discard distances (i.e. remove edges)
     before conversion to similarities
    :param float gamma: The gamma coefficient in the exponential decay transform (s = e^(-d * gamma)
    :return:
    :raise TypeError: on unrecognised method.
    """
    if method == 'max_minus_d':
        distance[distance > threshold] = distance.max()
        return distance.max() - distance

    elif method == 'reciprocal':
        distance += 1
        return 1. / distance

    elif method == 'exp_decay':
        if gamma <= 0:
            raise ValueError('Invalid gamma: must be >= 0')
        return np.exp(-distance * gamma)

    else:
        raise TypeError('Unrecognised method parameters to distance to similarity.')


def mcl_clustering(similarity):
    """
    Apply MCL clustering on a similarity matrix.

    :param np.array similarity: similarity matrix
    :return pd.Series: series of cluster nodes. Indexed for join with Spatial Variant Table
    (index = PDB_dbResNum)
    """
    size = similarity.shape[0]
    with NamedTemporaryFile() as input_file, NamedTemporaryFile() as output_file:
        for i in xrange(size):
            for j in xrange(size):
                if similarity[i, j] != 0.:
                    input_file.write("{} {} {}\n".format(i, j, similarity[i, j]))
        try:
            # todo redirect std_out to logger.info()
            std_out = subprocess.check_output(
                ['mcl', input_file.name, '--abc', '-o', output_file.name])
        except subprocess.CalledProcessError as e:
            log.error(e, exc_info=True)

        return [line.rstrip().split('\t') for line in output_file.readlines()]


def clustering(table, column, clustering_method=mcl_clustering):
    """ Wrapper to support multiple clustering methods

    :param function clustering_method: clustering method
    :param pd.DataFrame table: table from ProteoFAV.merge_tables
    :param str column: name of the column to be clustered
    :return:
    """
    variant_3D = table.copy().loc[
        table['Cartn_x'].notnull() & table[column].notnull(),
        ['Cartn_x', 'Cartn_y', 'Cartn_z']]

    # First calculated the pairwise euclidean distance
    distance = pdist(variant_3D)
    similarity = squareform(distance_to_similarity(distance))
    clusters = clustering_method(similarity)
    # this wont work if a residues participates of two clusters, so we do it reversed to keep
    # residues in their biggest cluster
    for i, cluster_i in enumerate(reversed(clusters)):
        cluster_i = [int(x) for x in cluster_i]
        variant_3D.loc[variant_3D.index[cluster_i], 'cluster_id'] = i
    return variant_3D


if __name__ == '__main__':
    table = cluster_significance('O15294', column='ensembl_somatic')
