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

        clusters = [line.rstrip().split('\t') for line in output_file.readlines()]

    # this could not be more obfuscated
    clusters = pd.DataFrame(clusters)
    clusters = clusters.stack()
    clusters.index = clusters.index.droplevel(1)
    clusters = pd.Series(clusters.index.astype(int), index=clusters.astype(int).values)
    clusters.name = 'cluster_id'
    return clusters


def clustering(table, column):
    """

    :param pd.DataFrame table: Merged Table
    :param str column: Name of the column to be clustered
    :return:
    """
    table = table.copy()

    # First calculated the pairwise euclidean distance
    distance = pdist(table.loc[table['Cartn_x'].notnull() & table[column].notnull(),
                               ['Cartn_x', 'Cartn_y', 'Cartn_z']])
    similarity = squareform(distance_to_similarity(distance))
    clusters = mcl_clustering(similarity)
    return table.join(clusters)


def cluster_significance(uniprot_id):
    # variants_cols=None cluster which variants
    # n_samples=10 number of times to bootstrap
    # project_name where to save the files

    # create a directory


    table = proteofav.merge_tables(uniprot_id=uniprot_id)  # TODO use atoms=None and centroid posit
    variants = proteofav.variants.select_variants(uniprot_id)
    table.UniProt_dbResNum = table.UniProt_dbResNum.astype(float)
    table = table.merge(variants, left_on='UniProt_dbResNum', right_index=True, how='left')
    return table


#     for i in range(n_samples):
#         temp = variants.copy().sample(frac=1, replace=True)
#         temp.index = variants.index




if __name__ == '__main__':
    pass
