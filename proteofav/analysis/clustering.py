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
from scipy.cluster.hierarchy import fcluster
from scipy.cluster.hierarchy import linkage
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


def cluster_significance(uniprot_id, cif_path=None, column='somatic_variants'):
    # cluster significance experiment
    # variants_cols=None cluster which variants
    # n_samples=10 number of times to bootstrap
    # project_name where to save the files

    # create a directory

    # load the structure table and variants

    # TODO use atoms=None and centroid position
    """

    :param str or None cif_path: path for external cif file
    :param str uniprot_id: uniprot id
    :param str column: location to be clustered
    :return:
    """
    if cif_path is not None:
        table = proteofav.merge_tables(cif_path=cif_path, atoms='centroid', ignore_DSSP=True)
    else:
        table = proteofav.merge_tables(uniprot_id=uniprot_id)
    variants = proteofav.variants.select_variants(uniprot_id, column)
    table['UniProt_dbResNum '] = table['UniProt_dbResNum'].astype(float)
    table = table.merge(variants, left_on='UniProt_dbResNum', right_index=True, how='left')

    # fall back to PDB reference
    table.set_index('PDB_dbResNum', inplace=True)

    # multiple variants in the same residue are represented as multiple nodes in the same location
    # TODO past version also supported clustering residues with/out variants
    temp = pd.DataFrame(table[column].tolist(), index=table.index)
    temp = temp.stack()
    temp.index = temp.index.droplevel(1)
    temp.name = column
    table = table.drop([column], axis=1).join(temp)

    return clustering(table, column)
    #     for i in range(n_samples):
    #         temp = variants.copy().sample(frac=1, replace=True)
    #         temp.index = variants.index


def exp_decay(distance, gamma=1):
    if gamma <= 0:
        raise ValueError('Invalid gamma: must be >= 0')
    return np.exp(-distance * gamma)


def max_minus_distance(distance, threshold=10):
    distance[distance > threshold] = distance.max()
    return distance.max() - distance


def complete_clustering(distance, method='complete', threshold=4):
    return fcluster(linkage(distance, method=method), threshold, 'distance').tolist()


def pv_clustering(distance_matrix, alpha=.95, n_boots=1000):
    """

    :param distance_matrix:
    :param alpha:
    :param n_boots:
    :return:
    """

    import rpy2.robjects.numpy2ri
    from rpy2.interactive.packages import importr

    rpy2.robjects.numpy2ri.activate()
    pv_clust = importr('pvclust')

    pv_clust_obj = pv_clust.pvclust(squareform(distance_matrix).T,
                                    method_dist='euclidean',
                                    method_hclust='complete',
                                    nboot=n_boots)
    pv_clust.pvrect(pv_clust_obj, alpha=alpha)

    clusters_r = pv_clust.pvpick(pv_clust_obj, alpha=alpha)
    return [[label[1:] for label in _] for _ in clusters_r.rx2('clusters')]


def mcl_clustering(distance, similarity_fun=max_minus_distance, **kwargs):
    """
    Apply MCL clustering on a similarity matrix.

    :param np.array similarity: similarity matrix
    :return pd.Series: series of cluster nodes. Indexed for join with Spatial Variant Table
    (index = PDB_dbResNum)
    """
    similarity = squareform(similarity_fun(distance, **kwargs))
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


def clustering(table, column, clustering_method, shuffle=False, **kwargs):
    """ Wrapper to support multiple clustering methods

    :param shuffle:
    :param kwargs: object
    :param function similarity_fun:
    :param function clustering_method: clustering method
    :param pd.DataFrame table: table from ProteoFAV.merge_tables
    :param str column: name of the column to be clustered
    :return:
    """
    if shuffle:
        table[column] = table[column].sample(frac=1, replace=True).values

    variant_3D = table.copy().loc[
        table['Cartn_x'].notnull() & table[column].notnull(),
        ['Cartn_x', 'Cartn_y', 'Cartn_z']]

    # First calculated the pairwise euclidean distance
    distance = pdist(variant_3D)
    clusters = clustering_method(distance, **kwargs)
    # this wont work if a residues participates of two clusters, so we do it reversed to keep
    # residues in their biggest cluster
    if not clusters:
        raise ValueError(' No clusters with giving paramenters')
    try:
        for i, cluster_i in enumerate(reversed(clusters)):
            cluster_i = [int(x) for x in cluster_i]
            variant_3D.loc[variant_3D.index[cluster_i], 'cluster_id'] = i
    except TypeError:
        variant_3D['cluster_id'] = clusters
    return variant_3D['cluster_id']


def write_chimera_cluster_as_sphere(table, column):
    """

    :param table: table of indexed by the residue position in the 3D structure
    :param column: column with cluster ids
    """
    for cluster_id in table[column].unique():
        if np.isnan(cluster_id) or cluster_id == -1:
            continue
        current = table[table[column] == cluster_id]
        # TODO: make sure the sphere comprises all variants with the radius parameter?
        # TODO: vary the transparency with the number of clusters
        parameters = {'radius': 10,
                      'model': 0,
                      'residues': ",".join(current.index.unique().astype(str)),
                      'color': "1,0,0,0.5"}

        print("shape sphere radius {radius} center #0:{residues} color {color} mesh true".format(
            **parameters))

if __name__ == '__main__':
    table = cluster_significance('O15294', column='ensembl_somatic')
