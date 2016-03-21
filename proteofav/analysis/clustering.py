#!/usr/bin/env python
# -*- coding: utf-8

"""
"Starting out with the example from
http://stackoverflow.com/questions/21638130/tutorial-for-scipy-cluster-hierarchy
"""

import csv
import os
import sys
from operator import itemgetter
from subprocess import call
from time import strftime

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.cluster.hierarchy as hac
from proteofav.main import merge_tables
from mcl.mcl_clustering import mcl
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from scipy.spatial import ConvexHull
from scipy.spatial.distance import pdist, squareform, euclidean
from scipy.spatial.qhull import QhullError

from proteofav.analysis.random_annotations import add_random_disease_variants
from proteofav.analysis.utils import _get_colors, _autoscale_axes, _fractional_to_cartesian
from proteofav.analysis.utils import delete_file, _get_colors, _fractional_to_cartesian, \
    _autoscale_axes

__author__ = 'smacgowan'


##############################################################################
# Distance Matrix Helpers
##############################################################################
def atom_dist(table, mask, cartesian=False):
    """
    Find the intervariant distances for a given PDB chain.

    :param table:
    :param cartesian:
    :return: The reduced distance matrix produced by `pdist` and the XYZ
     coordinates of the mapped variants.
    """
    included = table[mask]
    excluded = table[~mask]

    # Build array distance matrix
    included, included_xyz = select_valid_coordinates(included)
    excluded, excluded_xyz = select_valid_coordinates(excluded)

    # Convert to fractional coordinates
    if cartesian:
        pdb_id = table.PDB_dbAccessionId.unique()[0]
        included_xyz = _fractional_to_cartesian(included_xyz, pdb_id)
        excluded_xyz = _fractional_to_cartesian(excluded_xyz, pdb_id)

    d = pdist(included_xyz)

    return d, included_xyz, included[['UniProt_dbResNum', 'chain_id']], excluded_xyz


def select_valid_coordinates(table):
    """
    Get coordinates out of a structure table.

    :param table: A ProteoFAV structure table.
    :return: A tuple containing a filtered structure table [0] and the coordinate array [1]
    """
    real_cartn_mask = np.isfinite(table['Cartn_x'])
    table = table[real_cartn_mask]
    coord_array = np.array([table.Cartn_x, table.Cartn_y, table.Cartn_z])
    coord_array = np.transpose(coord_array)
    return table, coord_array


def dist_to_sim(d, method='max_minus_d', threshold=float('inf'), gamma=1./3, **kwargs):
    """
    Covert a scipy distance matrix into a similarity matrix.

    :param d: A reduced distance matrix, such as produced by `pdist` or the
     expanded, squarematrix form.
    :param method: The distance-to-similarity conversion to employ before MCL analysis.
      Choose from: 'max_minus_d', 'reciprocal', 'exp_decay' or 'alt_reciprocal'
    :param threshold: The threshold to discard distances (i.e. remove edges)
     before conversion to similarities
    :param gamma: The gamma coefficient in the exponential decay transform (s = e^(-d * gamma)
    :return:

    NB. **kwargs is just so that this can be called with extraneous arguments.
    """
    if method == 'max_minus_d':
        d[d > threshold] = d.max()
        s = d.max() - d

    if method == 'reciprocal':
        d = d + 1
        s = 1. / d

    if method == 'exp_decay':
        if gamma <= 0:
            raise ValueError('Invalid gamma: must be >= 0')
        s = np.exp(-d * gamma)
        s[d > threshold] = 0

    if method == 'alt_reciprocal':
        ##d[d > threshold] = float('inf')  #TODO: Any threshold here hides cluster structure. Why?
        s = 1. / (d / d.max())

    return s


##############################################################################
# MCL Program Interface
##############################################################################
def write_mcl_input(s):
    """
    Re-format a similarity matrix into a list of edge weights and write to a file.

    :param s: A similarity matrix
    :return: Writes a CSV file in the current directory
    """
    with open('mcl_interactions.csv', 'wb') as file:
        writer = csv.writer(file, delimiter=' ')
        length = len(s)
        for i in xrange(length):
            for j in xrange(length):
                if s[i, j] != 0.:
                    writer.writerow([i, j, s[i, j]])


def read_mcl_clusters(file='mcl_results.txt'):
    """
    Parse a results file produced by the MCL program.

    :param file: The name of an MCL results file
    :return: A list containing the nodes found in each cluster
    """
    with open(file, 'rb') as csvfile:
        reader = csv.reader(csvfile, delimiter='\t')
        clusters = [row for row in reader]

    return clusters


def launch_mcl(s, format='partition', inflate=None, silent=True, **kwargs):
    """
    Perform MCL analysis using an external MCL implementation.

    :param s: A similarity matrix
    :param format: 'partition' returns the clusters in the same format produced
      by `fcluster` on a linkage object;
    any other value leaves the clusters in the format provided by the MCL program
    :return: The clusters found by the MCL program

    NB. **kwargs is just so that this can be called with extraneous arguments.
    """
    # Write the MCL input graph to a file
    write_mcl_input(s)

    # Build and submit the call
    mcl_call = ['mcl', 'mcl_interactions.csv', '--abc', '-o', 'mcl_results.txt']
    if inflate:
        mcl_call.append('-I')
        mcl_call.append(str(inflate))
    if silent:
        call(mcl_call, stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
    else:
        call(mcl_call)

    # Read results
    clusters = read_mcl_clusters('mcl_results.txt')

    # Cleanup
    delete_file('mcl_results.txt')
    delete_file('mcl_interactions.csv')

    if format == 'partition':
        part = ['unassigned'] * len(s)
        cluster_id = 0
        for i in clusters:
            for j in i:
                part[int(j)] = cluster_id
            cluster_id += 1
        return part
    else:
        return clusters


##############################################################################
# Cluster results formatting
##############################################################################
def cluster_dict_to_partitions(cluster_dict):
    """
    Convert the cluster output format from the MCL function to a partition list;
     like that produced by `fcluster`

    :param cluster_dict: A dictionary where the keys indicate cluster_dict and the
     list elements indicate node membership
    :return: A list where the positioning indexes correspond to node numbers and the
     values correspond to its cluster
    """
    new_dict = {e: k for k, v in cluster_dict.items() for e in v}
    part = [new_dict[k] + 1 for k in sorted(new_dict.keys())]
    part = np.array(part)
    return part


def merge_clusters_to_table(residue_ids, partition, target_table, multichain=False):
    """
    Annotate a structure table with cluster membership.

    :param residue_ids:
    :param partition:
    :param target_table:
    :param multichain:
    :return:
    """
    df = pd.DataFrame(residue_ids)
    df['Cluster'] = pd.Series(partition)
    target_table.UniProt_dbResNum = target_table.UniProt_dbResNum.astype('float')
    if multichain:
        merge_columns = ['UniProt_dbResNum', 'chain_id']
    else:
        merge_columns = 'UniProt_dbResNum'
    merged = pd.merge(target_table, df, on=merge_columns)
    return merged


##############################################################################
# Comparing clustering routines
##############################################################################
def linkage_cluster(d, methods=['single', 'complete'], similarity='max_minus_d',
                    **kwargs):
    """

    :param d: A reduced distance matrix, such as produced by `pdist`
    :param methods: A list of strings indicating the cluster methods to employ.
      Choose from: 'single', 'complete', 'average' and 'mcl'
    :param invert_method: See `dist_to_sim`
    :param threshold: See `dist_to_sim`
    :return:
    """
    linkages = []
    for method in methods:
        if not method.startswith('mcl'):
            z = hac.linkage(d, method=method)
            linkages.append([z, method])
        else:
            sq = squareform(d)
            s = dist_to_sim(sq, method=similarity, **kwargs)
            if method.startswith('mcl_program'):
                clusters = launch_mcl(s, **kwargs)
                linkages.append([clusters, method])
            else:
                if not 'inflate' in kwargs:
                    inflate_factor = 2
                else:
                    inflate_factor = kwargs['inflate']
                M, clusters = mcl(s, max_loop=50, inflate_factor=inflate_factor)
                linkages.append([[M, clusters], method])

    return linkages


def compare_clustering(linkages, xyz, title=None, addn_points=None):
    """
    Generate a plots comparing the cluster analyses provided.

    :param linkages: A list of cluster analyses as produced by `linkage_cluster`
    :param xyz: The coordinates of the original clustered points
    :return:
    """
    nrows = len(linkages)
    fig, axes23 = plt.subplots(nrows, 3)
    links, methods = zip(*linkages)
    offsets = [0, 3, 6, 9][:nrows]

    for z, method, axes, offset in zip(links, methods, axes23, offsets):

        # Plotting
        if not method.startswith('mcl'):
            axes[0].plot(range(1, len(z) + 1), z[::-1, 2])
            knee = np.diff(z[::-1, 2], 2)
            axes[0].plot(range(2, len(z)), knee)

            num_clust1 = knee.argmax() + 2
            knee[knee.argmax()] = 0
            num_clust2 = knee.argmax() + 2

            axes[0].text(num_clust1, z[::-1, 2][num_clust1 - 1], 'possible\n<- knee point')

            part1 = hac.fcluster(z, num_clust1, 'maxclust')
            part2 = hac.fcluster(z, num_clust2, 'maxclust')

        elif method.startswith('mcl_program'):
            part1 = np.array(z)
            num_clust1 = max(part1) + 1
            part2 = part1
            num_clust2 = num_clust1

        else:
            # Convert the MCL clusters into a partition vector
            part1 = cluster_dict_to_partitions(z[1])
            num_clust1 = max(part1)
            part2 = part1
            num_clust2 = num_clust1

        clr = _get_colors(max([num_clust1, num_clust2]))

        for part, i in zip([part1, part2], [2, 3]):
            ax = fig.add_subplot(nrows, 3, i + offset, projection='3d')
            for cluster in set(part):
                ax.scatter(xyz[part == cluster, 0], xyz[part == cluster, 1],
                           xyz[part == cluster, 2], c=clr[cluster - 1])
                # Compute and plot the Point Cloud Complex Hull
                points = xyz[part == cluster]
                if len(points) >= 4:
                    try:
                        hull = ConvexHull(points)
                        for simplex in hull.simplices:
                            simplex = np.append(simplex, simplex[0])  # Closes facet
                            ax.plot(points[simplex, 0], points[simplex, 1],
                                    points[simplex, 2], color=clr[cluster - 1])
                            tri = Poly3DCollection([zip(points[simplex, 0], points[simplex, 1],
                                                        points[simplex, 2])],
                                                   alpha=0.2)
                            tri.set_color(clr[cluster - 1])
                            tri.set_edgecolor('k')
                            ax.add_collection3d(tri)
                    except QhullError:
                        pass
                else:
                    points = np.vstack([points, points[0]])  # Gives triangle for 3-points
                    ax.plot(points[:, 0], points[:, 1], points[:, 2], color=clr[cluster - 1])
            if addn_points is not None:
                ax.scatter(addn_points[:, 0], addn_points[:, 1], addn_points[:, 2], c='grey', s=5,
                           alpha=0.3)
                x = np.append(xyz[:, 0], addn_points[:, 0])
                y = np.append(xyz[:, 1], addn_points[:, 1])
                z = np.append(xyz[:, 2], addn_points[:, 2])
                all_points = np.array([x, y, z]).T
                x, y, z = _autoscale_axes(all_points)
                ax.auto_scale_xyz(x, y, z)
            else:
                x, y, z = _autoscale_axes(xyz)
                ax.auto_scale_xyz(x, y, z)

        m = '\n(method: {})'.format(method)
        plt.setp(axes[0], title='Screeplot{}'.format(m), xlabel='partition',
                 ylabel='{}\ncluster distance'.format(m))
        plt.setp(axes[1], title='{} Clusters'.format(num_clust1))
        plt.setp(axes[2], title='{} Clusters'.format(num_clust2))
        axes[1].axis('off')
        axes[2].axis('off')

    plt.tight_layout()
    plt.suptitle(title)
    file = 'cluster_figs/cluster_fig_' + strftime("%Y%m%d_%H%M%S") + '.png'
    plt.savefig(file, format='png')
    # plt.show()


##############################################################################
# Cluster statistic functions
##############################################################################
def n_clusters(partition):
    """
    Determine the number of clusters from a partition list.

    :param part: A partition list
    :type part: List
    :return: The number of clusters in the partition list
    """
    part = list(partition)  # Copy to avoid sorting outside scope
    part.sort()

    if part[0] == 0:
        return part[-1] + 1
    else:
        return part[-1]


def partition_to_sizes(part):
    """
    Determine the number of elements in each cluster from a partition list.

    :param part: A partition list
    :type part: List
    :return:
    """
    offset = 1  # Need this if 1-indexed
    if 0 in part:
        offset = 0

    sizes = [part.count(i + offset) for i in xrange(n_clusters(part))]

    return sizes


def top_k_clusters(sizes, k=2):
    """
    Get the number of sites clustered by the k biggest clusters.

    :param part:
    :return:
    """
    sizes.sort(reverse=True)
    n_members = sum(sizes[:k])

    return n_members


def n_isolated(sizes):
    """
    Get the number of singleton clusters.

    :param sizes:
    :return:
    """
    return sizes.count(1)


def n_50_clusters(sizes):
    return n_x_clusters(sizes, 50)


def n_x_clusters(sizes, percent):
    """
    Get the number of clusters required to account for a given percentage of residues.

    :param sizes:
    :return: Minimum number of clusters that contains 'percent'% the observations.
    """
    n_required = sum(sizes) / (100. / percent)  ## TODO: Not perfect for odd numbers
    sizes = list(sizes)
    sizes.sort()
    n_obs = 0
    n_clusters = 0
    while n_obs < n_required:
        n_obs += sizes.pop()
        n_clusters += 1
    return n_clusters


def cluster_size_stats(part, statistics=(np.mean, np.median, np.std, min, max, len,
                                         top_k_clusters, n_isolated, n_50_clusters),
                       names=False):
    """
    Calculate a bunch of statistics for the clusters sizes in a particular clustering.

    :param part:
    :param statistics:
    :return:
    """
    sizes = partition_to_sizes(part)
    results = [stat(sizes) for stat in statistics]
    if not names:
        return results
    else:
        stat_functions = [stat.__name__ for stat in statistics]
        results = {k: v for k, v in zip(stat_functions, results)}
        return results


def centroids(partition, observations):
    """
    Return cluster centroids.

    :param part: Partition vector
    :type part: List
    :param obs: Observation array
    :type obs: NumPy array
    :return: A dictionary containing the cluster centroids.
    """
    cluster_ids = list(set(partition))
    cluster_ids.sort()
    cluster_centroids = {}
    for cid in cluster_ids:
        members = member_indexes(partition, cid)
        member_observations = observations[members, :]
        cluster_centroids[cid] = np.mean(member_observations, axis=0)
    return cluster_centroids


def member_indexes(partition, cluster_id):
    """
    Get the partition list indexes for residues in a given cluster.

    :param partition:
    :param cluster_id:
    :return:
    """
    return [i for i, x in enumerate(partition) if x == cluster_id]


def davies_bouldin(partition, observations):
    """
    Return the Davies-Bouldin index for the given clustering.

    :param part: Partition vector
    :type part: List
    :param obs: Observation array
    :type obs: NumPy array
    :return: The Davies-Bouldin index
    """

    n = n_clusters(partition)
    cluster_ids = list(set(partition))
    cluster_ids.sort()
    cluster_centroids = centroids(partition, observations)
    terms = []
    for i in cluster_ids:
        ci = cluster_centroids[i]
        i_indexes = member_indexes(partition, i)
        sigma_i = np.mean([euclidean(u, ci) for u in observations[i_indexes, :]])
        Sij = []
        for j in cluster_ids:
            if i != j:
                cj = cluster_centroids[j]
                d = euclidean(ci, cj)
                j_indexes = member_indexes(partition, j)
                sigma_j = np.mean([euclidean(u, cj) for u in observations[j_indexes, :]])
                Sij.append((sigma_i + sigma_j) / d)
        if len(Sij) > 0:
            terms.append(max(Sij))
    return sum(terms) / n


def dunn(partition, observations):
    """
    Return the Dunn index for a clustering.

    :param part: Partition vector
    :type part: List
    :param obs: Observation array
    :type obs: NumPy array
    :return: The Dunn index
    """
    cluster_ids = list(set(partition))

    # NA if only one cluster
    if len(cluster_ids) == 1:
        return np.nan

    cluster_ids.sort()
    min_inter = min(pdist(centroids(partition, observations).values()))
    intras = []
    for cid in cluster_ids:
        member_ids = member_indexes(partition, cid)
        d = pdist(observations[member_ids, :])
        if len(d) != 0:  ## Handles singleton clusters
            intras.append(max(d))

    if len(intras) == 0:
        return np.nan
    else:
        max_intra = max(intras)

    return(min_inter / max_intra)


##############################################################################
# Cluster geometry
##############################################################################
def cluster_hull_volumes(partition, points):
    """
    Calculate the volumes of the convex hulls defined by a clustering.

    :param partition:
    :param points:
    :return:

    See http://stackoverflow.com/questions/24733185/volume-of-convex-hull-with-qhull-from-scipy
    for some discussion (and some of the code used here).
    """
    def tetrahedron_volume(a, b, c, d):
        return np.abs(np.einsum('ij,ij->i', a-d, np.cross(b-d, c-d))) / 6

    def convex_hull_volume_bis(pts):
        ch = ConvexHull(pts)

        simplices = np.column_stack((np.repeat(ch.vertices[0], ch.nsimplex),
                                     ch.simplices))
        tets = ch.points[simplices]
        return np.sum(tetrahedron_volume(tets[:, 0], tets[:, 1],
                                         tets[:, 2], tets[:, 3]))

    cluster_ids = list(set(partition))
    cluster_ids.sort()
    volumes = {}
    for cid in cluster_ids:
        member_ids = member_indexes(partition, cid)
        member_pts = points[member_ids, :]
        if len(member_pts) >= 4:  ## Ignore 'plane' clusters, TODO: handle duplicates
            volume = convex_hull_volume_bis(member_pts)
        else:
            volume = 0
        volumes[cid] = 0 if volume is None else volume
    return volumes


##############################################################################
# Cluster bootstrapping
##############################################################################
def bootstrap(table, methods, n_residues, n_phenotypes,
              samples=10, **kwargs):
    """
    Cluster a random selection of residues a given number of times and return all the cluster partitions.

    :param table:
    :param methods:
    :param statistics:
    :return:
    """
    partitions = []
    for i in xrange(samples):
        # Build random data
        test_table = add_random_disease_variants(table, n_residues, n_phenotypes)

        # Compute the clustering
        d, points, resids, unmapped_points = atom_dist(test_table)
        links = linkage_cluster(d, methods, **kwargs)
        partitions.append(links[0][0])

    return partitions


def bootstrap_stats(partitions, **kwargs):
    """
    Collect cluster size stats from a list of cluster partitions.

    :param partitions:
    :param statistics:
    :return:
    """
    stats = [cluster_size_stats(part, **kwargs) for part in partitions]

    return stats


def boot_pvalue(sample_stats, test):
    """
    Calculate a p-Value for observing various cluster statistics against a random sample.

    :param sample_stats:
    :param test:
    :return:
    """
    n_samples = len(sample_stats)
    p_value = sum(map(test, sample_stats)) / float(n_samples)
    p_value = round(p_value * n_samples) / n_samples
    if p_value == 0.:
        p_value = '< ' + str(1. / n_samples)

    return p_value


def tail_thresholds(alpha, samples, stats):
    """

    :param alpha:
    :param samples:
    :param stats:
    :return:
    """
    thresholds = []
    for i in zip(*stats):
        i = np.sort(i)
        thresholds.append((i[int((alpha / 2.0) * samples)], i[int((1 - alpha / 2.0) * samples)]))

    return thresholds


def plot_sample_distributions(results, names):
    for i in xrange(len(names)):
        plt.subplot(3, 4, i + 1)
        plt.title(names[i])
        data = zip(*results['sample_stats'])[i]
        if np.nan in data:
            data = [x for x in data if not np.isnan(x)]
        if all(isinstance(x, int) for x in data):
            plt.hist(data, color='c', bins=range(min(data), max(data) + 1, 1))
        else:
            plt.hist(data, color='c')
        plt.axvline(results['obs_stats'][names[i]], color='b', linestyle='dashed', linewidth=2)
        plt.xticks(rotation=45)


##############################################################################
# Convenience wrappers
##############################################################################
def cluster_table(table, mask, method, sites_only=True, similarity='max_minus_d',
                  mask_column='resn', **kwargs):
    """
    Cluster certain residues in a structure table.

    :param table:
    :param mask:
    :param test_significance:
    :param kwargs:
    :return:
    """

    # Create string describing clustering for metadata
    cluster_parameters_whitelist = ['method', 'threshold', 'similarity']
    cluster_parameters = []
    for k, v in locals().iteritems():
        if k in cluster_parameters_whitelist:
            cluster_parameters.append(str(k))
            cluster_parameters.append(str(v))
    cluster_meta_tag = '_'.join(cluster_parameters)

    # # Apply mask
    # table = table[mask]
    # mask = np.array([True] * len(table))  #TODO: Remove the need for this hack needed for `atom_dist`

    # Dedupe table by UniProt and chain if only clustering sites
    if sites_only:
        table = table.drop_duplicates(subset=['UniProt_dbResNum', 'chain_id'])
        # mask = np.array([True] * len(table))  #TODO: See above

    # Perform clustering
    d, points, resids, unmapped_points = atom_dist(table, table[mask_column].notnull())
    links = linkage_cluster(d, methods=method, similarity=similarity, **kwargs)
    part = links[0][0]

    # Format the results
    labelled_points = add_clusters_to_points(part, points)
    annotated_table = add_clusters_to_table(labelled_points, table)
    annotated_table.loc[annotated_table.cluster_id.notnull(), 'cluster_info'] = cluster_meta_tag  ## New column until pandas reliably stores metadata

    return annotated_table


def test_cluster_significance(test_table, method, similarity, table, show_progress, n_samples, return_samples,
                              mask_column='resn', **kwargs):
    """
    Test the significance of the clusters in a clustered structure table by bootstrapping random residue selections.

    :param test_table:
    :param method:
    :param similarity:
    :param table:
    :param show_progress:
    :param n_samples:
    :param return_samples:
    :param mask_column:
    :param kwargs:
    :return:
    """
    n_variants = sum(table[mask_column].notnull() & np.isfinite(table['Cartn_x']))

    # Diferent columns are present depending on what variant annotation was used
    if mask_column == 'resn':
        n_phenotypes = len(table.disease.dropna().unique())
        variant_columns = ['resn', 'mut', 'disease']
    elif mask_column =='variant_id':
        n_phenotypes = 1
        variant_columns = ['translation', 'variant_id', 'start', 'residues', 'from_aa', 'to_aa']
    clean_table = table.drop_duplicates(subset=['UniProt_dbResNum', 'chain_id']).drop(variant_columns,
                                                                                      axis=1)
    clean_table = clean_table[np.isfinite(clean_table['Cartn_x'])]
    # Bootstrap samples
    annotated_tables = bootstrap_residue_clusters(clean_table, method,
                                                  n_phenotypes,
                                                  n_samples,
                                                  n_variants,
                                                  show_progress,
                                                  similarity,
                                                  **kwargs)
    # Drop unneccesary data
    column_whitelist = ['Cartn_x', 'Cartn_y', 'Cartn_z', 'cluster_id']
    annotated_tables[:] = [i.loc[i.cluster_id.notnull(), column_whitelist] for i in annotated_tables]

    part, points = clustered_table_to_partition_and_points(test_table)
    bs_stats, p_values, sample_cluster_sizes, stats = collect_cluster_sample_statistics(part, points,
                                                                                        annotated_tables)

    results = {'p': p_values, 'obs_stats': stats, 'sample_stats': bs_stats}
    if return_samples:
        results.update({'sample_size_dist': sample_cluster_sizes, 'random_samples': annotated_tables})

    return results


def collect_cluster_sample_statistics(test_part, test_points, sample_tables):
    """

    :param test_part:
    :param test_points:
    :param sample_tables:
    :return:
    """
    # Parse the random sample cluster tables
    parsed_samples = [clustered_table_to_partition_and_points(i) for i in sample_tables]

    # Gather the random sample cluster statistics
    sample_davies_bouldins = [davies_bouldin(part, points) for part, points in parsed_samples]
    sample_dunns = [dunn(part, points) for part, points in parsed_samples]
    sample_largest_cluster_volume = [max(cluster_hull_volumes(part, points)) for part, points in parsed_samples]
    random_sample_cluster_statistics = bootstrap_stats(zip(*parsed_samples)[0], names=True)  # Partition metrics

    # Combine the random sample cluster statistics
    for i in xrange(len(sample_tables)):
        random_sample_cluster_statistics[i].update({'Davies-Bouldin': sample_davies_bouldins[i]})
        random_sample_cluster_statistics[i].update({'Dunn_index': sample_dunns[i]})
        random_sample_cluster_statistics[i].update({'largest_cluster_volume': sample_largest_cluster_volume[i]})

    # Collect the observed cluster statistics
    observed_stats = {}
    observed_stats.update(cluster_size_stats(test_part, names=True))
    observed_stats.update(cluster_spatial_statistics(test_part, test_points))

    # Calculate p for randomly seeing a value smaller, equal to or larger than observed statistic
    p_values = {}
    for stat_key in observed_stats:
        sampled = map(itemgetter(stat_key), random_sample_cluster_statistics)
        obs = observed_stats[stat_key]
        left = boot_pvalue(sampled, lambda x: x < obs)
        mid = boot_pvalue(sampled, lambda x: x == obs)
        right = boot_pvalue(sampled, lambda x: x > obs)
        p_values.update({stat_key: (left, mid, right)})

    ## For complete cluster size distribution
    sample_cluster_sizes = []
    for i in zip(*parsed_samples)[0]:
        sample_cluster_sizes.append(partition_to_sizes(i))
    sample_cluster_sizes = [e for sublist in sample_cluster_sizes for e in sublist]  # Flatten

    return random_sample_cluster_statistics, p_values, sample_cluster_sizes, observed_stats


def cluster_spatial_statistics(test_part, test_points):
    """

    :param test_part:
    :param test_points:
    :return:
    """
    obs_stats = []
    obs_stats.append(davies_bouldin(test_part, test_points))
    obs_stats.append(dunn(test_part, test_points))
    obs_stats.append(max(cluster_hull_volumes(test_part, test_points).values()))
    names = []
    names.append('Davies-Bouldin')
    names.append('Dunn_index')
    names.append('largest_cluster_volume')
    results = {k: v for k, v in zip(names, obs_stats)}
    return results


##############################################################################
# Bootstrapping
##############################################################################
def bootstrap_residue_clusters(clean_table, method, n_phenotypes, n_samples, n_variants, show_progress, similarity,
                               **kwargs):
    """

    :param clean_table: A structure table
    :param method: Method to use for clustering.
    :param n_phenotypes: Number of phenotypes for `add_random_disease_variants`.
    :param n_samples: Number of random samples to produce and assess.
    :param n_variants: Number of variants for `add_random_disease_variants`.
    :param show_progress: Show a progress bar.
    :param kwargs: Arguments passed to `cluster_table`.
    :return:
    """
    annotated_tables = []
    for i in xrange(n_samples):
        sample_table = add_random_disease_variants(clean_table, n_variants, n_phenotypes)
        sample_mask = sample_table.resn.notnull()
        annotated_table = cluster_table(sample_table, sample_mask, method, n_samples=0,
                                        similarity=similarity, **kwargs)
        annotated_tables.append(annotated_table)
        if show_progress:
            pc_complete = (i + 1) / float(n_samples) * 100
            if pc_complete % 10 == 0:
                sys.stdout.write("\r%d%%" % (pc_complete))
                sys.stdout.flush()

    return annotated_tables


##############################################################################
# Data handling
##############################################################################
def add_clusters_to_points(cluster_partition, clustered_points):
    """

    :param cluster_partition:
    :param clustered_points:
    :return:
    """
    points = pd.DataFrame(clustered_points)
    points.loc[:, 'cluster_id'] = pd.Series(cluster_partition, index=points.index)
    points = points.rename(columns={0: 'Cartn_x', 1: 'Cartn_y', 2: 'Cartn_z'})
    return points


def add_clusters_to_table(labelled_points, clustered_table):
    """

    :param cluster_partition:
    :param clustered_table:
    :return:
    """
    index_name = clustered_table.index.name
    if index_name is not None:
        merged_table = clustered_table.reset_index().merge(labelled_points, how='left').set_index(index_name)
    else:
        merged_table = clustered_table.merge(labelled_points, how='left')
    return merged_table


def clustered_table_to_partition_and_points(table):
    """

    :param table: A structure table with clustered sites.
    :return: A 2-tuple with the partition vector (list) and point positions (np.array).
    """
    partition = list(pd.Series(table.cluster_id.dropna(), dtype=int))
    points = np.array(table.loc[table.cluster_id.notnull(), ['Cartn_x', 'Cartn_y', 'Cartn_z']])
    return partition, points


if __name__ == '__main__':
    # Porphobilinogen deaminase example
    table = merge_tables(pdb_id='3ecr', chain='A', uniprot_variants=True)
    mask = table.resn.notnull()
    d, points, resids, unmapped_points = atom_dist(table, mask)
    links = linkage_cluster(d, methods=['single', 'complete'])
    compare_clustering(links, points, '3ecr(a) P08397')
    links = linkage_cluster(d, methods=['average', 'mcl_program'], threshold=10)
    compare_clustering(links, points, '3ecr(a) P08397')
    links = linkage_cluster(d, methods=['mcl_program', 'mcl'], threshold=10)
    compare_clustering(links, points, '3ecr(a) P08397', addn_points=unmapped_points)

    # Comparing MCL inflation factors and distance to similarity conversion
    sq = squareform(d)
    links = []
    for inf_fact in [2., 6.]:
        s = dist_to_sim(sq, method='max_minus_d', threshold=10)
        links.append([mcl(s, max_loop=50, inflate_factor=inf_fact),
                      'mcl_IF=' + str(inf_fact) + '\nmax_minus_d(t=10)'])
    compare_clustering(links, points, '3ecr(a) P08397')

    # No min. distance for connected threshold (saturated network)
    links = []
    for inf_fact in [2., 6.]:
        s = dist_to_sim(sq, method='max_minus_d')
        links.append([mcl(s, max_loop=50, inflate_factor=inf_fact),
                      'mcl_IF=' + str(inf_fact) + '\nmax_minus_d(t=inf)'])
    compare_clustering(links, points, '3ecr(a) P08397')

    # Using reciprocal distance for similarity
    links = []
    for inf_fact in [2., 6.]:
        s = dist_to_sim(sq, method='reciprocal')
        links.append([mcl(s, max_loop=50, inflate_factor=inf_fact),
                      'mcl_IF=' + str(inf_fact) + '\nreciprocal'])
    compare_clustering(links, points, '3ecr(a) P08397')

    # KRT14 from K5/14 dimer example (multichain)
    table = merge_tables(pdb_id='3tnu', chain='all', uniprot_variants=True)
    mask = table.resn.notnull()
    d, points, resids, unmapped_points = atom_dist(table, mask)
    links = linkage_cluster(d, methods=['average', 'mcl_program'], threshold=15)
    compare_clustering(links, points, '3tnu(a/b) P02533/P13647', addn_points=unmapped_points)

    # Now look at the same structure with random variants added
    n_variants = sum(table.resn.notnull())
    table = merge_tables(pdb_id='3tnu', chain='all')
    table = add_random_disease_variants(table, n_variants, 1)
    mask = table.resn.notnull()
    d, points, resids, unmapped_points = atom_dist(table, mask)
    links = linkage_cluster(d, methods=['average', 'mcl_program'], threshold=15)
    compare_clustering(links, points, '3tnu(a/b) P02533/P13647', addn_points=unmapped_points)

    # Serine/threonine-protein kinase receptor R3, Telangiectasia example
    table = merge_tables(pdb_id='3my0', chain='A', uniprot_variants=True)
    mask = table.resn.notnull()
    d, points, resids, unmapped_points = atom_dist(table, mask)
    links = linkage_cluster(d, methods=['average', 'mcl_program'], threshold=10)
    compare_clustering(links, points, '3my0(a) P37023', addn_points=unmapped_points)

    # Alpha-galactosidase A, Fabry disease example
    table = merge_tables(pdb_id='3s5z', chain='A', uniprot_variants=True)
    mask = table.resn.notnull()
    d, points, resids, unmapped_points = atom_dist(table, mask)
    links = linkage_cluster(d, methods=['average', 'mcl_program'], threshold=10)
    compare_clustering(links, points, '3s5z(a) P06280', addn_points=unmapped_points)

    # Cholinesterase, BChE deficiency example
    table = merge_tables(pdb_id='4tpk', chain='A', uniprot_variants=True)
    mask = table.resn.notnull()
    d, points, resids, unmapped_points = atom_dist(table, mask)
    links = linkage_cluster(d, methods=['average', 'mcl_program'], threshold=10)
    compare_clustering(links, points, '4tpk(a) P06276', addn_points=unmapped_points)

    # UDP-glucose 4-epimerase, EDG example
    table = merge_tables(pdb_id='1ek6', chain='A', uniprot_variants=True)
    mask = table.resn.notnull()
    d, points, resids, unmapped_points = atom_dist(table, mask)
    links = linkage_cluster(d, methods=['average', 'mcl_program'], threshold=10)
    compare_clustering(links, points, '1ek6(a) Q14376', addn_points=unmapped_points)
