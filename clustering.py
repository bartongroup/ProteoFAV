__author__ = 'smacgowan'

# Starting out with the example from http://stackoverflow.com/questions/21638130/tutorial-for-scipy-cluster-hierarchy

import numpy as np
import scipy.cluster.hierarchy as hac
import matplotlib.pyplot as plt

from main import merge_tables
from to_table import _fetch_uniprot_variants
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial import ConvexHull
from scipy.spatial.qhull import QhullError
from mcl.mcl_clustering import mcl
from time import strftime
from utils import _get_colors
import csv
from subprocess import call


def variant_distances(pdb_id, chain, uniprot_id):
    """
    Find the intervariant distances for a given PDB chain
    :param pdb_id:
    :param chain:
    :param uniprot_id:
    :return: The reduced distance matrix produced by `pdist` and the XYZ coordinates of the mapped variants.
    """
    # Fetch raw data
    structure = merge_tables(pdb_id=pdb_id, chain=chain)  ## Don't add variants yet!
    variants = _fetch_uniprot_variants(uniprot_id)
    # Merge these
    structure.UniProt_dbResNum = structure.UniProt_dbResNum.astype('float')
    merged = pd.merge(structure, variants, on='UniProt_dbResNum')
    # Build array distance matrix
    merged_real = merged[np.isfinite(merged['Cartn_x'])]
    xyz = np.array([merged_real.Cartn_x, merged_real.Cartn_y, merged_real.Cartn_z])
    xyz = np.transpose(xyz)
    return pdist(xyz), xyz


def invert_distances(d, method, threshold=float('inf')):
    """

    :param d: A reduced distance matrix, such as produced by `pdist` or the expanded, squarematrix form.
    :param method: The distance-to-similarity conversion to employ before MCL analysis. Choose from: 'max_minus_d' or 'reciprocal'
    :param threshold: The threshold to discard distances (i.e. remove edges) before conversion to similarities
    :return:
    """
    if method == 'max_minus_d':
        d[d > threshold] = d.max()
        s = d.max() - d
    if method == 'reciprocal':
        d = d + 1
        s = 1. / d
    return s


def linkage_cluster(a, methods=['single', 'complete'], invert_method='max_minus_d', threshold=float('inf')):
    """

    :param a: A reduced distance matrix, such as produced by `pdist`
    :param methods: A list of strings indicating the cluster methods to employ. Choose from: 'single', 'complete', 'average' and 'mcl'
    :param invert_method: See `invert_distances`
    :param threshold: See `invert_distances`
    :return:
    """
    linkages = []
    for method in methods:
        if method != 'mcl':
            z = hac.linkage(a, method=method)
            linkages.append([z, method])
        else:
            a = squareform(a)
            a = invert_distances(a, invert_method, threshold)
            M, clusters = mcl(a, max_loop=50)
            linkages.append([[M, clusters], method])

    return linkages


def cluster_dict_to_partitions(clusters):
    """
    Convert the cluster output format from the MCL function to a partition list; like that produced by `fcluster`
    :param clusters: A dictionary where the keys indicate clusters and the list elements indicate node membership
    :return: A list where the positioning indexes correspond to node numbers and the values correspond to its cluster
    """
    d = {i: k for k, v in clusters.items() for i in v}
    part1 = []
    for k in sorted(d.keys()):
        part1.append(d[k] + 1)
    part1 = np.array(part1)
    return part1


def write_mcl_input(s):
    """
    Re-format a similarity matrix into a list of edge weights and write to a file.
    :param s: A similarity matrix
    :return: Writes a CSV file in the current directory
    """
    with open('mcl_interactions.csv', 'wb') as csvfile:
        writer = csv.writer(csvfile, delimiter=' ')
        length = len(s)
        for i in xrange(length):
            for j in xrange(length):
                if s[i, j] != 0.:
                    writer.writerow([i, j, s[i,j]])


def read_mcl_clusters(file='mcl_results.txt'):
    """
    Parse a results file produced by the MCL program.
    :param file: The name of an MCL results file
    :return: A list containing the nodes found in each cluster
    """
    clusters = []
    with open('mcl_results.txt', 'rb') as csvfile:
        reader = csv.reader(csvfile, delimiter='\t')
        for row in reader:
            clusters.append(row)
    return clusters


def launch_mcl(s, format='partition'):
    """
    Perform MCL analysis using an external MCL implementation.
    :param s: A similarity matrix
    :param format: 'partition' returns the clusters in the same format produced by `fcluster` on a linkage object;
    any other value leaves the clusters in the format provided by the MCL program
    :return: The clusters found by the MCL program
    """
    write_mcl_input(s)
    call(['mcl', 'mcl_interactions.csv', '--abc', '-o', 'mcl_results.txt'])
    clusters = read_mcl_clusters('mcl_results.txt')
    if format == 'partition':
        part = ['unassigned'] * len(s)
        cluster_id = 0
        for i in clusters:
            for j in i:
                part[int(j)] = cluster_id
            cluster_id += 1
        return part

    return clusters


def compare_clustering(linkages, xyz):
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
            axes[0].plot(range(1, len(z)+1), z[::-1, 2])
            knee = np.diff(z[::-1, 2], 2)
            axes[0].plot(range(2, len(z)), knee)

            num_clust1 = knee.argmax() + 2
            knee[knee.argmax()] = 0
            num_clust2 = knee.argmax() + 2

            axes[0].text(num_clust1, z[::-1, 2][num_clust1-1], 'possible\n<- knee point')

            part1 = hac.fcluster(z, num_clust1, 'maxclust')
            part2 = hac.fcluster(z, num_clust2, 'maxclust')

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
                try:
                    points = xyz[part == cluster]
                    if len(points) >= 4:
                        hull = ConvexHull(points)
                        for simplex in hull.simplices:
                            ax.plot_wireframe(points[simplex, 0], points[simplex, 1],
                                              points[simplex, 2], color=clr[cluster - 1])
                except QhullError:
                    pass

        m = '\n(method: {})'.format(method)
        plt.setp(axes[0], title='Screeplot{}'.format(m), xlabel='partition',
                 ylabel='{}\ncluster distance'.format(m))
        plt.setp(axes[1], title='{} Clusters'.format(num_clust1))
        plt.setp(axes[2], title='{} Clusters'.format(num_clust2))
        axes[1].axis('off')
        axes[2].axis('off')

    plt.tight_layout()
    file = 'cluster_figs/cluster_fig_' + strftime("%Y%m%d_%H%M%S") + '.png'
    plt.savefig(file, format='png')
    #plt.show()


if __name__ == '__main__':
    # Porphobilinogen deaminase example
    d, points = variant_distances(pdb_id='3ecr', chain='A', uniprot_id='P08397')
    links = linkage_cluster(d, methods=['single', 'complete'])
    compare_clustering(links, points)
    links = linkage_cluster(d, methods=['average', 'mcl'])
    compare_clustering(links, points)

    # Comparing MCL inflation factors and distance to similarity conversion
    sq = squareform(d)
    links = []
    for inf_fact in [2., 6.]:
        s = invert_distances(sq, method='max_minus_d', threshold=10)
        links.append([mcl(s, max_loop=50, inflate_factor=inf_fact),
                      'mcl_IF=' + str(inf_fact) + '\nmax_minus_d(t=10)'])
    compare_clustering(links, points)

    # No min. distance for connected threshold (saturated network)
    links = []
    for inf_fact in [2., 6.]:
        s = invert_distances(sq, method='max_minus_d')
        links.append([mcl(s, max_loop=50, inflate_factor=inf_fact),
                      'mcl_IF=' + str(inf_fact) + '\nmax_minus_d(t=inf)'])
    compare_clustering(links, points)

    # Using reciprocal distance for similarity
    links = []
    for inf_fact in [2., 6.]:
        s = invert_distances(sq, method='reciprocal')
        links.append([mcl(s, max_loop=50, inflate_factor=inf_fact),
                      'mcl_IF=' + str(inf_fact) + '\nreciprocal'])
    compare_clustering(links, points)


    # KRT14 from K5/14 dimer example
    d, points = variant_distances(pdb_id='3tnu', chain='A', uniprot_id='P02533')
    links = linkage_cluster(d, methods=['average', 'mcl'])
    compare_clustering(links, points)