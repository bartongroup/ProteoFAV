__author__ = 'smacgowan'

# Starting out with the example from http://stackoverflow.com/questions/21638130/tutorial-for-scipy-cluster-hierarchy

import csv
from subprocess import call
from time import strftime

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.cluster.hierarchy as hac
from mcl.mcl_clustering import mcl
from scipy.spatial import ConvexHull
from scipy.spatial.distance import pdist, squareform
from scipy.spatial.qhull import QhullError

from main import merge_tables
from to_table import _fetch_uniprot_variants
from utils import _get_colors, autoscale_axes, fractional_to_cartesian
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection


def variant_distances(pdb_id, chains, uniprot_ids, cartesian=False):
    """
    Find the intervariant distances for a given PDB chain
    :param pdb_id:
    :param chain:
    :param uniprot_id:
    :return: The reduced distance matrix produced by `pdist` and the XYZ coordinates of the mapped variants.
    """
    merged = pd.DataFrame()
    unmapped = pd.DataFrame()
    for chain, uniprot_id in zip(chains, uniprot_ids):
        # Fetch raw data
        structure = merge_tables(pdb_id=pdb_id, chain=chain)  ## Don't add variants yet!
        variants = _fetch_uniprot_variants(uniprot_id)
        # Merge these
        structure.UniProt_dbResNum = structure.UniProt_dbResNum.astype('float')
        table = pd.merge(structure, variants, on='UniProt_dbResNum', how='left')  ## TODO: Optionally uniquify variants at residue level
        merged = merged.append(table[table.resn.notnull()])
        unmapped = unmapped.append(table[table.resn.isnull()])

    # Build array distance matrix
    merged_real, xyz = extract_coords(merged)
    _, unmapped_xyz = extract_coords(unmapped)

    # Convert to fractional coordinates
    if cartesian:
        xyz = fractional_to_cartesian(xyz, pdb_id)
        unmapped_xyz = fractional_to_cartesian(unmapped_xyz, pdb_id)

    return pdist(xyz), xyz, merged_real.UniProt_dbResNum, unmapped_xyz


def extract_coords(merged):
    merged_real = merged[np.isfinite(merged['Cartn_x'])]
    xyz = np.array([merged_real.Cartn_x, merged_real.Cartn_y, merged_real.Cartn_z])
    xyz = np.transpose(xyz)
    return merged_real, xyz


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


def linkage_cluster(a, methods=['single', 'complete'], invert_method='max_minus_d', threshold=float('inf'),
                    inflate=None):
    """

    :param a: A reduced distance matrix, such as produced by `pdist`
    :param methods: A list of strings indicating the cluster methods to employ. Choose from: 'single', 'complete', 'average' and 'mcl'
    :param invert_method: See `invert_distances`
    :param threshold: See `invert_distances`
    :return:
    """
    linkages = []
    for method in methods:
        if not method.startswith('mcl'):
            z = hac.linkage(a, method=method)
            linkages.append([z, method])
        else:
            sq = squareform(a)
            s = invert_distances(sq, invert_method, threshold)
            if method.startswith('mcl_program'):
                clusters = launch_mcl(s, inflate=inflate)
                linkages.append([clusters, method])
            else:
                if not inflate:
                    inflate_factor = 2
                else:
                    inflate_factor = inflate
                M, clusters = mcl(s, max_loop=50, inflate_factor=inflate_factor)
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
                    writer.writerow([i, j, s[i, j]])


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


def launch_mcl(s, format='partition', inflate=None):
    """
    Perform MCL analysis using an external MCL implementation.
    :param s: A similarity matrix
    :param format: 'partition' returns the clusters in the same format produced by `fcluster` on a linkage object;
    any other value leaves the clusters in the format provided by the MCL program
    :return: The clusters found by the MCL program
    """
    write_mcl_input(s)
    mcl_call = ['mcl', 'mcl_interactions.csv', '--abc', '-o', 'mcl_results.txt']
    if inflate:
        mcl_call.append('-I')
        mcl_call.append(str(inflate))
    call(mcl_call)
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


def merge_clusters_to_table(residue_ids, partition, target_table):
    df = pd.DataFrame(residue_ids)
    df['Cluster'] = pd.Series(partition)
    target_table.UniProt_dbResNum = target_table.UniProt_dbResNum.astype('float')
    merged = pd.merge(target_table, df, on = 'UniProt_dbResNum')
    return merged


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
                            tri = Poly3DCollection([zip(points[simplex, 0], points[simplex, 1], points[simplex, 2])],
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
                ax.scatter(addn_points[:,0], addn_points[:,1], addn_points[:,2], c='grey', s=5, alpha=0.3)
                x = np.append(xyz[:, 0], addn_points[:,0])
                y = np.append(xyz[:, 1], addn_points[:,1])
                z = np.append(xyz[:, 2], addn_points[:,2])
                all_points = np.array([x, y, z]).T
                x, y, z = autoscale_axes(all_points)
                ax.auto_scale_xyz(x, y, z)
            else:
                x, y, z = autoscale_axes(xyz)
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


if __name__ == '__main__':
    # Porphobilinogen deaminase example
    d, points, resids, unmapped_points = variant_distances(pdb_id='3ecr', chains=['A'], uniprot_ids=['P08397'])
    links = linkage_cluster(d, methods=['single', 'complete'])
    compare_clustering(links, points, '3ecr(a) P08397')
    links = linkage_cluster(d, methods=['average', 'mcl'])
    compare_clustering(links, points, '3ecr(a) P08397')
    links = linkage_cluster(d, methods=['mcl_program', 'mcl'], threshold=10)
    compare_clustering(links, points, '3ecr(a) P08397', addn_points=unmapped_points)

    # Comparing MCL inflation factors and distance to similarity conversion
    sq = squareform(d)
    links = []
    for inf_fact in [2., 6.]:
        s = invert_distances(sq, method='max_minus_d', threshold=10)
        links.append([mcl(s, max_loop=50, inflate_factor=inf_fact),
                      'mcl_IF=' + str(inf_fact) + '\nmax_minus_d(t=10)'])
    compare_clustering(links, points, '3ecr(a) P08397')

    # No min. distance for connected threshold (saturated network)
    links = []
    for inf_fact in [2., 6.]:
        s = invert_distances(sq, method='max_minus_d')
        links.append([mcl(s, max_loop=50, inflate_factor=inf_fact),
                      'mcl_IF=' + str(inf_fact) + '\nmax_minus_d(t=inf)'])
    compare_clustering(links, points, '3ecr(a) P08397')

    # Using reciprocal distance for similarity
    links = []
    for inf_fact in [2., 6.]:
        s = invert_distances(sq, method='reciprocal')
        links.append([mcl(s, max_loop=50, inflate_factor=inf_fact),
                      'mcl_IF=' + str(inf_fact) + '\nreciprocal'])
    compare_clustering(links, points, '3ecr(a) P08397')

    # KRT14 from K5/14 dimer example (multichain)
    d, points, resids, unmapped_points = variant_distances(pdb_id='3tnu', chains=['A', 'B'], uniprot_ids=['P02533', 'P13647'])
    links = linkage_cluster(d, methods=['average', 'mcl'], threshold=10)
    compare_clustering(links, points, '3tnu(a/b) P02533/P13647', addn_points=unmapped_points)

    # Serine/threonine-protein kinase receptor R3, Telangiectasia example
    d, points, resids, unmapped_points = variant_distances(pdb_id='3my0', chains=['A'], uniprot_ids=['P37023'])
    links = linkage_cluster(d, methods=['average', 'mcl_program'], threshold=10)
    compare_clustering(links, points, '3my0(a) P37023', addn_points=unmapped_points)

    # Alpha-galactosidase A, Fabry disease example
    d, points, resids, unmapped_points = variant_distances(pdb_id='3s5z', chains=['A'], uniprot_ids=['P06280'])
    links = linkage_cluster(d, methods=['average', 'mcl_program'], threshold=10)
    compare_clustering(links, points, '3s5z(a) P06280', addn_points=unmapped_points)

    # Cholinesterase, BChE deficiency example
    d, points, resids, unmapped_points = variant_distances(pdb_id='4tpk', chains=['A'], uniprot_ids=['P06276'])
    links = linkage_cluster(d, methods=['average', 'mcl_program'], threshold=10)
    compare_clustering(links, points, '4tpk(a) P06276', addn_points=unmapped_points)

    # UDP-glucose 4-epimerase, EDG example
    d, points, resids, unmapped_points = variant_distances(pdb_id='1ek6', chains=['A'], uniprot_ids=['Q14376'])
    links = linkage_cluster(d, methods=['average', 'mcl_program'], threshold=10)
    compare_clustering(links, points, '1ek6(a) Q14376', addn_points=unmapped_points)