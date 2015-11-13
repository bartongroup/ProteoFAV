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


def variant_distances(pdb_id, chain, uniprot_id):
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
    if method == 'max_minus_d':
        d[d > threshold] = max(d)
        s = max(d) - d
    if method == 'reciprocal':
        try:
            s = 1. / d
        except ZeroDivisionError:
            d = d + 1
            s = 1. / d
    return s


def linkage_cluster(a, methods=['single', 'complete'], invert_method='max_minus_d', threshold=float('inf')):

    linkages = []
    for method in methods:
        if method != 'mcl':
            z = hac.linkage(a, method=method)
            linkages.append([z, method])
        else:
            a = invert_distances(a, invert_method, threshold)
            M, clusters = mcl(squareform(a), max_loop=50)
            linkages.append([[M, clusters], method])

    return linkages


def compare_clustering(linkages, xyz):

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
            d = {i: k for k, v in z[1].items() for i in v}
            part1 = []
            for k in sorted(d.keys()):
                part1.append(d[k] + 1)
            part1 = np.array(part1)
            num_clust1 = max(d.values())
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
    d, points = variant_distances(pdb_id='3ecr', chain='A', uniprot_id='P08397')
    links = linkage_cluster(d, methods=['single', 'complete'])
    compare_clustering(links, points)
    links = linkage_cluster(d, methods=['average', 'mcl'])
    compare_clustering(links, points)

    links = []
    for inf_fact in [2., 6.]:
        d[d > 10] = max(d)
        d = max(d) - d
        links.append([mcl(squareform(d), max_loop=50, inflate_factor=inf_fact), 'mcl_IF=' + str(inf_fact)])
    compare_clustering(links, points)

    d, points = variant_distances(pdb_id='3tnu', chain='A', uniprot_id='P02533')
    links = linkage_cluster(d, methods=['average', 'mcl'])
    compare_clustering(links, points)