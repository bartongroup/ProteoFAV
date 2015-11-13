__author__ = 'smacgowan'

# Starting out with the example from http://stackoverflow.com/questions/21638130/tutorial-for-scipy-cluster-hierarchy

import numpy as np
import scipy.cluster.hierarchy as hac
import matplotlib.pyplot as plt

from main import merge_tables
from to_table import _fetch_uniprot_variants
import pandas as pd
from scipy.spatial.distance import pdist
from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial import ConvexHull, Delaunay
from scipy.spatial.qhull import QhullError

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


def linkage_cluster(a, xyz):
    fig, axes23 = plt.subplots(2, 3)

    for method, axes, offset in zip(['single', 'complete'], axes23, [0, 3]):
        z = hac.linkage(a, method=method)

        # Plotting
        axes[0].plot(range(1, len(z)+1), z[::-1, 2])
        knee = np.diff(z[::-1, 2], 2)
        axes[0].plot(range(2, len(z)), knee)

        num_clust1 = knee.argmax() + 2
        knee[knee.argmax()] = 0
        num_clust2 = knee.argmax() + 2

        axes[0].text(num_clust1, z[::-1, 2][num_clust1-1], 'possible\n<- knee point')

        if num_clust1 > 12:
            num_clust1 = 12
        if num_clust2 > 12:
            num_clust2 = 12

        part1 = hac.fcluster(z, num_clust1, 'maxclust')
        part2 = hac.fcluster(z, num_clust2, 'maxclust')

        clr = ['#2200CC' ,'#D9007E' ,'#FF6600' ,'#FFCC00' ,'#ACE600' ,'#0099CC' ,
        '#8900CC' ,'#FF0000' ,'#FF9900' ,'#FFFF00' ,'#00CC01' ,'#0055CC']

        for part, i in zip([part1, part2], [2, 3]):
            ax = fig.add_subplot(2, 3, i + offset, projection='3d')
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
    plt.show()


if __name__ == '__main__':
    d, points = variant_distances(pdb_id='3ecr', chain='A', uniprot_id='P08397')
    linkage_cluster(d, points)