

import analysis.clustering
import logging
import requests
import main
import matplotlib.pyplot as plt
import os.path
import time

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

def query_uniprot(search_terms=('keyword:Disease', 'reviewed:yes', 'organism:human', 'database:(type:pdb)')):
    """
    Query the UniProt API for proteins that have particualar characteristics.

    :param search_terms: A tuple of UniProt Query search terms.
    :return: A list of UniProt IDs
    """
    url = 'http://www.uniprot.org/uniprot'
    params = {'query': ' AND '.join(search_terms),
              'format': 'tab', 'columns': 'id'}
    logging.info('Querying UniProt DB...')
    r = requests.get(url, params=params)
    uniprots = r.content.split('\n')[1:]
    logging.info('Retreived {} UniProt IDs matching query.'.format(len(uniprots)))

    return uniprots


if __name__ == '__main__':

    logging.disable(logging.DEBUG)
    logging.getLogger("requests").setLevel(logging.WARNING)

    # Get suitable list of proteins and their structure / variant data
    protein_set = query_uniprot()[:50]  # Results from default query terms
    logging.info('Processing {} UniProt IDs'.format(len(protein_set)))
    structure_tables = []
    for prot in protein_set:
        plot_file = 'cluster_figs/sample_stats/cluster_stats_' + prot + '.png'
        if not os.path.isfile(plot_file):
            # TODO: Figure a way to complete analysis for as many proteins as possible
            logging.info('Processing UniProt ID {} out of {}...'.format(protein_set.index(prot), len(protein_set)))
            try:
                structure_tables.append((prot, main.merge_tables(uniprot_id=prot, chain='all', uniprot_variants=True)))
            except:
                log.warning('Cannot get structure table for {}... skipping.'.format(prot))
        else:
            log.info('Results for {} already exist.'.format(prot))

    # Run analysis
    results = []
    for prot, table in structure_tables:
        deduped = table.drop_duplicates(subset=['UniProt_dbResNum', 'chain_id'])
        mask = deduped.resn.notnull()
        n_variants = sum(mask)
        try:
            log.info('Running cluster analysis for {}.'.format(prot))
            results.append((prot, analysis.clustering.cluster_table(deduped, mask=mask, method=['mcl_program'],
                                                                    n_samples=50, threshold=7.5,
                                                                    return_samples=True))
                           )
        except:
            log.warning('Cannot complete cluster analysis for {}... skipping.'.format(prot))

    # Create results plots

    names = ['mean', 'median', 'std', 'min', 'max', 'len', 'top_k_clusters', 'n_isolated', 'n_50_clusters']
    names.append('Davies-Bouldin')
    names.append('Dunn_index')
    names.append('largest_cluster_volume')

    for prot, stats in results:
        try:
            log.info('Plotting results for {}.'.format(prot))
            analysis.clustering.plot_sample_distributions(stats, names)
            plt.suptitle(prot)
            plt.tight_layout()
            plt.savefig(plot_file, format='png')
            plt.close()
        except:
            log.warning('Cannot provide plot for {}... skipping.'.format(prot))

