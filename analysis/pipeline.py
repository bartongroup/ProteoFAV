
import logging
logging.basicConfig(filename='variant_clustering.log', level=logging.DEBUG)
logger = logging.getLogger(__name__)

import sys
sys.path.extend(['/Users/smacgowan/PycharmProjects/ProteoFAV'])

import analysis.clustering
import argparse
import cPickle as pickle
import logging
import requests
import main
import matplotlib.pyplot as plt
import os
import time
from config import defaults

def query_uniprot(search_terms=('keyword:Disease', 'reviewed:yes', 'organism:human', 'database:(type:pdb)')):
    """
    Query the UniProt API for proteins that have particualar characteristics.

    :param search_terms: A tuple of UniProt Query search terms.
    :return: A list of UniProt IDs
    """
    url = 'http://www.uniprot.org/uniprot'
    params = {'query': ' AND '.join(search_terms),
              'format': 'tab', 'columns': 'id'}
    logger.info('Querying UniProt DB...')
    r = requests.get(url, params=params)
    uniprots = r.content.split('\n')[1:]
    logger.info('Retreived {} UniProt IDs matching query.'.format(len(uniprots)))

    return uniprots


if __name__ == '__main__':

    # Parameters
    parser = argparse.ArgumentParser(description='Execute disease variant structure clustering pipeline')
    parser.add_argument('RESULTS_DIR', type=str, help='Directory to save results')
    parser.add_argument('--IF', dest='inflate', type=float, help='MCL inflation factor. Will use default if ommitted')
    parser.add_argument('--threshold', dest='threshold', type=float, default=7.5,
                        help='Distance threshold to break edges of positional similarity graph')

    # Setup results dir
    args = parser.parse_args()
    if not os.path.isdir(args.RESULTS_DIR):
        os.makedirs(args.RESULTS_DIR)

    # Logging setup and log pipeline configuration

    # define a Handler which writes INFO messages or higher to the sys.stderr
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    logger.addHandler(console)

    logging.getLogger("requests").setLevel(logging.WARNING)
    logger.info('Starting disease variant clustering pipeline.')
    for arg, value in sorted(vars(args).items()):
        logger.info("Pipeline argument %s: %r", arg, value)


    # Get suitable list of proteins and their structure / variant data
    protein_set = query_uniprot()[:10]  # Results from default query terms
    logger.info('Processing {} UniProt IDs'.format(len(protein_set)))
    structure_tables = []
    for prot in protein_set:
        table_pickle_name = 'structure_table_' + prot + '.pkl'
        table_file_name = os.path.join(defaults.db_analysis, table_pickle_name)
        if not os.path.isfile(table_file_name):
            # TODO: Figure a way to complete analysis for as many proteins as possible
            logger.info('Processing UniProt ID {} out of {}...'.format(protein_set.index(prot), len(protein_set)))
            try:
                structure_table = main.merge_tables(uniprot_id=prot, chain='all', uniprot_variants=True)
                structure_tables.append((prot, structure_table))
                with open(table_file_name, 'wb') as output:
                    pickle.dump(structure_table, output, -1)
            except:
                logger.warning('Cannot get structure table for {}... skipping.'.format(prot))
        else:
            structure_table = pickle.load(open(table_file_name, 'rb'))
            structure_tables.append((prot, structure_table))
            logger.info('Reloaded structure table for {}.'.format(prot))

    # Run analysis
    results = []
    for prot, table in structure_tables:
        deduped = table.drop_duplicates(subset=['UniProt_dbResNum', 'chain_id'])
        mask = deduped.resn.notnull()
        n_variants = sum(mask)
        cluster_pickle_name = 'cluster_results_' + prot + '.pkl'
        cluster_file_name = os.path.join(defaults.db_analysis, cluster_pickle_name)
        if not os.path.isfile(cluster_file_name):
            try:
                logger.info('Running cluster analysis for {}.'.format(prot))
                cluster_table = analysis.clustering.cluster_table(deduped, mask=mask, method=['mcl_program'],
                                                                  n_samples=50, return_samples=True, **vars(args))
                with open(cluster_file_name, 'wb') as output:
                    pickle.dump(cluster_table, output, -1)
                results.append((prot, n_variants, cluster_table))
            except:
                logger.warning('Cannot complete cluster analysis for {}... skipping.'.format(prot))
        else:
            cluster_table = pickle.load(open(cluster_file_name, 'rb'))
            results.append((prot, n_variants, cluster_table))
            logger.info('Reloaded clustering results for {}.'.format(prot))

    # Create results plots

    names = ['mean', 'median', 'std', 'min', 'max', 'len', 'top_k_clusters', 'n_isolated', 'n_50_clusters']
    names.append('Davies-Bouldin')
    names.append('Dunn_index')
    names.append('largest_cluster_volume')

    for prot, n_variants, stats in results:
        plot_file_name = 'cluster_metrics_' + prot + '.png'
        plot_file = os.path.join(args.RESULTS_DIR, plot_file_name)
        if not os.path.isfile(plot_file):
            try:
                logger.info('Plotting results for {}.'.format(prot))
                analysis.clustering.plot_sample_distributions(stats, names)
                plt.suptitle(prot)
                plt.tight_layout()
                plt.savefig(plot_file, format='png')
                plt.close()
            except:
                logger.warning('Cannot provide plot for {}... skipping.'.format(prot))
        else:
            logger.info('Results already plotted for {}... skipping.'.format(prot))

    # Write some summary stats
    results_file = os.path.join(args.RESULTS_DIR, 'results_summary.txt')
    with open(results_file, 'w+') as summary:
        summary.write('# p-values indicate the proportion of random samples are smaller/equal/larger than observed.\n')
        summary.write('UniProtID\tN_variants\tmax_cluster_size\tp(small/equal/large)\n')
        for i in results:
            summary.write('\t'.join([i[0], str(i[1]), str(i[2]['obs_stats']['max']),
                                     str(i[2]['p']['max']) + '\n']))

