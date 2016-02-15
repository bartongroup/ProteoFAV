import logging

logging.basicConfig(filename='variant_clustering.log', level=logging.DEBUG)
logger = logging.getLogger(__name__)

import sys

sys.path.extend(['/Users/smacgowan/PycharmProjects/ProteoFAV'])

from proteofav.analysis import query_uniprot
import argparse
import cPickle as pickle
import logging
from proteofav import main
import os
import time
from proteofav.utils import create_directory
from proteofav.analysis.utils import is_valid_file, create_directory

if __name__ == '__main__':

    # Parameters
    parser = argparse.ArgumentParser(description='Execute disease variant structure clustering pipeline')
    parser.add_argument('RESULTS_DIR', type=str, help='Directory to save results')
    parser.add_argument('--proteins', type=lambda x: is_valid_file(parser, x), help='File containing UniProt IDs')
    parser.add_argument('--IF', dest='inflate', type=float, help='MCL inflation factor. Will use default if ommitted')
    parser.add_argument('--threshold', dest='threshold', type=float, default=7.5,
                        help='Distance threshold to break edges of positional similarity graph')
    parser.add_argument('--similarity', dest='similarity', type=str, default='max_minus_d',
                        help='Method to transform interatomic distances into a graph')
    parser.add_argument('--retry_failed', dest='retry_failed', action='store_true')

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

    # Setup directory structure ----------------------------------------------------------------------------------------
    # Cluster analyses are stored according to the parameter sets used
    cluster_parameters_whitelist = ['inflate', 'threshold', 'similarity']
    cluster_parameters = []
    for k, v in vars(args).iteritems():
        if k in cluster_parameters_whitelist:
            cluster_parameters.append(str(k))
            cluster_parameters.append(str(v))
    cluster_dir = os.path.join('pickled', 'clusters_' + '_'.join(cluster_parameters))

    # Create required folders
    if not os.path.isdir('pickled'):
        create_directory('pickled')

    structure_table_dir = os.path.join('pickled', 'structure_tables')
    for datadir in [structure_table_dir, cluster_dir]:
        if not os.path.isdir(datadir):
            create_directory(datadir)

    # Get UniProt IDs --------------------------------------------------------------------------------------------------

    # Get suitable list of proteins
    if not args.proteins:
        protein_set = query_uniprot()[:10]  # Results from default query terms
    else:
        protein_set = [x.strip() for x in args.proteins]

    logger.info('Processing {} UniProt IDs'.format(len(protein_set)))

    # Get structure and variant data -----------------------------------------------------------------------------------

    # Get structure and uniprot variant table. Will be reloaded if previously pickled,
    # skipped if a previous attempt failed and 'retry_failed' is not set and retrieved
    # in any other case
    structure_tables = []
    for prot in protein_set:
        table_failed_placeholder = os.path.join(structure_table_dir, 'structure_table_' + prot + '.failed')
        table_pickle_name = 'structure_table_' + prot + '.pkl'
        table_file_name = os.path.join(structure_table_dir, table_pickle_name)
        if not os.path.isfile(table_file_name):
            if not os.path.isfile(table_failed_placeholder) or args.retry_failed:
                # TODO: Figure a way to complete analysis for as many proteins as possible
                logger.info('Processing UniProt ID {} out of {}...'.format(protein_set.index(prot) + 1,
                                                                           len(protein_set)))
                try:
                    structure_table = main.merge_tables(uniprot_id=prot, chain='all', uniprot_variants=True)
                    structure_tables.append((prot, structure_table))
                    with open(table_file_name, 'wb') as output:
                        pickle.dump(structure_table, output, -1)
                    if os.path.isfile(table_failed_placeholder):
                        os.remove(table_failed_placeholder)
                except:
                    logger.warning('Cannot get structure table for {}... skipping.'.format(prot))
                    with open(table_failed_placeholder, 'w') as output:
                        output.write('Last attempted at {}\n'.format(time.strftime('%a %d %b - %H:%M:%S')))
            else:
                logger.info('Structure table for {} recorded as unavailable... skipping.'.format(prot))
        else:
            structure_table = pickle.load(open(table_file_name, 'rb'))
            structure_tables.append((prot, structure_table))
            logger.info('Reloaded structure table for {}.'.format(prot))

    # Cluster analysis -------------------------------------------------------------------------------------------------
    results = []
    for prot, table in structure_tables:

        # Basic processing
        deduped = table.drop_duplicates(subset=['UniProt_dbResNum', 'chain_id'])
        mask = deduped.resn.notnull()
        n_variants = sum(mask)

        # File names
        cluster_failed_placeholder = os.path.join(cluster_dir, 'cluster_results_' + prot + '.failed')
        cluster_pickle_name = 'cluster_results_' + prot + '.pkl'
        cluster_file_name = os.path.join(cluster_dir, cluster_pickle_name)

        # Check if the results are already available or marked as failed
        if not os.path.isfile(cluster_file_name):
            if not os.path.isfile(cluster_failed_placeholder) or args.retry_failed:

                # Try to complete the cluster analysis and store the results. Otherwise, record the failure and move on.
                try:
                    logger.info('Running cluster analysis for {}.'.format(prot))
                    annotated_table = proteofav.analysis.clustering.cluster_table(deduped, mask=mask, method=['mcl_program'],
                                                                                  **vars(args))
                    cluster_table = proteofav.analysis.clustering.test_cluster_significance(annotated_table, method=['mcl_program'],
                                                                                            table=deduped, show_progress=False,
                                                                                            n_samples=50, return_samples=True,
                                                                                            **vars(args))
                    with open(cluster_file_name, 'wb') as output:
                        pickle.dump(cluster_table, output, -1)
                    results.append((prot, n_variants, cluster_table))
                except:
                    logger.warning('Cannot complete cluster analysis for {}... skipping.'.format(prot))
                    with open(cluster_failed_placeholder, 'w') as output:
                        output.write('Last attempted at {}\n'.format(time.strftime('%a %d %b - %H:%M:%S')))

            else:
                # Log and skip
                logger.info('Clustering for {} recorded as failed... skipping.'.format(prot))

        else:
            # Reload the clustering
            cluster_table = pickle.load(open(cluster_file_name, 'rb'))
            results.append((prot, n_variants, cluster_table))
            logger.info('Reloaded clustering results for {}.'.format(prot))

    # Create results plots
    names = ['mean', 'median', 'std', 'min', 'max', 'len', 'top_k_clusters', 'n_isolated', 'n_50_clusters']
    names.append('Davies-Bouldin')
    names.append('Dunn_index')
    names.append('largest_cluster_volume')
    # for prot, n_variants, stats in results:
    #     plot_failed_placeholder = os.path.join(args.RESULTS_DIR, 'plot_file_' + prot + '.failed')
    #     plot_file_name = 'cluster_metrics_' + prot + '.png'
    #     plot_file = os.path.join(args.RESULTS_DIR, plot_file_name)
    #     if not os.path.isfile(plot_file):
    #         if not os.path.isfile(plot_failed_placeholder) or args.retry_failed:
    #             try:
    #                 logger.info('Plotting results for {}.'.format(prot))
    #                 analysis.clustering.plot_sample_distributions(stats, names)
    #                 plt.suptitle(prot)
    #                 plt.tight_layout()
    #                 plt.savefig(plot_file, format='png')
    #                 plt.close()
    #             except:
    #                 logger.warning('Cannot provide plot for {}... skipping.'.format(prot))
    #                 with open(plot_failed_placeholder, 'w') as output:
    #                     output.write('Last attempted at {}\n'.format(time.strftime('%a %d %b - %H:%M:%S')))
    #     else:
    #         logger.info('Results already plotted for {}... skipping.'.format(prot))

    # Write some summary stats
    results_file = os.path.join(args.RESULTS_DIR, 'results_summary.txt')
    with open(results_file, 'w+') as summary:
        # Results Header
        names_pvalues = []
        names_psmall = []
        names_plarge = []
        for var in names:
            names_pvalues.append('p(s/e/l)_' + var)
            names_psmall.append('p_small_' + var)
            names_plarge.append('p_large_' + var)
        summary.write('# p-values indicate the proportion of random samples are smaller/equal/larger than observed.\n')
        summary.write('UniProtID\tN_variants\t' + '\t'.join(names) + '\t' + '\t'.join(names_psmall) + '\t' +
                      '\t'.join(names_plarge) + '\t' + '\t'.join(names_pvalues) + '\n')

        for i in results:
            uniprot_id = i[0]
            n_variants = str(i[1])

            p_values = []
            p_small = []
            p_large = []
            observed_stats = []
            for var in names:
                observed_stats.append(str(i[2]['obs_stats'][var]))
                p_tuple = i[2]['p'][var]
                p_values.append(str(p_tuple))

                p_num = [0 if isinstance(x, str) else x for x in p_tuple]
                p_small.append(str(sum(p_num[:2])))
                p_large.append(str(sum(p_num[1:])))


            summary.write('\t'.join([uniprot_id, n_variants] + observed_stats + p_small + p_large + p_values) + '\n')
