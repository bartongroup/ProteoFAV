

import analysis.clustering
import requests
import main


def query_uniprot(search_terms=('keyword:Disease', 'reviewed:yes', 'organism:human', 'database:(type:pdb)')):
    """
    Query the UniProt API for proteins that have particualar characteristics.

    :param search_terms: A tuple of UniProt Query search terms.
    :return: A list of UniProt IDs
    """
    url = 'http://www.uniprot.org/uniprot'
    params = {'query': ' AND '.join(search_terms),
              'format': 'tab', 'columns': 'id'}
    r = requests.get(url, params=params)
    uniprots = r.content.split('\n')[1:]

    return uniprots


if __name__ == '__main__':

    # Get suitable list of proteins and their structure / variant data
    protein_set = query_uniprot()[:50]  # Results from default query terms
    structure_tables = []
    for prot in protein_set:
        structure_tables.append((prot, main.merge_tables(uniprot_id=prot, chain='all', uniprot_variants=True)))

    # Run analysis
    results = []
    for prot, table in structure_tables:
        deduped = table.drop_duplicates(subset=['UniProt_dbResNum', 'chain_id'])
        mask = deduped.resn.notnull()
        n_variants = sum(mask)
        if n_variants >= 10:  # For now, only look at structures with this number of variants
            results.append((prot, analysis.clustering.cluster_table(deduped, mask=mask, method=['mcl_program'],
                                                                    n_samples=50, threshold=7.5,
                                                                    return_samples=True))
                           )
