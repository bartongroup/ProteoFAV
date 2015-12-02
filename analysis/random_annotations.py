__author__ = 'smacgowan'

"""
Created on 25/11/2015
Functions to help with random selection of PDB atoms/residues and significance testing via bootstrapping.
"""

from numpy.random import choice, random_integers
from numpy import array
from pandas import Series, DataFrame
from string import letters
from random import shuffle


def random_uniprot_patho_table(merge_table, n_residues, n_phenotypes=1):

    # TODO: if I read a variant table, I can get the actual AA composition and dssp and keep these constant

    # We are going to pick residues from the supplied merge_table, so first get rid of unobserved residues
    merge_table = merge_table[merge_table.aa.notnull()]

    # Pick out the columns we want
    columns = []
    for col in ['UniProt_dbResNum', 'aa']:
        columns.append(merge_table.columns.get_loc(col))

    # Create a random selection
    nrow = len(merge_table)
    rows = random_integers(0, nrow - 1, n_residues)
    table = merge_table.iloc[rows, columns]

    # Now generate random variant residues and phenotypes
    selection = random_integers(0, n_phenotypes - 1, n_residues)
    disease = [letters[i] for i in selection]
    aa1 = array(list('ACDEFGHIKLMNPQRSTVWY'))
    mut = choice(aa1, n_residues, replace=True)  # TODO: set p to match a variant distribution

    # And add them to table and make it exactly like a real uniprot_variant table
    table['mut'] = mut
    table['disease'] = disease
    table = table.reset_index()
    del table['PDB_dbResNum']
    table = table.rename(columns={'aa': 'resn'})

    # TODO: Fix any entries where resn == mut
    # eq = table.resn == table.mut
    # vars.mut[eq] = 'K'  ## has warning but works

    return table




def add_random_disease_variants(merge_table, n_residues, n_phenotypes):
    """
    Add phenotype annotations to random entries on a PDB merged table.
    :param merge_table:
    :param n_residues:
    :param n_phenotypes:
    :return:
    """
    # If we have more than one chain then we need to use split/apply/combine
    n_chains = len(merge_table.chain_id.dropna().unique())
    if n_chains == 1:
        variants = random_uniprot_patho_table(merge_table, n_residues, n_phenotypes)
        table = merge_table.merge(variants, how='left')
        return table

    else:
        # Randomly distribute the number of variants to each chain
        n_vars_per_chain = []
        remaining = n_residues
        for _ in xrange(n_chains - 1):
            n_vars_per_chain.append(random_integers(0, remaining, 1))
            remaining = remaining - n_vars_per_chain[-1]
        n_vars_per_chain.append(remaining)
        shuffle(n_vars_per_chain)

        # Now split, apply, combine by iteration so we can vary the number of residues
        grouped = merge_table.groupby('chain_id')  # TODO: This drops any where chain_id is NaN, fix or log
        newtable = DataFrame()
        for i, (name, group) in enumerate(grouped):
            newtable = newtable.append(add_random_disease_variants(group, n_vars_per_chain[i], n_phenotypes))
        return newtable



