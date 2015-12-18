#!/usr/bin/env python
# -*- coding: utf-8

"""
Created on 25/11/2015
Functions to help with random selection of PDB atoms/residues and significance
testing via bootstrapping.
"""

import logging
from numpy.random import choice, random_integers, permutation
from numpy import array
from pandas import Series, DataFrame, crosstab
from string import letters
from random import shuffle

__author__ = 'smacgowan'


def random_uniprot_patho_table(merge_table, n_residues, n_phenotypes=1,
                               proportions=None, repeats_allowed=True, force=False):
    """

    :param merge_table:
    :param n_residues:
    :param n_phenotypes:
    :param proportions: A Pandas Series denoting the number of each residue classification required
    as produced by Series.value_counts on an amino acid column OR a Pandas DataFrame containing
    a cross tabulation of two attributes giving the proportion of each residue class required as
    produced by `pandas.crosstab`.
    :return:
    """

    # TODO: if I read a variant table, I can get the actual AA composition and dssp and
    # keep these constant
    # TODO: currently different chains are handled by the `add_random_disease_variants`
    # wrapper, is this ideal?

    # We are going to pick residues from the supplied merge_table, so first get rid of
    # unobserved residues
    merge_table = merge_table[merge_table.aa.notnull()].reset_index()

    # Create a random selection
    ## TODO: Were rand indexes better as allowed repeats
    shuffled_table = merge_table.reindex(permutation(merge_table.index))

    if repeats_allowed:
        # Copy all rows 21 times so that we have p > 0 of all possible variants on the same residue
        row_replicator = range(len(shuffled_table)) * 21
        shuffled_table = shuffled_table.iloc[row_replicator, :].reset_index()
        shuffled_table = shuffled_table.reindex(permutation(shuffled_table.index))

    # By now index has been reset as neccessary so can pick columns
    columns = []
    for col in ['UniProt_dbResNum', 'aa']:
        columns.append(shuffled_table.columns.get_loc(col))

    if proportions is None:
        table = shuffled_table.iloc[:n_residues, columns]

    elif type(proportions) == Series:
        table = DataFrame()
        n_residues = sum(proportions)  ## TODO: Fix this!
        column_name = proportions.name
        for value, count in proportions.iteritems():
            mask = shuffled_table[[column_name]] == value
            rows = range(count)
            try:
                table = table.append(shuffled_table[mask].iloc[rows, columns])
            except IndexError:
                msg = 'Fewer residues with {}: {} available than requested including repeats!'
                logging.error(msg.format(column_name, value))
                if force:
                    msg = 'Falling back to giving all variants on this residue type.'
                    logging.warning(msg)
                    table = table.append(shuffled_table[mask].iloc[:, columns])
                    n_residues = n_residues - count + sum(mask)
                else:
                    msg = 'Could not generate variants with requested proportions. Try changing ' \
                          'proportions or use force=True.'
                    logging.error(msg)
                    raise


            # If number drawn forces replicates, report it
            n_choices = sum(merge_table[column_name] == value)
            ratio = count / float(n_choices)
            if count > n_choices:
                msg = 'FORCED REPEATS: Fewer residues with {}: {} available than requested.'
                logging.warning(msg.format(column_name, value))
            elif ratio > 0.8:
                msg = 'RESTRICTED SELECTION: {}% of available residues with {}: {} requested.'
                logging.warning(msg.format(round(ratio * 100), column_name, value))

    elif type(proportions) == DataFrame:
        # TODO: need to handle NaNs in `ss`
        # In this case it is a cross tab and we want the proportions fixed for each pair
        table = DataFrame()
        n_residues = proportions.sum().sum()  ## TODO: Fix this!
        column_name = proportions.columns.name
        row_name = proportions.index.name
        for column, series in proportions.iteritems():
            for row, count in series.iteritems():
                mask = (shuffled_table[column_name] == column) & (shuffled_table[row_name] == row)
                rows = range(count)
                try:
                    table = table.append(shuffled_table[mask].iloc[rows, columns])
                except IndexError:
                    msg = 'Fewer residues with {}: {} and {}: {} available than requested ' \
                          'including repeats!'
                    logging.error(msg.format(column_name, column, row_name, row))
                    if force:
                        msg = 'Falling back to giving all variants on this residue type.'
                        logging.warning(msg)
                        table = table.append(shuffled_table[mask].iloc[:, columns])
                        n_residues = n_residues - count + sum(mask)
                    else:
                        msg = 'Could not generate variants with requested proportions. ' \
                              'Try changing proportions or use force=True.'
                        logging.error(msg)
                        raise

                # If number drawn forces replicates, report it
                orig_mask = (merge_table[column_name] == column) & (merge_table[row_name] == row)
                n_choices = sum(orig_mask)
                ratio = count / float(n_choices)
                if count > n_choices:
                    msg = 'FORCED REPEATS: Fewer residues with {}: {} and {}: {} available ' \
                          'than requested.'
                    logging.warning(msg.format(column_name, column, row_name, row))
                elif ratio > 0.8:
                    msg = 'RESTRICTED SELECTION: {}% of available residues with {}: {} ' \
                          'and {}: {} requested.'
                    logging.warning(msg.format(round(ratio * 100), column_name, column, row_name, row))

    # Now generate random variant residues and phenotypes
    selection = random_integers(0, n_phenotypes - 1, n_residues)
    disease = [letters[i] for i in selection]
    aa1 = array(list('ACDEFGHIKLMNPQRSTVWY'))
    # TODO: set p to match a variant distribution
    mut = choice(aa1, n_residues, replace=True)

    # And add them to table and make it exactly like a real uniprot_variant table
    table['mut'] = mut
    table['disease'] = disease
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
    n_prots = len(merge_table.UniProt_dbAccessionId.dropna().unique())
    n_chains = len(merge_table.chain_id.dropna().unique())
    if n_prots == 1 and n_chains == 1:
        variants = random_uniprot_patho_table(merge_table, n_residues, n_phenotypes)
        table = merge_table.merge(variants, how='left')
        return table

    elif n_chains == 1:
        # This means we have one chain that maps to multiple UniProts so we can
        # randomly distribute the number of variants to each unique protein
        n_vars_per_prot = []
        remaining = n_residues
        for _ in xrange(n_prots - 1):
            n_vars_per_prot.append(random_integers(0, remaining, 1))
            remaining = remaining - n_vars_per_prot[-1]
        n_vars_per_prot.append(remaining)
        shuffle(n_vars_per_prot)

        # Now split, apply, combine by iteration so we can vary the number of residues
        # TODO: This drops any where uniprot is NaN, fix or log
        grouped = merge_table.groupby('UniProt_dbAccessionId')
        table = DataFrame()
        for i, (name, group) in enumerate(grouped):
            table = table.append(add_random_disease_variants(group, n_vars_per_prot[i],
                                                             n_phenotypes))
        return table

    elif n_prots == 1:
        # This means there is one protein represented by multiple chains. Depending
        # on the degree of overlaping UniProt coverage we have to avoid doubling the
        # total number of variants.

        mask = merge_table.chain_id.notnull()
        n_dupe = sum(merge_table[mask].duplicated(subset='UniProt_dbResNum', keep=False)) / n_chains
        pc_duplicated = float(n_dupe) / len(merge_table[mask].UniProt_dbResNum.dropna().unique())
        adjustment = (1. / n_chains) * pc_duplicated
        adj_n_residues = int(round(n_residues * adjustment))

        # Now fetch variants and merge as usual
        discrepancy = 2
        while abs(discrepancy) > 1:
            variants = random_uniprot_patho_table(merge_table, adj_n_residues, n_phenotypes)
            table = merge_table.merge(variants, how='left')
            discrepancy = sum(table.resn.notnull()) - n_residues
            adj_n_residues = adj_n_residues - (discrepancy / 2)

        return table

    else:
        grouped = merge_table.groupby(['UniProt_dbAccessionId', 'chain_id'])

        # Work out the divisor for even distribution
        cross = crosstab(merge_table.UniProt_dbAccessionId, merge_table.chain_id)
        num = cross._get_numeric_data()
        num[num > 0] = 1
        divisor = cross.sum().sum()

        n_residues = int(round(n_residues / divisor))
        table = DataFrame()
        for name, group in grouped:
            table = table.append(add_random_disease_variants(group, n_residues, n_phenotypes))
        return table
