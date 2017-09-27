#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import absolute_import

import logging

import pandas as pd

from proteofav.fetchers import fetch_uniprot_sequence
from proteofav.fetchers import fetch_uniprot_formal_specie
from proteofav.fetchers import _uniprot_info
from proteofav.fetchers import map_gff_features_to_sequence
from proteofav.fetchers import _uniprot_to_ensembl_xref
from proteofav.fetchers import _fetch_ensembl_variants
from proteofav.fetchers import _fetch_sequence_from_ensembl_protein
from proteofav.fetchers import _fetch_uniprot_ensembl_mapping
from proteofav.library import valid_ensembl_species

__all__ = ["_match_uniprot_ensembl_seq",
           "_apply_sequence_index_map",
           "_compare_sequences",
           "_count_mismatches",
           "select_variants",
           "parse_uniprot_variants",
           "select_uniprot_variants"]
log = logging.getLogger('proteofav.config')


##############################################################################
# Private methods
##############################################################################
def _match_uniprot_ensembl_seq(uniprot_id):
    """Matches a Uniprot sequence to a identical Ensembl sequence making the
    residue mapping trivial.

    :param uniprot_id:
    :return:
    :raises ValueError: when no match is found between uniprot and ensembl

        ..info:
      ensembl of UNIPROT ACCESSION with more than one gene P04637
    """

    species = fetch_uniprot_formal_specie(uniprot_id)
    species = species.lower().replace(' ', '_')
    if species not in valid_ensembl_species:
        raise ValueError('{} is not a valid Ensembl specie ({}).'.format(
            species, uniprot_id))
    # we remove the isoform identifier to xref to ensembl, but keep it while checking the seq
    ensembl_xref = _uniprot_to_ensembl_xref(uniprot_id.split('-')[0],
                                            species=species)

    if ensembl_xref.empty:
        raise ValueError('No cross-reference found for {} Ensembl mapping.'.format(
            uniprot_id))

    ensembl_ptn_ids = ensembl_xref.loc[ensembl_xref.type == 'translation', 'id']

    uniprot_sequence = fetch_uniprot_sequence(uniprot_id)
    for ensembl_ptn_id in ensembl_ptn_ids:
        ensembl_ptn_seq = _fetch_sequence_from_ensembl_protein(ensembl_ptn_id)
        if _compare_sequences(uniprot_sequence, ensembl_ptn_seq, permissive=False):
            return ensembl_ptn_id
    raise ValueError('No protein with the same sequence was retrieved from Ensembl {}'.format(
        uniprot_id))


def _apply_sequence_index_map(indexes, imap):  # TODO revise
    """

    :param indexes:
    :param map:
    :return:
    """

    # perform the raw translation
    translation = []
    for i in indexes:
        equivalent = imap.get(i)
        translation.append(equivalent)

    return translation


def _compare_sequences(sequence1, sequence2, permissive=True, n_mismatches=0):
    """Compares two given sequences in terms of length and sequence content.

    :param sequence1: First sequence
    :param sequence2: Second sequence
    :param permissive: if True it allow sequences with different sizes to return True
    :param n_mismatches: number of allowed mismatches
    :return: simply a true or false
    :rtype: boolean
    """
    if not permissive and len(sequence1) != len(sequence2):
        return False

    if _count_mismatches(sequence1, sequence2) > n_mismatches:
        return False

    return True


def _count_mismatches(sequence1, sequence2):
    """
    Counts the number of mismatches between two sequences
    of the same length.

    :param sequence1: sequence 1
    :param sequence2: sequence 2
    :return: The number of mismatches between sequences 1 and 2.
    """
    return sum(i != j for i, j in zip(sequence1, sequence2))


##############################################################################
# Public methods
##############################################################################
def select_variants(uniprot_id, features=('uniprot', 'ensembl_somatic', 'ensembl_germline')):
    """

    :param features:
    :param uniprot_id:
    """
    ensembl_ptn_id = None
    tables = []
    # use the uniprot natural variants as our reference (left) table
    if 'uniprot' in features:
        table_uni = parse_uniprot_variants(uniprot_id)
        if not table_uni.empty:
            table_uni = table_uni['ids']
            table_uni.name = 'uniprot_variants'
        else:
            table_uni['uniprot_variants'] = None

        tables.append(table_uni)

    try:
        ensembl_ptn_id = _match_uniprot_ensembl_seq(uniprot_id)

    except ValueError as e:
        # No Ensembl mapping
        log.error(e)
        return pd.concat(tables, axis=1)

    if 'ensembl_somatic' in features:
        table_som = _fetch_ensembl_variants(
            ensembl_ptn_id, feature='somatic_transcript_variation')

        if not table_som.empty:
            table_som = table_som[table_som['type'] == 'coding_sequence_variant']
            table_som = table_som.groupby('start')['id'].apply(list)
            table_som.name = 'somatic_variants'
        else:
            table_som['somatic_variants'] = None

        tables.append(table_som)

    if 'ensembl_germline' in features:
        table_ger = _fetch_ensembl_variants(ensembl_ptn_id, feature='transcript_variation')

        if not table_ger.empty:
            table_ger = table_ger[table_ger['type'] == 'missense_variant']
            table_ger = table_ger.groupby('start')['id'].apply(list)
            table_ger.name = 'germline_variants'
        else:
            table_ger['germline_variants'] = None

        tables.append(table_ger)

    return pd.concat(tables, axis=1)


def parse_uniprot_variants(uniprot_id):
    """

    :param uniprot_id: Uniprot accession to a protein P04637
    :return: table[disease: list of str,
                   transition: list of str (residue pairs),
                   ids: list of str]
    :rtype: pandas.DataFrame
    """

    disease_group = '\[\'In ([?P<disease>a-zA-Z0-9_ ]+)[.;]'
    res_transition_group = '(?P<ref>[A-Z]+)->(?P<new>[A-Z]+)'
    ids_group = "\(\[\'([a-zA-Z0-9_]+)\'\]\)"

    table = map_gff_features_to_sequence(uniprot_id, query_type='Natural variant')

    if table.empty:
        return table

    table['disease'] = table['annotation'].str.findall(disease_group)
    table['transition'] = table['annotation'].str.findall(res_transition_group)
    table['ids'] = table['annotation'].str.findall(ids_group)

    del table['annotation']

    return table


# TODO we need to figure out a manual function to map sequences that have
# small differences
def select_uniprot_variants(identifier, align_transcripts=False):
    """
    Summarise variants for a protein in the UniProt

    :type align_transcripts: object
    :param identifier: UniProt ID
    :return: table with variants, rows are residues
    """
    # TODO FIX docstring
    # get organism and sequence for the provided identifier

    # TODO use gff?
    uni = _uniprot_info(identifier, cols=['organism', 'sequence'])
    org = ('_'.join(uni.loc[0, 'Organism'].split()[-3:-1])).lower()
    seq = uni.loc[0, 'Sequence']

    # get the ensembl ids: this also validate this species as available
    # through ensembl
    ens = _fetch_uniprot_ensembl_mapping(identifier, species=org)

    # get the ensembl protein ids
    ens_pros = ens.loc[0, 'TRANSLATION']
    if not isinstance(ens_pros, list):
        ens_pros = [ens_pros, ]

    # get the sequence of the ensembl protein
    usable_indexes = []
    aligned_indexes = []
    seq_maps = []

    for i, enspro in enumerate(ens_pros):
        seq_pro = _fetch_sequence_from_ensembl_protein(enspro, protein=True)

        # validate if the sequence of uniprot and ensembl protein matches
        if _compare_sequences(seq, seq_pro, permissive=False):
            usable_indexes.append(i)
        elif _compare_sequences(seq, seq_pro, permissive=True):
            usable_indexes.append(i)
            seq_maps.append(None)
            n_mismatches = _count_mismatches(seq, seq_pro)
            message = "{0}: Sequences are of same length but have {1} mismatch(s)".format(
                enspro, n_mismatches)
            logging.warning(message)
        elif align_transcripts:
            message = "Sequences don't match! Will attempt alignment... {}".format(enspro)
            logging.warning(message)
            aligned_indexes.append(i)
            ensembl_to_uniprot = _apply_sequence_index_map(seq_pro, seq)
            seq_maps.append(ensembl_to_uniprot)
        else:
            message = "Sequences don't match! skipping... {}".format(enspro)
            logging.warning(message)
            seq_maps.append(None)

    # get the variants for the ensembl proteins that match the uniprot
    tables = []
    for i in usable_indexes:
        vars = _fetch_ensembl_variants(ens_pros[i], feature='transcript_variation')
        muts = _fetch_ensembl_variants(ens_pros[i], feature='somatic_transcript_variation')

        # TODO: From... TO... mutated residues as different columns in the table
        # TODO: Shouldn't the default behaviour return all the columns?
        if not vars.empty:
            tables.append(vars[['translation', 'id', 'start', 'residues']])
        if not muts.empty:
            tables.append(muts[['translation', 'id', 'start', 'residues']])

    # Get variants from aligned sequences
    if align_transcripts:
        for i in aligned_indexes:
            vars = _fetch_ensembl_variants(ens_pros[i], feature='transcript_variation')
            muts = _fetch_ensembl_variants(ens_pros[i], feature='somatic_transcript_variation')

            var_table = vars[['translation', 'id', 'start', 'residues']]
            mut_table = muts[['translation', 'id', 'start', 'residues']]

            var_table.start = _apply_sequence_index_map(var_table.start, seq_maps[i])
            mut_table.start = _apply_sequence_index_map(mut_table.start, seq_maps[i])

            if not var_table.empty:
                tables.append(var_table)
            if not var_table.empty:
                tables.append(mut_table)

    # to_unique = lambda series: series.unique()
    # return table.groupby('start').apply(to_unique)
    table = pd.concat(tables)
    return table


if __name__ == '__main__':
    pass
