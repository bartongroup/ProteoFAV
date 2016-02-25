#!/usr/bin/env python
# -*- coding: utf-8 -*-


import cPickle as pickle
import logging
import os

import pandas as pd
import requests
from Bio import pairwise2

from proteofav.config import defaults
from proteofav.uniprot import (map_gff_features_to_sequence, _uniprot_to_ensembl_xref,
                               _uniprot_info, get_uniprot_formal_specie, get_uniprot_sequence)
from proteofav.utils import (get_url_or_retry, check_local_or_fetch)

log = logging.getLogger(__name__)


##############################################################################
# Private methods
##############################################################################
def _fetch_ensembl_variants(ensembl_ptn_id, feature=None):
    """Queries the Ensembl API for germline variants (mostly dbSNP) and somatic
    (mostly COSMIC) based on Ensembl Protein identifiers (e.g. ENSP00000326864).

    :param ensembl_ptn_id: Ensembl accession to a protein: ENSP00000XXXXXX
    :return: table[Parent: str,
                   allele§: str,
                   clinical_significance: list,
                   codons: str,
                   end: int,
                   feature_type: str,
                   id: str,
                   minor_allele_frequency: float ,
                   polyphen: float,
                   residues: str,
                   seq_region_name: str,
                   sift: float,
                   start: int,
                   translation: str,
                   type: str ]
    :rtype: pandas.DataFrame
    """
    ensembl_endpoint = "overlap/translation/"
    suported_feats = ['transcript_variation', 'somatic_transcript_variation']
    if feature is None:
        raise NotImplementedError('Use two functions call to get both somatic'
                                  ' and germline variants.')
        # params = {'feature': supported_feats,
        #           'type': 'missense_variant'}
    elif feature not in suported_feats:
        raise NotImplementedError(
                'feature argument should be one of {} or None for all'.format(
                        ', '''.join(suported_feats)))
    else:
        params = {'feature': feature}
    url = defaults.api_ensembl + ensembl_endpoint + ensembl_ptn_id

    rows = get_url_or_retry(url, json=True, **params)
    return pd.DataFrame(rows)


def _fetch_uniprot_variants(identifier, _format='tab'):
    """
    Request human curated variants from UniProt.

    :param identifier: UniProt ID
    :param _format: request output from the webserver
    :return: pandas dataframe
    """
    # Check if query has already been saved

    query_file_name = defaults.db_germline_variants + 'uniprot_variants_' + identifier + '.pkl'
    if not os.path.isfile(query_file_name):
        url = defaults.api_uniprot + '?query=accession:' + identifier
        url += '&format=' + _format
        url += '&columns=feature(NATURAL+VARIANT)'
        result = get_url_or_retry(url)
        with open(query_file_name, 'wb') as output:
            pickle.dump(result, output, -1)
    else:
        # Read stored response
        log.debug('Reloading stored response.')
        result = pickle.load(open(query_file_name, 'rb'))

    # Complicated parsing
    records = result.replace('Natural variant\n', '').split('.; ')
    variants = [['UniProt_dbResNum', 'resn', 'mut', 'disease']]
    for record in records:

        # Position and mutation
        entry = record.split(' ')
        resi = entry[1]
        resn = entry[3]
        if resn == 'Missing':
            mut = 'NA'
        else:
            mut = entry[5]

        # Label disease if specified
        try:
            ind = entry.index('(in')
            if entry[ind + 1].startswith('dbSNP'):
                disease = "Not present"
            else:
                disease = entry[ind + 1].replace(';', '').replace(').', '')
        except ValueError:
            disease = "Not found"

        variants.append([resi, resn, mut, disease])

    table = pd.DataFrame(variants, columns=variants.pop(0))
    table.UniProt_dbResNum = table.UniProt_dbResNum.astype('float')

    return table


def _fetch_variant_characteristics_from_identifier(identifier, species='human'):
    """
    Queries the Ensembl API for variant IDs (e.g rs376845802 or COSM302853).

    :param identifier: variant ID
    :param species: Ensembl species
    :return: pandas table dataframe
    """
    ensembl_endpoint = "variation/{}/".format(species)
    # params = {'pops': '1', 'phenotypes': '1', 'genotypes': '1'}
    url = defaults.api_ensembl + ensembl_endpoint + str(identifier)
    data = get_url_or_retry(url, json=True)

    rows = []
    information = {}
    for parent in data:
        if parent == "mappings":
            for entry in data[parent]:
                for key in entry:
                    try:
                        if entry[key] in information[key]:
                            continue
                        information[key].append(entry[key])
                    except KeyError:
                        information[key] = entry[key]
                    except AttributeError:
                        information[key] = [information[key]]
                        information[key].append(entry[key])
        else:
            information[parent] = data[parent]

    rows.append(information)
    return pd.DataFrame(rows)


def _fetch_variant_characteristics_from_identifiers(variant_ids, use_vep=False):
    """
    Retrieves variant annotation from ENSEMBL.

    :param variant_ids: Ensembl Variant identifier
    :param use_vep: whether to use predicted variants from VEP
    :return:
    """

    # POST if given a list of ids
    if isinstance(variant_ids, (list, pd.Series)):
        # Remove any nans from the list
        variant_ids = [i for i in variant_ids if not str(i) == 'nan']

        ensembl_endpoint = "variation/homo_sapiens"
        if use_vep:
            ensembl_endpoint = "vep/human/id"
        url = defaults.api_ensembl + ensembl_endpoint
        headers = {"Content-Type": "application/json",
                   "Accept": "application/json"}
        data = '{ "ids" : ' + str(variant_ids).replace("u'", "'") \
               + ', "phenotypes" : 1 }'  # FIXME
        data = data.replace("'", "\"")
        result = requests.post(url, headers=headers, data=data)
        return result
    # GET if given single id
    if isinstance(variant_ids, str):
        ensembl_endpoint = "variation/homo_sapiens/" + variant_ids
        if use_vep:
            ensembl_endpoint = "vep/human/id/" + variant_ids
        headers = {"Content-Type": "application/json"}
        params = {"phenotypes": 1}
        url = defaults.api_ensembl + ensembl_endpoint
        # TODO should use
        result = get_url_or_retry(url, json=True, header=headers, **params)
        return result


def _sequence_from_ensembl_protein(identifier, protein=True):
    """
    Gets the sequence for an Ensembl identifier.

    :param identifier: Ensembl ID
    :return: sequence
    """
    ensembl_endpoint = "sequence/id/"
    url = defaults.api_ensembl + ensembl_endpoint + str(identifier)
    header = {'content-type': 'text/plain'}
    if protein:
        params = {'type': 'protein'}
    else:
        params = {}
    sequence = get_url_or_retry(url, header=header, **params)
    return sequence


def _uniprot_ensembl_mapping(identifier, species='human'):
    """
    Uses the UniProt mapping service to try and get Ensembl IDs for the
    UniProt accession identifier provided

    :param identifier: UniProt accession identifier
    :param species: Ensembl species
    :return: pandas table dataframe
    """
    from library import valid_ensembl_species
    if species not in valid_ensembl_species:
        raise ValueError('Provided species {} is not valid'.format(species))

    information = {}
    rows = []

    ensembl_endpoint = "xrefs/symbol/{}/".format(species)
    url = defaults.api_ensembl + ensembl_endpoint + str(identifier)
    data = get_url_or_retry(url, json=True)
    for entry in data:
        typ = entry['type'].upper()
        eid = entry['id']
        try:
            if eid in information[typ]:
                continue
            information[typ].append(eid)
        except KeyError:
            information[typ] = eid
        except AttributeError:
            information[typ] = [information[typ]]
            information[typ].append(eid)

    rows.append(information)
    return pd.DataFrame(rows)


def _match_uniprot_ensembl_seq(uniprot_id):
    """Matches a Uniprot sequence to a identical Ensembl sequence making the
    residue mapping trivial.

    :param uniprot_id:
    :return:
    :raises ValueError: when no match is found between uniprot and ensembl

        ..info:
      ensembl of UNIPROT ACCESSION with more than one gene P04637
    """

    species = get_uniprot_formal_specie(uniprot_id)
    ensembl_xref = _uniprot_to_ensembl_xref(uniprot_id, species=species.lower().replace(' ', '_'))

    if ensembl_xref.empty:
        raise ValueError('Not reference to {} was found in Ensembl not found.'.format(uniprot_id))

    ensembl_ptn_ids = ensembl_xref.loc[ensembl_xref.type == 'translation', 'id']

    uniprot_sequence = get_uniprot_sequence(uniprot_id)
    for ensembl_ptn_id in ensembl_ptn_ids:
        ensembl_ptn_seq = _sequence_from_ensembl_protein(ensembl_ptn_id)
        if _compare_sequences(uniprot_sequence, ensembl_ptn_seq, permissive=False):
            return ensembl_ptn_id


# TODO improve this: length / matching sequence
# TODO: should this be coming from GFF?


def _map_sequence_indexes(from_seq, to_seq):
    """
    Gets a map between sequences.

    :param from_seq: input sequence
    :param to_seq: input sequence
    :return: a map between sequences.
    :rtype: dict
    """

    def aligned_seq_indexes(seq):
        seq_indexes = []
        i = 0
        for res in seq:
            if res != '-':
                seq_indexes.append(i)
                i += 1
            else:
                seq_indexes.append('-')
        return seq_indexes

    # build the local alignment
    alignments = pairwise2.align.localxx(from_seq, to_seq)
    scores = zip(*alignments)[2]
    message = "Alignment score(s): {}".format(scores)
    logging.info(message)
    message = "First alignment:\n" + pairwise2.format_alignment(*alignments[0])
    logging.debug(message)
    if len(scores) > 1:
        message = "Found multiple alignments, arbitrarily proceeding with the first."
        logging.warning(message)

    # create the index mapping
    seq_one = aligned_seq_indexes(alignments[0][0])
    seq_two = aligned_seq_indexes(alignments[0][1])
    outmap = dict(zip(seq_one, seq_two))

    return outmap


def _apply_sequence_index_map(indexes, imap):
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
    # let's use the uniprot natural variants as our reference (left) table
    if 'uniprot' in features:
        table_uni = parse_uniprot_variants(uniprot_id)
        table_uni = table_uni['ids']
        table_uni.name = 'uniprot_variants'

        tables.append(table_uni)

    if 'ensembl_somatic' in features:
        ensembl_ptn_id = _match_uniprot_ensembl_seq(uniprot_id)
        table_som = _fetch_ensembl_variants(ensembl_ptn_id, feature='somatic_transcript_variation')
        table_som = table_som[table_som['type'] == 'coding_sequence_variant']
        table_som = table_som.groupby('start')['id'].apply(list)
        table_som.name = 'somatic_variants'

        tables.append(table_som)

    if 'ensembl_germline' in features:
        if ensembl_ptn_id is None:
            ensembl_ptn_id = _match_uniprot_ensembl_seq(uniprot_id)
        table_ger = _fetch_ensembl_variants(ensembl_ptn_id, feature='transcript_variation')
        table_ger = table_ger[table_ger['type'] == 'missense_variant']
        table_ger = table_ger.groupby('start')['id'].apply(list)
        table_ger.name = 'germline_variants'

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
    check_local_or_fetch(uniprot_id)

    disease_group = '\[\'In ([?P<disease>a-zA-Z0-9_ ]+)[.;]'
    res_transition_group = '(?P<ref>[A-Z]+)->(?P<new>[A-Z]+)'
    ids_group = "\(\[\'([a-zA-Z0-9_]+)\'\]\)"

    table = map_gff_features_to_sequence(uniprot_id, query_type='Natural variant')
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
    ens = _uniprot_ensembl_mapping(identifier, species=org)

    # get the ensembl protein ids
    ens_pros = ens.loc[0, 'TRANSLATION']
    if not isinstance(ens_pros, list):
        ens_pros = [ens_pros, ]

    # get the sequence of the ensembl protein
    usable_indexes = []
    aligned_indexes = []
    seq_maps = []

    for i, enspro in enumerate(ens_pros):
        seq_pro = _sequence_from_ensembl_protein(enspro, protein=True)

        # validate if the sequence of uniprot and ensembl protein matches
        if _compare_sequences(seq, seq_pro, permissive=False):
            usable_indexes.append(i)
        elif _compare_sequences(seq, seq_pro, permissive=True):
            usable_indexes.append(i)
            seq_maps.append(None)
            n_mismatches = _count_mismatches(seq, seq_pro)
            message = "{0}: Sequences are of same length but have {1} mismatch(s)".format(enspro,
                                                                                          n_mismatches)
            logging.warning(message)
        elif align_transcripts:
            message = "Sequences don't match! Will attempt alignment... {}".format(enspro)
            logging.warning(message)
            aligned_indexes.append(i)
            ensembl_to_uniprot = _map_sequence_indexes(seq_pro, seq)
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
