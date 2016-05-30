#!/usr/bin/env python
# -*- coding: utf-8 -*-
import logging

import pandas as pd
import requests
from pandas.io.json import json_normalize

from proteofav.config import defaults
from proteofav.library import valid_ensembl_species
from proteofav.uniprot import (map_gff_features_to_sequence, _uniprot_to_ensembl_xref,
                               fetch_uniprot_formal_specie, fetch_uniprot_sequence)
from proteofav.utils import (get_url_or_retry, check_local_or_fetch)

log = logging.getLogger(__name__)


##############################################################################
# Private methods
##############################################################################
def _fetch_icgc_variants(identifier):
    """
    Queries the ICGC data portal for the PanCancer variants based on Ensembl
    transcript identifier.
    :param str identifier: Ensembl transcript
    :return pandas.DataFrame: pandas table dataframe
    """
    transition_regex = '(?P<ref>[A-Z])(?P<position>[0-9]+)(?P<new>[A-Z])'
    variation_endpoint = "protein/"
    url = defaults.api_icgc + variation_endpoint + identifier
    # fetches a nested json
    data = get_url_or_retry(url, json=True)
    # normalise the data, making it flat
    data = pd.io.json.json_normalize(
        data['hits'],
        ['transcripts'],
        ['affectedDonorCountTotal', 'id', 'mutation'], meta_prefix='_')

    data = data[data.id == identifier]
    consequence = data.pop('consequence')
    if consequence.index.duplicated().any():
        log.warn('Joining ICGC variant data with its consequences data aborted:'
                 ' Duplicated index for {}.'.format(identifier))
    else:
        data = data.join(consequence.apply(pd.Series), rsuffix='_protein')
        transition = data.aaMutation.str.extract(transition_regex)
        data = data.join(transition)

    return data


def _fetch_ebi_variants(identifier, flat_xrefs=True):
    """Fetchs the varition data from EBI also used in the UniProt feature viwer
    :param identifier: UniProt identifier
    :return pandas.DataFrame: the table where each row is associated with a variation entry
    """
    endpoint = "variation/"
    url = defaults.api_ebi_uniprot + endpoint + identifier
    data = get_url_or_retry(url, json=True)
    table = json_normalize(data, ['features'], meta=['accession', 'entryName'])

    # flatten the xref field, which has the id column.
    # ideally this could be normalised with:
    # table = json_normalize(data, ['features', 'xref'] ...
    # but this field is not present in all entries,
    if flat_xrefs:
        flat_xref = table['xrefs'].apply(pd.Series).stack().apply(pd.Series)
        flat_xref.reset_index(level=1, drop=True, inplace=True)
        table = table.join(flat_xref)
    return table


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
    supported_feats = ['transcript_variation', 'somatic_transcript_variation']
    if feature is None:
        raise NotImplementedError('Use two functions call to get both somatic'
                                  ' and germline variants.')
        # params = {'feature': supported_feats,
        #           'type': 'missense_variant'}
    elif feature not in supported_feats:
        raise NotImplementedError('feature argument should be one of {} or None for all'
                                  ''.format(', '''.join(supported_feats)))
    else:
        params = {'feature': feature}
    url = defaults.api_ensembl + ensembl_endpoint + ensembl_ptn_id

    rows = get_url_or_retry(url, json=True, **params)
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


def _uniprot_ensembl_mapping(identifier, species=None):
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
        ensembl_ptn_seq = _sequence_from_ensembl_protein(ensembl_ptn_id)
        if _compare_sequences(uniprot_sequence, ensembl_ptn_seq, permissive=False):
            return ensembl_ptn_id
    raise ValueError('No protein with the same sequence was retrivied from Ensembl {}'.format(
        uniprot_id))


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
    check_local_or_fetch(uniprot_id)

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


if __name__ == '__main__':
    pass
