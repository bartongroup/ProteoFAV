# -*- coding: utf-8 -*-

import json
import logging
import requests
import pandas as pd
from io import BytesIO
from io import StringIO
from pandas.io.json import json_normalize

from proteofav.annotation import select_annotation
from proteofav.utils import fetch_from_url_or_retry, row_selector
from proteofav.library import valid_ensembl_species

from proteofav.config import defaults

log = logging.getLogger('proteofav.config')

__all__ = ["fetch_uniprot_variants",
           "fetch_ensembl_variants",
           "fetch_ensembl_sequence_from_id",
           "fetch_uniprot_ensembl_mapping",
           "fetch_ensembl_uniprot_mapping",
           "fetch_uniprot_species_from_id",
           "fetch_pdb_uniprot_mapping",
           "fetch_uniprot_pdb_mapping",
           "get_ensembl_species_from_uniprot",
           "get_uniprot_id_from_mapping",
           "get_ensembl_protein_id_from_mapping",
           "get_preferred_uniprot_id_from_mapping",
           "get_preferred_ensembl_id_from_mapping",
           "_match_uniprot_ensembl_seq",
           "_apply_sequence_index_map",
           "_compare_sequences",
           "_count_mismatches",
           "fetch_uniprot_sequence",
           "fetch_uniprot_formal_specie",
           "select_variants",
           "parse_uniprot_variants",
           "select_uniprot_variants",
           "icgc_missense_variant",
           "_fetch_icgc_variants"]


def fetch_uniprot_variants(identifier, retry_in=(429,)):
    """
    Queries the EBI Proteins API for variants based on
    UniProt identifiers (e.g. O15294).

    JSON structure flattening is sorted out in another helper
    function.

    :param identifier: UniProt ID
    :param retry_in: http code for retrying connections
    :return: Requests response object
    """

    url_root = defaults.api_proteins
    url_endpoint = "variation/{}".format(identifier)
    url = url_root + url_endpoint
    response = fetch_from_url_or_retry(url=url, json=True, retry_in=retry_in)
    return response


def fetch_ensembl_variants(identifier, feature=None, retry_in=(429,)):
    """
    Queries the Ensembl API for germline variants (mostly dbSNP) and somatic
    (mostly COSMIC) based on Ensembl Protein identifiers (e.g. ENSP00000326864).

    :param identifier: Ensembl Protein ID
    :param feature: variation type
    :param retry_in: http code for retrying connections
    :return: Requests response object
    """

    url_root = defaults.api_ensembl
    url_endpoint = "overlap/translation/{}".format(identifier)
    supported_feats = ['transcript_variation', 'somatic_transcript_variation']
    if feature is None:
        raise NotImplementedError('Use two functions call to get both somatic'
                                  ' and germline variants.')
    elif feature not in supported_feats:
        raise NotImplementedError("Feature argument should be one of '{}'"
                                  "".format("', '".join(supported_feats)))
    else:
        params = {'feature': feature}
        # params = {'feature': feature,
        #           'type': 'missense_variant'}
    url = url_root + url_endpoint
    response = fetch_from_url_or_retry(url, json=True, retry_in=retry_in, **params)
    return response


def fetch_ensembl_sequence_from_id(identifier, protein=True, retry_in=(429,)):
    """
    Queries the Ensembl REST API for the sequence of a Ensembl Protein ID.

    :param identifier: Ensembl Protein ID
    :param protein: boolean
    :param retry_in: http code for retrying connections
    :return: Requests response object
    """
    url_root = defaults.api_ensembl
    url_endpoint = "sequence/id/{}".format(identifier)
    url = url_root + url_endpoint
    if protein:
        params = {'type': 'protein'}
    else:
        params = {}
    response = fetch_from_url_or_retry(url=url, json=True,
                                       retry_in=retry_in, **params)
    return response


def fetch_uniprot_ensembl_mapping(identifier, species='homo_sapiens', retry_in=(429,)):
    """
    Uses the Ensembl REST mapping service to try and get Ensembl IDs for
    the UniProt accession identifier provided.

    :param identifier: UniProt accession identifier
    :param retry_in: http code for retrying connections
    :param species: Ensembl species
    :return: Requests response object
    """

    if species not in valid_ensembl_species:
        raise ValueError('Provided species {} is not valid'.format(species))

    url_root = defaults.api_ensembl
    url_endpoint = "xrefs/symbol/{}/{}".format(species, identifier)
    url = url_root + url_endpoint

    response = fetch_from_url_or_retry(url=url, json=True, retry_in=retry_in)
    return response


def fetch_ensembl_uniprot_mapping(identifier, retry_in=(429,)):
    """
    Uses the Ensembl REST mapping service to try and get UniProt IDs for
    the Ensembl Protein accession identifier provided.

    :param identifier: Ensembl Protein accession identifier
    :param retry_in: http code for retrying connections
    :return: Requests response object
    """

    url_root = defaults.api_ensembl
    url_endpoint = "xrefs/id/{}".format(identifier)
    url = url_root + url_endpoint
    # params = {'external_db': "Uniprot/SWISSPROT"}
    # params = {'external_db': "Uniprot/SPTREMBL"}

    response = fetch_from_url_or_retry(url=url, json=True, retry_in=retry_in)
    return response


def fetch_uniprot_species_from_id(identifier, retry_in=(429,)):
    """
    Retrieve Species from UniProt ID.

    :param identifier: UniProt accession identifier
    :param retry_in: http code for retrying connections
    :return: Requests response object
    """

    url_root = defaults.api_uniprot
    url_endpoint = "?query={}&columns=organism&format=tab".format(identifier)
    url = url_root + url_endpoint
    response = fetch_from_url_or_retry(url=url, json=False, retry_in=retry_in)
    return response


def fetch_pdb_uniprot_mapping(identifier, retry_in=(429,)):
    """
    Queries the PDBe API for SIFTS mapping between PDB - UniProt.
    One to many relationship expected.

    :param identifier: PDB accession identifier
    :param retry_in: http code for retrying connections
    :return: Requests response object
    """

    sifts_endpoint = 'mappings/uniprot/{}'.format(identifier)
    url = defaults.api_pdbe + sifts_endpoint
    response = fetch_from_url_or_retry(url, json=True, retry_in=retry_in)
    return response


def fetch_uniprot_pdb_mapping(identifier, retry_in=(429,)):
    """
    Queries the PDBe API for SIFTS mapping between UniProt - PDB entries.
    One to many relationship expected.

    :param identifier: UniProt accession identifier
    :param retry_in: http code for retrying connections
    :return: Requests response object
    """
    sifts_endpoint = 'mappings/best_structures/{}'.format(identifier)
    url = defaults.api_pdbe + sifts_endpoint
    response = fetch_from_url_or_retry(url, json=True, retry_in=retry_in)
    return response


def get_ensembl_species_from_uniprot(data):
    """
    Gets a Species Name from a UniProt organism lookup.

    :param data: Request object from the UniProt Query endpoint.
    :return: (str) formatted species name
    """
    organism = str(data.content, encoding='utf-8').split('\n')[1]
    species = '_'.join(organism.split()[0:2]).lower()
    return species


def get_ensembl_protein_id_from_mapping(data):
    """
    Gets a list of Ensembl IDs from a 'xrefs/symbol/' mapping.

    :param data: Requests object from the Ensembl-UniProt Mapping
    :return: list of Ensembl Protein IDs
    """
    ensps = []
    for entry in data:
        if 'type' in entry and 'id' in entry:
            if entry['type'] == 'translation':
                if entry['id'] not in ensps:
                    ensps.append(entry['id'])
    return ensps


def get_uniprot_id_from_mapping(data, full_entry=False, uniprot_id=None):
    """
    Gets a list of UniProt IDs from a '"xrefs/id/"' mapping.

    :param data: Requests object from the Ensembl-UniProt Mapping
    :param full_entry: (boolean) if True gets dictionary instead of just
        the UniProt IDs
    :param uniprot_id: if not None means that we want the data for a specific
        UniProt ID
    :return: list of UniProt IDs
    """
    uniprots = []
    for entry in data:
        if 'dbname' in entry and 'primary_id' in entry:
            if uniprot_id is not None and entry['primary_id'] == uniprot_id:
                if full_entry:
                    uniprots = [entry]
                else:
                    uniprots = [entry['primary_id']]
                break
            elif entry['dbname'] == 'Uniprot/SWISSPROT':
                if entry['primary_id'] not in uniprots:
                    if full_entry:
                        uniprots.append(entry)
                    else:
                        uniprots.append(entry['primary_id'])
            elif entry['dbname'] == 'Uniprot/SPTREMBL':
                if entry['primary_id'] not in uniprots:
                    if full_entry:
                        uniprots.append(entry)
                    else:
                        uniprots.append(entry['primary_id'])
    return uniprots


def get_preferred_uniprot_id_from_mapping(data):
    """
    Takes a list of Ensembl xrefs/ids mapped from a UniProt ID
    and gets the preferred entry (Many-to-one), based on seq
    identity and coverage.

    :param data: list of dictionaries
    :return: (str) preferred UniProt ID
    """

    best_match = None
    curr_ix = -1
    prev_identity = 0
    prev_coverage = 0
    prev_id = "-" * 100
    for ix, entry in enumerate(data):
        if ('ensembl_identity' in entry and 'xref_identity' in entry and
                    'xref_start' in entry and 'xref_end' in entry):
            identity = entry['ensembl_identity'] + entry['xref_identity']
            coverage = entry['xref_end'] - entry['xref_start']
            if identity + coverage >= prev_identity + prev_coverage:
                prev_identity = identity
                prev_coverage = coverage
                # preferring the smallest UniProt ID (for getting variants)
                if len(entry['primary_id']) < len(prev_id):
                    prev_id = entry['primary_id']
                    curr_ix = ix
    if curr_ix != -1 and 'primary_id' in data[curr_ix]:
        best_match = data[curr_ix]['primary_id']
    return best_match


def get_preferred_ensembl_id_from_mapping(identifiers, uniprot_id=None):
    """
    Takes a list of Ensembl xrefs/ids mapped from a UniProt ID
    and gets the preferred entry (Many-to-one), based on seq
    identity and coverage.

    :param identifiers: list of Ensembl IDs
    :param uniprot_id: if not None means that we want the data for a specific
        UniProt ID
    :return: (str) preferred Ensembl ID
    """

    best_match = None
    curr_ix = -1
    prev_identity = 0
    prev_coverage = 0
    for ix, ensp in enumerate(identifiers):
        info = fetch_ensembl_uniprot_mapping(ensp).json()
        # gets the mapping for a specific uniprot
        data = get_uniprot_id_from_mapping(info, full_entry=True,
                                           uniprot_id=uniprot_id)
        for entry in data:
            if ('ensembl_identity' in entry and 'xref_identity' in entry and
                        'xref_start' in entry and 'xref_end' in entry):

                identity = entry['ensembl_identity'] + entry['xref_identity']
                coverage = entry['xref_end'] - entry['xref_start']
                if identity + coverage > prev_identity + prev_coverage:
                    prev_identity = identity
                    prev_coverage = coverage
                    curr_ix = ix
    if curr_ix != -1:
        best_match = identifiers[curr_ix]
    return best_match


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
    ensembl_xref = pd.DataFrame(fetch_uniprot_ensembl_mapping(uniprot_id.split('-')[0],
                                                              species=species).json())

    if ensembl_xref.empty:
        raise ValueError('No cross-reference found for {} Ensembl mapping.'.format(
            uniprot_id))

    ensembl_ptn_ids = ensembl_xref.loc[ensembl_xref.type == 'translation', 'id']

    uniprot_sequence = fetch_uniprot_sequence(uniprot_id)
    for ensembl_ptn_id in ensembl_ptn_ids:
        ensembl_ptn_seq = fetch_ensembl_sequence_from_id(ensembl_ptn_id)
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


def fetch_uniprot_sequence(uniprot_id):
    """
    Gets current sequence of a Uniprot entry.

    :param str uniprot_id: Uniprot accession
    :return str: the sequence

    >>> print(fetch_uniprot_formal_specie('P17612'))
    Homo sapiens
    """

    return _uniprot_info(uniprot_id, cols='sequence').iloc[0, 1]


def fetch_uniprot_formal_specie(uniprot_id, remove_isoform=True):
    """
    Gets the species name of an organism expressing a protein.

    :param Bool remove_isoform: whether to remove the isoform identifier.
    :param str uniprot_id: Uniprot accession
    :return: the species name (two words)
    :rtype: str or None

    >>> print(fetch_uniprot_sequence('P17612'))[:20]
    MGNAAAAKKGSEQESVKEFL
    """
    if remove_isoform:
        uniprot_id = uniprot_id.split('-')[0]

    full_specie = _uniprot_info(uniprot_id, cols='organism').iloc[0, 1]

    try:
        return " ".join(full_specie.split()[0:2])
    except AttributeError:
        log.error('Could not retrieve {} information. Maybe obsolete '
                  'UniProt accession?'.format(uniprot_id))
        return None


def _uniprot_info(uniprot_id, retry_in=(503, 500), cols=None):
    """
    Retrieve Uniprot information from the database.

    :param str uniprot_id: Uniprot accession identifier
    :param retry_in: iterable of status_code to be retried upon error.
    :type retry_in: list of [int]
    :return pandas.DataFrame: table from Uniprot.
    Default table columns:
    :raises requests.HTTPError: when hits a bad status_code
    """

    if not cols:
        cols = ('id', 'entry name', 'reviewed', 'protein names', 'genes', 'organism',
                'sequence', 'length')
    elif isinstance(cols, str):
        cols = ('id', cols)

    params = {'query': 'accession:' + str(uniprot_id),
              'columns': ",".join(cols),
              'format': 'tab',
              'contact': ""}
    url = defaults.api_uniprot
    response = fetch_from_url_or_retry(url=url, retry_in=retry_in, **params).content
    try:
        data = pd.read_table(StringIO(response))
    except TypeError:
        # python 3.5
        data = pd.read_table(BytesIO(response))
    except ValueError as e:
        log.error(e)
        return None
    # id column is called Entry in the table
    return row_selector(data, 'Entry', uniprot_id)


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
        table_som = fetch_ensembl_variants(
            ensembl_ptn_id, feature='somatic_transcript_variation')

        if not table_som.empty:
            table_som = table_som[table_som['type'] == 'coding_sequence_variant']
            table_som = table_som.groupby('start')['id'].apply(list)
            table_som.name = 'somatic_variants'
        else:
            table_som['somatic_variants'] = None

        tables.append(table_som)

    if 'ensembl_germline' in features:
        table_ger = fetch_ensembl_variants(ensembl_ptn_id, feature='transcript_variation')

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

    table = select_annotation(uniprot_id, annotation_agg=True,
                              query_type='Natural variant')

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
    ens = fetch_uniprot_ensembl_mapping(identifier, species=org)

    # get the ensembl protein ids
    ens_pros = ens.loc[0, 'TRANSLATION']
    if not isinstance(ens_pros, list):
        ens_pros = [ens_pros, ]

    # get the sequence of the ensembl protein
    usable_indexes = []
    aligned_indexes = []
    seq_maps = []

    for i, enspro in enumerate(ens_pros):
        seq_pro = fetch_ensembl_sequence_from_id(enspro, protein=True)

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
        vars = fetch_ensembl_variants(ens_pros[i], feature='transcript_variation')
        muts = fetch_ensembl_variants(ens_pros[i], feature='somatic_transcript_variation')

        # TODO: From... TO... mutated residues as different columns in the table
        # TODO: Shouldn't the default behaviour return all the columns?
        if not vars.empty:
            tables.append(vars[['translation', 'id', 'start', 'residues']])
        if not muts.empty:
            tables.append(muts[['translation', 'id', 'start', 'residues']])

    # Get variants from aligned sequences
    if align_transcripts:
        for i in aligned_indexes:
            vars = fetch_ensembl_variants(ens_pros[i], feature='transcript_variation')
            muts = fetch_ensembl_variants(ens_pros[i], feature='somatic_transcript_variation')

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


def raise_if_not_ok(response):
    """

    :param response:
    """
    if not response.ok:
        response.raise_for_status()


def icgc_missense_variant(ensembl_gene_id):
    """Fetch a gene missense variants from ICGC data portal.

    :param str ensembl_gene_id: ensembl gene accession
    :return pd.DataFrame: DataFrame with one mutation per row.
    :example:

        >>> table = icgc_missense_variant('ENSG00000012048')
        >>> table.loc[0, ['id', 'mutation', 'type', 'start']]
        id                          MU601299
        mutation                         G>A
        type        single base substitution
        start                       41245683
        Name: 0, dtype: object

    .. note:: ICGC doc https://dcc.icgc.org/docs/#!/genes/findMutations_get_8
    """
    base_url = defaults.api_icgc + ensembl_gene_id + "/mutations/"
    headers = {'content-type': 'application/json'}
    filt = json.dumps({"mutation": {"consequenceType": {"is": ['missense_variant']}}})
    params = {"filters": filt}

    # First counts the number of mutation entries
    counts_resp = requests.get(base_url + "counts/", headers=headers, params=params)
    raise_if_not_ok(counts_resp)
    total = counts_resp.json()['Total']

    # then iterate the pages of entries, max 100 entries per page
    hits = []
    params['size'] = 100
    for i in range(total // 100 + 1):
        params['from'] = i * 100 + 1
        mutation_resp = requests.get(base_url, headers=headers, params=params)
        raise_if_not_ok(mutation_resp)
        hits.extend(mutation_resp.json()['hits'])

    return pd.DataFrame(hits)


def _fetch_icgc_variants(identifier):
    """
    Queries the ICGC data portal for the PanCancer variants based on Ensembl
    transcript identifier.
    :param str identifier: Ensembl transcript
    :return pandas.DataFrame: pandas table dataframe
    """
    transition_regex = '(?P<ref>[A-Z])(?P<position>[0-9]+)(?P<new>[A-Z\*])?'
    variation_endpoint = "protein/"
    url = defaults.api_icgc + variation_endpoint + identifier
    # fetches a nested json
    # normalise the data, making it flat
    data = fetch_from_url_or_retry(url, json=True).json()
    data = pd.io.json.json_normalize(
        data['hits'],
        ['transcripts'],
        ['affectedDonorCountTotal', 'id', 'mutation'], meta_prefix='_')
    data = data[data['id'] == identifier]
    data.drop(['id'], axis=1, inplace=True)
    data.rename(columns={'_id': 'id'}, inplace=True)

    consequence = data.pop('consequence')
    if consequence.index.duplicated().any():  # pragma: no cover
        log.warning('Joining ICGC variant data with its consequences data aborted:'
                    ' Duplicated index for {}.'.format(identifier))
    else:
        data = data.join(consequence.apply(pd.Series), rsuffix='_protein')
        transition = data.aaMutation.str.extract(transition_regex)
        data = data.join(transition)

    return data

