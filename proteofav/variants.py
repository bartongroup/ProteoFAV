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

__all__ = ["_fetch_icgc_variants",
           "_fetch_ebi_variants",
           "_fetch_ensembl_variants",
           "_sequence_from_ensembl_protein",
           "_uniprot_ensembl_mapping",
           "_match_uniprot_ensembl_seq",
           "_apply_sequence_index_map",
           "_compare_sequences",
           "_count_mismatches",
           "fetch_uniprot_sequence",
           "fetch_uniprot_formal_specie",
           "select_variants",
           "parse_uniprot_variants",
           "select_uniprot_variants"]


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


def _fetch_ebi_variants(uniprot_idd, flat_xrefs=True):
    """
    Fetchs the variant data from EBI. This datasource is also used in the UniProt feature viewer
    :param bool flat_xrefs: whether to parse or not the cross-reference field that indicates
    data provenance.
    :param uniprot_idd: UniProt accession
    :return pd.DataFrame: the table where each row is associated with a variation entry
    .. note::
    if flat_xrefs is true multiple rows are produced with the same index to
    """
    endpoint = "variation/"
    url = defaults.api_proteins + endpoint + uniprot_idd

    data = fetch_from_url_or_retry(url, json=True).json()
    data = json_normalize(data, ['features'], meta=['accession', 'entryName'])

    # flatten the xref field, which has the id column.
    # ideally this could be normalised with:
    # table = json_normalize(data, ['features', 'xref'] ...
    # but this field is not present in all entries,
    if flat_xrefs:
        flat_xref = data['xrefs'].apply(pd.Series).stack().apply(pd.Series)
        flat_xref.reset_index(level=1, drop=True, inplace=True)
        data = data.join(flat_xref)
    return data


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
        raise NotImplementedError('feature argument should be one of {} '.format(', '''.join(
            supported_feats)))
    else:
        params = {'feature': feature}
    url = defaults.api_ensembl + ensembl_endpoint + ensembl_ptn_id

    rows = fetch_from_url_or_retry(url, json=True, **params).json()
    return pd.DataFrame(rows)


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

    return fetch_from_url_or_retry(url, header=header, **params)


def _uniprot_ensembl_mapping(identifier, species=None):
    """
    Uses the UniProt mapping service to try and get Ensembl IDs for the
    UniProt accession identifier provided

    :param identifier: UniProt accession identifier
    :param species: Ensembl species
    :return: pandas table dataframe
    """

    if species not in valid_ensembl_species:
        raise ValueError('Provided species {} is not valid'.format(species))

    information = {}
    rows = []

    ensembl_endpoint = "xrefs/symbol/{}/".format(species)
    url = defaults.api_ensembl + ensembl_endpoint + str(identifier)
    data = fetch_from_url_or_retry(url, json=True).json()
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
        log.error('Could not retrieve {} information. Maybe obsolete UniProt accession?'.format(
            uniprot_id))
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


def _uniprot_to_ensembl_xref(uniprot_id, species='homo_sapiens'):
    """
    Return Gene, transcripts and translational ids from Ensembl to Uniprot.
    Ensembl -> Uniprot reference is ideal because Ensembl database change more
    often the Uniprot'smove quicker than Uniprot.

    :param str uniprot_id: Uniprot accession
    :param str species: species name
    :return pandas.DataFrame: table with columns
    """

    url = "{}xrefs/symbol/{}/{}?content-type=application/json".format(
        defaults.api_ensembl, species, uniprot_id)

    return pd.read_json(url)


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


def _pdb_uniprot_sifts_mapping(identifier):
    """
    Queries the PDBe API for SIFTS mapping between PDB - UniProt.
    One to many relationship expected.

    :param identifier: PDB id
    :return: pandas table dataframe
    """

    sifts_endpoint = 'mappings/uniprot/'
    url = defaults.api_pdbe + sifts_endpoint + identifier
    information = fetch_from_url_or_retry(url, json=True).json()

    rows = []
    for uniprot in information[identifier]['UniProt']:
        uniprots = {'uniprot_id': uniprot}
        rows.append(uniprots)
    return pd.DataFrame(rows)


def _uniprot_pdb_sifts_mapping(identifier):
    """
    Queries the PDBe API for SIFTS mapping between UniProt - PDB entries.
    One to many relationship expected.

    :param identifier: UniProt ID
    :return: pandas table dataframe
    """
    sifts_endpoint = 'mappings/best_structures/'
    url = defaults.api_pdbe + sifts_endpoint + str(identifier)
    information = fetch_from_url_or_retry(url, json=True).json()

    rows = []
    for entry in information[identifier]:
        rows.append(entry)
    return pd.DataFrame(rows)


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
