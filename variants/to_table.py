#!/usr/bin/env python
# -*- coding: utf-8 -*-
import cPickle as pickle
import json
import logging
import os.path
import sys
from StringIO import StringIO

import pandas as pd
import requests

from config import defaults
from utils import (get_url_or_retry, is_valid, compare_sequences, count_mismatches,
                   fetch_uniprot_gff, raise_if_not_ok, map_sequence_indexes,
                   apply_sequence_index_map)

log = logging.getLogger(__name__)


##############################################################################
# Private methods
##############################################################################
def _fetch_uniprot_variants(identifier, _format='tab'):
    """
    Request human curated variants from UniProt.

    :param identifier: UniProt ID
    :param _format: request output from the webserver
    :return: pandas dataframe
    """
    # Check if query has already been saved

    query_file_name = defaults.db_variants + 'uniprot_variants_' + identifier + '.pkl'
    if not os.path.isfile(query_file_name):
        url = defaults.http_uniprot + '?query=accession:' + identifier
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


def _variant_characteristics_from_identifiers(variant_ids, use_vep=False):
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
        r = requests.post(url, headers=headers, data=data)

    # GET if given single id
    if isinstance(variant_ids, str):
        ensembl_endpoint = "variation/homo_sapiens/" + variant_ids
        if use_vep:
            ensembl_endpoint = "vep/human/id/" + variant_ids
        headers = {"Content-Type": "application/json"}
        params = {"phenotypes": 1}
        url = defaults.api_ensembl + ensembl_endpoint
        r = requests.get(url, headers=headers, params=params)

    if not r.ok:
        r.raise_for_status()
        sys.exit()

    decoded = r.json()
    return decoded


def _fetch_ensembl_variants(ensembl_ptn_id, feature=None):
    """Queries the Ensembl API for germline variants (mostly dbSNP) and somatic
    (mostly COSMIC) based on Ensembl Protein identifiers (e.g. ENSP00000326864).

    :param ensembl_ptn_id: Ensembl acession to a protein: ENSP00000XXXXXX
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
        # params = {'feature': suported_feats,
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


def _ensembl_variant(identifier, species='human'):
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


def _sequence_from_ensembl_protein(identifier, species='human', protein=True):
    """
    Gets the sequence for an Ensembl identifier.

    :param identifier: Ensembl ID
    :param species: Ensembl species
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


def _uniprot_to_ensembl_xref(uniprot_id, specie='homo_sapiens'):
    """Return Gene, transcripts and tranlational ids from Ensembl to Uniprot.
    Ensembl -> Uniprot reference if better than otherwise due Ensembl move quicker than Uniprot.

    :param uniprot_id:
    :param specie:
    :return:
    :rtype: pandas.DataFrame
    """
    # TODO add increase specie support
    url = "{}xrefs/symbol/{}/{}?content-type=application/json".format(
            defaults.api_ensembl, specie, uniprot_id)

    return pd.read_json(url)


def _match_uniprot_ensembl_seq(uniprot_id):
    """

    :param uniprot_id:
    :return:
    """
    uniprot_seq_url = "{}{}.fasta".format(defaults.http_uniprot, str(uniprot_id))
    uniprot_seq = ''.join(get_url_or_retry(uniprot_seq_url).splitlines()[1:])

    ensembl_xref = _uniprot_to_ensembl_xref(uniprot_id)
    # TODO is human only for now

    # ensembl of UNIPROT ACCESSION with more than one gene P04637, which is rare
    # ensembl_gene_id = ensembl_xref.loc[ensembl_xref.type == 'gene', 'id'].first()

    ensembl_ptn_ids = ensembl_xref.loc[ensembl_xref.type == 'translation', 'id']
    uniprot_sequence = _uniprot_info(uniprot_id, cols='sequence').iloc[0, 1]
    for ensembl_ptn_id in ensembl_ptn_ids:
        ensembl_ptn_seq = _sequence_from_ensembl_protein(ensembl_ptn_id)
        if compare_sequences(uniprot_sequence, ensembl_ptn_seq, permissive=False):
            return ensembl_ptn_id


def _uniprot_info(identifier, retry_in=(503, 500), cols=None, check_id=False):
    """
    Retrive uniprot information from the database.

    :param identifier: UniProt accession identifier
    :return: pandas table dataframe
    """

    if check_id and not is_valid(identifier, database='uniprot'):
        raise ValueError("{} is not a valid UniProt identifier.".format(identifier))

    if not cols:
        cols = ('entry name', 'reviewed', 'protein names', 'genes', 'organism',
                'sequence', 'length')
    elif isinstance(cols, str):
        cols = ('entry name', cols)

    params = {'query': 'accession:' + str(identifier),
              'columns': ",".join(cols),
              'format': 'tab',
              'contact': ""}
    url = defaults.http_uniprot
    response = get_url_or_retry(url=url, retry_in=retry_in, **params)
    try:
        data = pd.read_table(StringIO(response))
    except ValueError as e:
        log.error(e)
        data = response
    return data


##############################################################################
# Public methods
##############################################################################


def select_missense_variants(uniprot_id):
    """Fetches uniprot variants from Uniprot. These are human curated and are
    related to disease causing mutations.

    Replaces _fetch_uniprot_variants
    :param uniprot_id:
    :return: table with mutations
    :rtype: pandas.DataFrame
    """
    # get ensembl accessions

    return


def icgc_missense_variant(ensembl_gene_id):
    """Fetch a gene missense variants from ICGC data portal.

    :param ensembl_gene_id: ensembl gene accession
    :type ensembl_gene_id: str
    :return: DataFrame with one mutation per row.
    :rtype : pandas.DataFrame

    :Example:


        >>> table = icgc_missense_variant('ENSG00000012048')
        >>> table.loc[0, ['id', 'mutation', 'type', 'start']]
        id                          MU601299
        mutation                         G>A
        type        single base substitution
        start                       41245683
        Name: 0, dtype: object


    .. note:: ICGC doc https://dcc.icgc.org/docs/#!/genes/findMutations_get_8

    """
    base_url = "https://dcc.icgc.org/api/v1/genes/" + ensembl_gene_id + "/mutations/"
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


def select_uniprot_variants(identifier, align_transcripts=False):
    """
    Summarise variants for a protein in the UniProt

    :type align_transcripts: object
    :param identifier: UniProt ID
    :return: table with variants, rows are residues
    """
    # TODO FIX docstring
    # get organism and sequence for the provided identifier

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
        seq_pro = _sequence_from_ensembl_protein(enspro, org, protein=True)

        # validate if the sequence of uniprot and ensembl protein matches
        if compare_sequences(seq, seq_pro, permissive=False):
            usable_indexes.append(i)
        elif compare_sequences(seq, seq_pro, permissive=True):
            usable_indexes.append(i)
            seq_maps.append(None)
            n_mismatches = count_mismatches(seq, seq_pro)
            message = "{0}: Sequences are of same length but have {1} mismatch(s)".format(enspro,
                                                                                          n_mismatches)
            logging.warning(message)
        elif align_transcripts:
            message = "Sequences don't match! Will attempt alignment... {}".format(enspro)
            logging.warning(message)
            aligned_indexes.append(i)
            ensembl_to_uniprot = map_sequence_indexes(seq_pro, seq)
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

            var_table.start = apply_sequence_index_map(var_table.start, seq_maps[i])
            mut_table.start = apply_sequence_index_map(mut_table.start, seq_maps[i])

            if not var_table.empty:
                tables.append(var_table)
            if not var_table.empty:
                tables.append(mut_table)

    # to_unique = lambda series: series.unique()
    # return table.groupby('start').apply(to_unique)
    table = pd.concat(tables)
    return table


def map_gff_features_to_sequence(identifier, query_type='', drop_types=(
        'Helix', 'Beta strand', 'Turn', 'Chain')):
    """Remaps features in the uniprot gff file to the sequence.

    :param query_type:
    :param identifier: UniProt-SP accession
    :param drop_types: Annotation type to be dropped
    :return: table read to be joined to main table
    """

    def annotation_writter(gff_row):
        """
        Establish a set of rules to annotates uniprot GFF.

        :param gff_row: each line in the GFF file
        :return:
        """
        if not gff_row.ID and not gff_row.Note:
            return gff_row.TYPE
        elif not gff_row.ID:
            return '{0.TYPE}: {0.Note}'.format(gff_row)
        elif not gff_row.Note:
            return '{0.TYPE} ({0.ID})'.format(gff_row)
        else:
            return '{0.TYPE}: {0.Note} ({0.ID})'.format(gff_row)

    table = fetch_uniprot_gff(identifier)
    if query_type:
        table = table[table.TYPE == query_type]
    elif drop_types:
        table = table[~table.TYPE.isin(drop_types)]

    lines = []
    for i, row in table.iterrows():
        lines.extend({'idx': i, 'annotation': annotation_writter(row)}
                     for i in range(row.START, row.END + 1))
    table = pd.DataFrame(lines)
    return table.groupby('idx').agg({'annotation': lambda x: ', '.join(x)})


if __name__ == '__main__':
    pass
