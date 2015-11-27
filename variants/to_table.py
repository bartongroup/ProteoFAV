#!/usr/bin/env python
# -*- coding: utf-8 -*-


import sys
import requests
import logging
import pandas as pd
from StringIO import StringIO
from urlparse import parse_qs

from config import defaults
from utils import get_url_or_retry, is_valid, compare_uniprot_ensembl_sequence

log = logging.getLogger(__name__)

__all__ = ["select_uniprot_gff", "select_uniprot_variants"]


##############################################################################
# Private methods
##############################################################################
def _fetch_uniprot_variants(identifier, _format='tab'):
    """Request human curated variants from UniProt

    :param identifier:
    :return:
    """
    url = defaults.http_uniprot + '?query=accession:' + identifier
    url += '&format=' + _format
    url += '&columns=feature(NATURAL+VARIANT)'

    r = requests.get(url)
    if not r.ok:
        r.raise_for_status()
        sys.exit()

    # Complicated parsing
    records = r.content.replace('Natural variant\n', '').split('.; ')
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
    """Retrieves variant annotation from ENSEMBL

    :param variant_ids: Enseble Variant identifier
    :param use_vep: wether to use predicted variants fro VEP
    :return:
    """

    # POST if given a list of ids
    if isinstance(variant_ids, list):
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


def _transcript_variants_ensembl(identifier, missense=True, species=None):
    """Queries the Ensembl API for transcript variants (mostly dbSNP)
    based on Ensembl Protein identifiers (e.g. ENSP00000326864).

    :param identifier: Ensembl Protein ID
    :param missense: if True only fetches missense variants
    :return: pandas table dataframe
    """
    # TODO not using organism
    ensembl_endpoint = "overlap/translation/"
    if missense:
        params = {'feature': 'transcript_variation',
                  'type': 'missense_variant'}
    else:
        params = {'feature': 'transcript_variation'}
    url = defaults.api_ensembl + ensembl_endpoint + str(identifier)
    rows = get_url_or_retry(url, json=True, **params)
    return pd.DataFrame(rows)


def _somatic_variants_ensembl(identifier, missense=True, species=None):
    """Queries the Ensembl API for somatic transcript variants (COSMIC) based on
     Ensembl Protein identifiers (e.g. ENSP00000326864).

    :param identifier: Ensembl Protein ID
    :param missense: if True only fetches missense variants
    :return: pandas table dataframe
    """
    ensembl_endpoint = "overlap/translation/"
    if missense:
        params = {'feature': 'somatic_transcript_variation',
                  'type': 'missense_variant'}
    else:
        params = {'feature': 'somatic_transcript_variation'}
    url = defaults.api_ensembl + ensembl_endpoint + identifier
    rows = get_url_or_retry(url, json=True, **params)
    return pd.DataFrame(rows)


def _ensembl_variant(identifier, species='human'):
    """Queries the Ensembl API for variant IDs (e.g rs376845802 or COSM302853).

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
    """    Gets the sequence for an Ensembl identifier.

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
    sequence = get_url_or_retry(url, json=False, header=header, **params)
    return sequence


def _uniprot_ensembl_mapping(identifier, species='human'):
    """Uses the UniProt mapping service to try and get Ensembl IDs for the
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


def _uniprot_info(identifier, retry_in=(503, 500), cols=None):
    """
    Retrive uniprot information from the database.

    :param identifier: UniProt accession identifier
    :return: pandas table dataframe
    """

    if not is_valid(identifier, database='uniprot'):
        raise ValueError(
            "{} is not a valid UniProt identifier.".format(identifier))

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
        log.errore(e)
        data = response
    return data


def _uniprot_gff(identifier):
    """Retrive UniProt data from the GFF file

    :param identifier: UniProt accession identifier
    :return: pandas table
    """
    url = defaults.http_uniprot + identifier + ".gff"
    cols = "NAME SOURCE TYPE START END SCORE STRAND FRAME GROUP empty".split()

    data = pd.read_table(url, skiprows=2, names=cols)
    groups = data.GROUP.apply(parse_qs)
    groups = pd.DataFrame.from_records(groups)
    data = data.merge(groups, left_index=True, right_index=True)

    return data


##############################################################################
# Public methods
##############################################################################
def select_uniprot_gff(identifier, drop_types=('Helix', 'Beta strand', 'Turn',
                                               'Chain')):
    """Sumarieses the GFF file for a UniProt acession.
    :param identifier: UniProt-SP acession
    :param drop_types: Annotation type to be droped
    :return: table read to be joined to main table
    """

    def annotation_writter(gff_row):
        """ Establish a set of rules to annotates uniprot GFF
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

    if not drop_types:
        drop_types = []
    data = _uniprot_gff(identifier)
    data = data[~data.TYPE.isin(drop_types)]

    lines = []
    for i, row in data.iterrows():
        lines.extend({'idx': i, 'annotation': annotation_writter(row)}
                     for i in range(row.START, row.END + 1))
    data = pd.DataFrame(lines)
    return data.groupby('idx').agg({'annotation': lambda x: ', '.join(x)})


def select_uniprot_variants(identifier):
    """Sumarise variants for a protein in the UniProt

    :param identifier: UniProt ID
    :return: table with variants, rows are residues
    """
    # get organism and sequence for the provided ientifier

    uni = _uniprot_gff(identifier, cols=['organism', 'sequence'])
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
    for i, enspro in enumerate(ens_pros):
        seq_pro = _sequence_from_ensembl_protein(enspro, org, protein=True)

        # validate if the sequence of uniprot and ensembl protein matches
        if compare_uniprot_ensembl_sequence(seq, seq_pro, permissive=False):
            usable_indexes.append(i)
        else:
            message = "Sequences don't match! skipping... {}".format(enspro)
            logging.warning(message)

    # get the variants for the ensembl proteins that match the uniprot
    tables = []
    for i in usable_indexes:
        vars = _transcript_variants_ensembl(ens_pros[i], missense=True)
        muts = _somatic_variants_ensembl(ens_pros[i], missense=True)

        # TODO: From... TO... mutated residues as different columns in the table

        tables.append(vars[['translation', 'id', 'start', 'residues']])
        tables.append(muts[['translation', 'id', 'start', 'residues']])

    # to_unique = lambda series: series.unique()
    # return table.groupby('start').apply(to_unique)
    table = pd.concat(tables)
    return table


if __name__ == '__main__':
    pass
