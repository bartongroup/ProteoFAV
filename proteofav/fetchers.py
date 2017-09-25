# -*- coding: utf-8

import logging

from proteofav.utils import fetch_from_url_or_retry
from proteofav.library import valid_ensembl_species_variation

from proteofav.config import defaults as config

log = logging.getLogger('proteofav')


class InvalidEnsemblSpecies(ValueError):
    pass


class Fetcher(object):
    def __init__(self, url, **kwargs):
        """
        :param url: (str) Full web-address
        :param cached: (boolean) if True, stores a pickle file locally
        :param cache_output: (str) file name if 'cached=True'
        """

        self.url = url
        self.kwargs = kwargs
        self.response = None
        self._fetch()

    def _fetch(self):
        self.response = fetch_from_url_or_retry(self.url, **self.kwargs)
        return self.response


def _fetch_summary_properties_pdbe(identifier, retry_in=(429,)):
    """
    Queries the PDBe API to get summary properties.

    :param identifier: PDB ID
    :param retry_in: http code for retrying connections
    :return: Requests response object
    """

    url_root = config.api_pdbe
    url_endpoint = "pdb/entry/summary/"
    url = url_root + url_endpoint + identifier
    b = Fetcher(url=url, json=True, retry_in=retry_in)
    return b.response


def _fetch_uniprot_id_from_name(identifier, retry_in=(429,)):
    """
    Retrieve UniProt ID from Name.

    :param identifier: UniProt accession identifier
    :param retry_in: http code for retrying connections
    :return: Requests response object
    """

    url_root = config.http_uniprot
    url_endpoint = "?query={}&columns=id&format=list".format(identifier)
    url = url_root + url_endpoint
    b = Fetcher(url=url, json=False, retry_in=retry_in)
    return b.response


def _fetch_uniprot_species_from_id(identifier, retry_in=(429,)):
    """
    Retrieve Species from UniProt ID.

    :param identifier: UniProt accession identifier
    :param retry_in: http code for retrying connections
    :return: Requests response object
    """

    url_root = config.http_uniprot
    url_endpoint = "?query={}&columns=organism&format=tab".format(identifier)
    url = url_root + url_endpoint
    b = Fetcher(url=url, json=False, retry_in=retry_in)
    return b.response


def _fetch_uniprot_variants_ebi(identifier, retry_in=(429,)):
    """
    Queries the EBI Proteins API for variants.
    based on UniProt identifiers (e.g. O15294).

    :param identifier: UniProt ID
    :param retry_in: http code for retrying connections
    :return: Requests response object
    """

    url_root = config.api_proteins
    url_endpoint = "variation/"
    url = url_root + url_endpoint + identifier
    b = Fetcher(url=url, json=True, retry_in=retry_in)
    return b.response


def _fetch_ensembl_uniprot_ensembl_mapping(identifier, retry_in=(429,),
                                           species='homo_sapiens'):
    """
    Uses the Ensembl REST mapping service to try and get Ensembl IDs for
    the UniProt accession identifier provided.

    :param identifier: UniProt accession identifier
    :param retry_in: http code for retrying connections
    :param species: Ensembl species
    :return: Requests response object
    """

    if species not in valid_ensembl_species_variation:
        raise InvalidEnsemblSpecies(
            'Provided species {} is not valid'.format(species))

    url_root = config.api_ensembl
    url_endpoint = "xrefs/symbol/{}/".format(species)
    url = url_root + url_endpoint + identifier

    b = Fetcher(url=url, json=True, retry_in=retry_in)
    return b.response


def _fetch_ensembl_ensembl_uniprot_mapping(identifier, retry_in=(429,)):
    """
    Uses the Ensembl REST mapping service to try and get UniProt IDs for
    the Ensembl Protein accession identifier provided.

    :param identifier: Ensembl Protein accession identifier
    :param retry_in: http code for retrying connections
    :return: Requests response object
    """

    url_root = config.api_ensembl
    url_endpoint = "xrefs/id/"
    url = url_root + url_endpoint + identifier
    # params = {'external_db': "Uniprot/SWISSPROT"}
    # params = {'external_db': "Uniprot/SPTREMBL"}

    b = Fetcher(url=url, json=True, retry_in=retry_in)
    return b.response


def _fetch_ensembl_transcript_variants(identifier, retry_in=(429,)):
    """
    Queries the Ensembl REST API for transcript variants (mostly from dbSNP).
    based on Ensembl Protein identifiers (e.g. ENSP00000275603).

    :param identifier: Ensembl Protein ID
    :param retry_in: http code for retrying connections
    :return: Requests response object
    """

    url_root = config.api_ensembl
    url_endpoint = "overlap/translation/"
    url = url_root + url_endpoint + identifier
    # params = {'feature': 'transcript_variation', 'type': 'missense_variant'}
    params = {'feature': 'transcript_variation'}
    b = Fetcher(url=url, json=True, retry_in=retry_in, **params)
    return b.response


def _fetch_ensembl_somatic_variants(identifier, retry_in=(429,)):
    """
    Queries the Ensembl REST API for somatic variants (mostly from COSMIC).
    based on Ensembl Protein identifiers (e.g. ENSP00000275603).

    :param identifier: Ensembl Protein ID
    :param retry_in: http code for retrying connections
    :return: Requests response object
    """

    url_root = config.api_ensembl
    url_endpoint = "overlap/translation/"
    url = url_root + url_endpoint + identifier
    params = {'feature': 'somatic_transcript_variation'}
    b = Fetcher(url=url, json=True, retry_in=retry_in, **params)
    return b.response


def _fetch_ensembl_variants_by_id(identifier, retry_in=(429,),
                                  species='homo_sapiens'):
    """
    Queries the Ensembl API for variant IDs (e.g rs376845802 or COSM302853).

    :param identifier: variant ID
    :param retry_in: http code for retrying connections
    :param species: Ensembl species
    :return: Requests response object
    """

    if isinstance(identifier, str):
        url_root = config.api_ensembl
        url_endpoint = "variation/{}/".format(species)
        url = url_root + url_endpoint + identifier
        # params = {'pops': '1', 'phenotypes': '1'} #, 'genotypes': '1', 'population_genotypes': '1'}
        b = Fetcher(url=url, json=True, retry_in=retry_in)
        return b.response
    elif isinstance(identifier, list):
        data = """{ "ids" : ["%s"]}""" % ('", "'.join(identifier))
        header = {"Accept": "application/json"}
        url_root = config.api_ensembl
        url_endpoint = "variation/{}".format(species)
        url = url_root + url_endpoint
        # params = {'pops': '1', 'phenotypes': '1'} #, 'genotypes': '1'}
        b = Fetcher(url=url, json=True, retry_in=retry_in,
                    post=True, data=data, header=header)
        return b.response


def _fetch_ensembl_sequence_from_id(identifier, retry_in=(429,)):
    """
    Queries the Ensembl REST API for the sequence of a Ensembl Protein ID.

    :param identifier: Ensembl Protein ID
    :param retry_in: http code for retrying connections
    :return: Requests response object
    """
    url_root = config.api_ensembl
    url_endpoint = "sequence/id/"
    url = url_root + url_endpoint + identifier
    params = {'type': 'protein'}
    b = Fetcher(url=url, json=True, retry_in=retry_in, **params)
    return b.response


def _fetch_best_structures_pdbe(identifier, retry_in=(429,)):
    """
    Queries the PDBe API sifts mappings best_structures endpoint.

    :param identifier: UniProt ID
    :param retry_in: http code for retrying connections
    :return: Requests response object
    """

    url_root = config.api_pdbe
    url_endpoint = "mappings/best_structures/"
    url = url_root + url_endpoint + identifier
    b = Fetcher(url=url, json=True, retry_in=retry_in)
    return b.response


def _get_preferred_assembly_id(identifier):
    """
    Gets the preferred assembly id for the given PDB ID, from the PDBe API.

    :param identifier: PDB ID
    :return: (str)
    """

    # getting the preferred biological assembly from the PDBe API
    pref_assembly = "1"
    try:
        data = _fetch_summary_properties_pdbe(identifier)
    except Exception as e:
        log.debug("Something went wrong for %s... %s", identifier, e)
    try:
        if data is not None:
            data = data.json()
            nassemblies = data[identifier][0]["assemblies"]
            if len(nassemblies) > 1:
                for entry in nassemblies:
                    if entry["preferred"]:
                        pref_assembly = entry["assembly_id"]
                        break
            else:
                pref_assembly = data[identifier][0]["assemblies"][0]["assembly_id"]
    except Exception as e:
        pref_assembly = "1"
        log.debug("Something went wrong for %s... %s", identifier, e)

    bio_best = str(pref_assembly)
    return bio_best


def _get_ensembl_protein_id_from_mapping(data):
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


def _get_uniprot_id_from_mapping(data, full_entry=False, uniprot_id=None):
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


def _get_preferred_uniprot_id_from_mapping(data):
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


def _get_preferred_ensembl_id_from_mapping(identifiers, cached=False,
                                           uniprot_id=None):
    """
    Takes a list of Ensembl xrefs/ids mapped from a UniProt ID
    and gets the preferred entry (Many-to-one), based on seq
    identity and coverage.

    :param identifiers: list of Ensembl IDs
    :param cached: (boolean) if True, stores a pickle file locally
    :param uniprot_id: if not None means that we want the data for a specific
        UniProt ID
    :return: (str) preferred Ensembl ID
    """

    best_match = None
    curr_ix = -1
    prev_identity = 0
    prev_coverage = 0
    for ix, ensp in enumerate(identifiers):
        info = _fetch_ensembl_ensembl_uniprot_mapping(ensp).json()
        # gets the mapping for a specific uniprot
        data = _get_uniprot_id_from_mapping(info, full_entry=True,
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


def _get_ensembl_species_from_uniprot(data):
    """
    Gets a Species Name from a UniProt organism lookup.

    :param data: Request object from the UniProt Query endpoint.
    :return: (str) formatted species name
    """
    organism = str(data.content, encoding='utf-8').split('\n')[1]
    species = '_'.join(organism.split()[0:2]).lower()
    return species
