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


def fetch_summary_properties_pdbe(identifier, retry_in=(429,)):
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


def get_preferred_assembly_id(identifier):
    """
    Gets the preferred assembly id for the given PDB ID, from the PDBe API.

    :param identifier: PDB ID
    :return: (str)
    """

    # getting the preferred biological assembly from the PDBe API
    pref_assembly = "1"
    try:
        data = fetch_summary_properties_pdbe(identifier)
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


def fetch_uniprot_id_from_name(identifier, cached=False, retry_in=(429,)):
    """
    Retrieve UniProt ID from Name.

    :param identifier: UniProt accession identifier
    :param cached: (boolean) if True, stores a pickle file locally
    :param retry_in: http code for retrying connections
    :return: Requests response object
    """

    url_root = config.http_uniprot
    url_endpoint = "?query={}&columns=id&format=list".format(identifier)
    url = url_root + url_endpoint
    b = Fetcher(url=url, cached=cached,
                cache_output="{}_id.pkl".format(identifier),
                json=False, retry_in=retry_in)
    return b.response


def fetch_uniprot_species_from_id(identifier, cached=False, retry_in=(429,)):
    """
    Retrieve Species from UniProt ID.

    :param identifier: UniProt accession identifier
    :param cached: (boolean) if True, stores a pickle file locally
    :param retry_in: http code for retrying connections
    :return: Requests response object
    """

    url_root = config.http_uniprot
    url_endpoint = "?query={}&columns=organism&format=tab".format(identifier)
    url = url_root + url_endpoint
    b = Fetcher(url=url, cached=cached,
                cache_output="{}_org.pkl".format(identifier),
                json=False, retry_in=retry_in)
    return b.response


def fetch_uniprot_variants_ebi(identifier, cached=False, retry_in=(429,)):
    """
    Queries the EBI Proteins API for variants.
    based on UniProt identifiers (e.g. O15294).

    :param identifier: UniProt ID
    :param cached: (boolean) if True, stores a pickle file locally
    :param retry_in: http code for retrying connections
    :return: Requests response object
    """

    url_root = config.api_proteins
    url_endpoint = "variation/"
    url = url_root + url_endpoint + identifier
    b = Fetcher(url=url, cached=cached,
                cache_output="{}_vars_uni.pkl".format(identifier),
                json=True, retry_in=retry_in)
    return b.response


def fetch_ensembl_uniprot_ensembl_mapping(identifier, cached=False, retry_in=(429,),
                                          species='homo_sapiens'):
    """
    Uses the Ensembl REST mapping service to try and get Ensembl IDs for
    the UniProt accession identifier provided.

    :param identifier: UniProt accession identifier
    :param cached: (boolean) if True, stores a pickle file locally
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

    b = Fetcher(url=url, cached=cached,
                cache_output="{}_uni_ens.pkl".format(identifier),
                json=True, retry_in=retry_in)
    return b.response


def fetch_ensembl_ensembl_uniprot_mapping(identifier, cached=False, retry_in=(429,)):
    """
    Uses the Ensembl REST mapping service to try and get UniProt IDs for
    the Ensembl Protein accession identifier provided.

    :param identifier: Ensembl Protein accession identifier
    :param cached: (boolean) if True, stores a pickle file locally
    :param retry_in: http code for retrying connections
    :return: Requests response object
    """

    url_root = config.api_ensembl
    url_endpoint = "xrefs/id/"
    url = url_root + url_endpoint + identifier
    # params = {'external_db': "Uniprot/SWISSPROT"}
    # params = {'external_db': "Uniprot/SPTREMBL"}

    b = Fetcher(url=url, cached=cached,
                cache_output="{}_ens_uni.pkl".format(identifier),
                json=True, retry_in=retry_in)
    return b.response


def fetch_ensembl_transcript_variants(identifier, cached=False, retry_in=(429,)):
    """
    Queries the Ensembl REST API for transcript variants (mostly from dbSNP).
    based on Ensembl Protein identifiers (e.g. ENSP00000275603).

    :param identifier: Ensembl Protein ID
    :param cached: (boolean) if True, stores a pickle file locally
    :param retry_in: http code for retrying connections
    :return: Requests response object
    """

    url_root = config.api_ensembl
    url_endpoint = "overlap/translation/"
    url = url_root + url_endpoint + identifier
    # params = {'feature': 'transcript_variation', 'type': 'missense_variant'}
    params = {'feature': 'transcript_variation'}
    b = Fetcher(url=url, cached=cached,
                cache_output="{}_vars_ens.pkl".format(identifier),
                json=True, retry_in=retry_in, **params)
    return b.response


def fetch_ensembl_somatic_variants(identifier, cached=False, retry_in=(429,)):
    """
    Queries the Ensembl REST API for somatic variants (mostly from COSMIC).
    based on Ensembl Protein identifiers (e.g. ENSP00000275603).

    :param identifier: Ensembl Protein ID
    :param cached: (boolean) if True, stores a pickle file locally
    :param retry_in: http code for retrying connections
    :return: Requests response object
    """

    url_root = config.api_ensembl
    url_endpoint = "overlap/translation/"
    url = url_root + url_endpoint + identifier
    params = {'feature': 'somatic_transcript_variation'}
    b = Fetcher(url=url, cached=cached,
                cache_output="{}_vars_cos.pkl".format(identifier),
                json=True, retry_in=retry_in, **params)
    return b.response


def fetch_ensembl_variants_by_id(identifier, cached=False, retry_in=(429,),
                                 species='homo_sapiens'):
    """
    Queries the Ensembl API for variant IDs (e.g rs376845802 or COSM302853).

    :param identifier: variant ID
    :param cached: (boolean) if True, stores a pickle file locally
    :param retry_in: http code for retrying connections
    :param species: Ensembl species
    :return: Requests response object
    """

    if isinstance(identifier, str):
        url_root = config.api_ensembl
        url_endpoint = "variation/{}/".format(species)
        url = url_root + url_endpoint + identifier
        # params = {'pops': '1', 'phenotypes': '1'} #, 'genotypes': '1', 'population_genotypes': '1'}
        b = Fetcher(url=url, cached=cached,
                    cache_output="{}_vars_ids.pkl".format(identifier),
                    json=True, retry_in=retry_in)
        return b.response
    elif isinstance(identifier, list):
        data = """{ "ids" : ["%s"]}""" % ('", "'.join(identifier))
        header = {"Accept": "application/json"}
        url_root = config.api_ensembl
        url_endpoint = "variation/{}".format(species)
        url = url_root + url_endpoint
        # params = {'pops': '1', 'phenotypes': '1'} #, 'genotypes': '1'}
        b = Fetcher(url=url, cached=cached,
                    cache_output="{}_vars_ids.pkl".format(identifier),
                    json=True, retry_in=retry_in,
                    post=True, data=data, header=header)
        return b.response


def fetch_ensembl_sequence_from_id(identifier, cached=False, retry_in=(429,)):
    """
    Queries the Ensembl REST API for the sequence of a Ensembl Protein ID.

    :param identifier: Ensembl Protein ID
    :param cached: (boolean) if True, stores a pickle file locally
    :param retry_in: http code for retrying connections
    :return: Requests response object
    """
    url_root = config.api_ensembl
    url_endpoint = "sequence/id/"
    url = url_root + url_endpoint + identifier
    params = {'type': 'protein'}
    b = Fetcher(url=url, cached=cached,
                cache_output="{}_ens_seq.pkl".format(identifier),
                json=True, retry_in=retry_in, **params)
    return b.response


def fetch_best_structures_pdbe(identifier, cached=False, retry_in=(429,)):
    """
    Queries the PDBe API SIFTS mappings best_structures endpoint.

    :param identifier: UniProt ID
    :param cached: (boolean) if True, stores a pickle file locally
    :param retry_in: http code for retrying connections
    :return: Requests response object
    """

    url_root = config.api_pdbe
    url_endpoint = "mappings/best_structures/"
    url = url_root + url_endpoint + identifier
    b = Fetcher(url=url, cached=cached,
                cache_output="{}_bs.pkl".format(identifier),
                json=True, retry_in=retry_in)
    return b.response
