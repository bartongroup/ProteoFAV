# -*- coding: utf-8

import logging
import pandas as pd
from io import BytesIO
from io import StringIO
from urllib.parse import parse_qs
from pandas.io.json import json_normalize

from proteofav.utils import fetch_from_url_or_retry
from proteofav.utils import row_selector
from proteofav.library import valid_ensembl_species
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


def fetch_uniprot_id_from_name(identifier, retry_in=(429,)):
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


def fetch_uniprot_species_from_id(identifier, retry_in=(429,)):
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


def fetch_uniprot_variants_ebi(identifier, retry_in=(429,)):
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


def fetch_ensembl_uniprot_ensembl_mapping(identifier, retry_in=(429,),
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


def fetch_ensembl_ensembl_uniprot_mapping(identifier, retry_in=(429,)):
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


def fetch_ensembl_transcript_variants(identifier, retry_in=(429,)):
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


def fetch_ensembl_somatic_variants(identifier, retry_in=(429,)):
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


def fetch_ensembl_variants_by_id(identifier, retry_in=(429,),
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


def fetch_ensembl_sequence_from_id(identifier, retry_in=(429,)):
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


def fetch_best_structures_pdbe(identifier, retry_in=(429,)):
    """
    Queries the PDBe API SIFTS mappings best_structures endpoint.

    :param identifier: UniProt ID
    :param retry_in: http code for retrying connections
    :return: Requests response object
    """

    url_root = config.api_pdbe
    url_endpoint = "mappings/best_structures/"
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


def get_preferred_ensembl_id_from_mapping(identifiers, cached=False,
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
        info = fetch_ensembl_ensembl_uniprot_mapping(ensp).json()
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


def get_ensembl_species_from_uniprot(data):
    """
    Gets a Species Name from a UniProt organism lookup.

    :param data: Request object from the UniProt Query endpoint.
    :return: (str) formatted species name
    """
    organism = str(data.content, encoding='utf-8').split('\n')[1]
    species = '_'.join(organism.split()[0:2]).lower()
    return species


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
    url = config.api_uniprot
    response = fetch_from_url_or_retry(url=url, retry_in=retry_in, **params)
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


def _fetch_uniprot_gff(uniprot_id):
    """
    Retrieve Uniprot data from the GFF file.

    :param str uniprot_id: Uniprot accession
    :return pandas.DataFrame: table
    :raises requests.HTTPError: when hits a bad status_code
    """
    url = config.api_uniprot + uniprot_id + ".gff"
    cols = "NAME SOURCE TYPE START END SCORE STRAND FRAME GROUP empty".split()

    data = pd.read_table(url, skiprows=2, names=cols)
    groups = data.GROUP.apply(parse_qs)
    groups = pd.DataFrame.from_records(groups)

    return data.merge(groups, left_index=True, right_index=True)


def map_gff_features_to_sequence(uniprot_id,
                                 query_type='',
                                 group_residues=True,
                                 drop_types=('Helix', 'Beta strand', 'Turn', 'Chain')):
    """
    Map Uniprot GFF features to the protein sequence.

    :param str uniprot_id: Uniprot accession
    :param str query_type: Select type of feature
    :param bool group_residues: by default each row in the resulting table,
        maps to a residue. When set to False, each row represent a feature
        per residue.

    :param tuple drop_types: Filter out some of the features, important to
        remove fetures that spam.

    :return pd.DataFrame: table. Columns will depend on paramenters.
    """

    def annotation_writer(gff_row):
        """
        Establish a set of rules to annotate Uniprot GFF.

        :param pd.Series gff_row: each line in the GFF file.
        :return str: template filled with type-specific fields.
        """
        if not gff_row.ID and not gff_row.Note:
            return gff_row.TYPE
        elif not gff_row.ID:

            return '{0.TYPE}: {0.Note}'.format(gff_row)
        elif not gff_row.Note:

            return '{0.TYPE} ({0.ID})'.format(gff_row)
        else:

            return '{0.TYPE}: {0.Note} ({0.ID})'.format(gff_row)

    table = _fetch_uniprot_gff(uniprot_id)
    if query_type:
        table = table[table.TYPE == query_type]
    elif drop_types:
        table = table[~table.TYPE.isin(drop_types)]

    lines = []
    for i, row in table.iterrows():
        lines.extend({'idx': i, 'annotation': annotation_writer(row)}
                     for i in range(row.START, row.END + 1))
    table = pd.DataFrame(lines)

    if table.empty:
        return table

    if group_residues:

        return table.groupby('idx').agg({'annotation': lambda x: ', '.join(x)})
    else:

        return table


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
        config.api_ensembl, species, uniprot_id)

    return pd.read_json(url)


def _fetch_icgc_variants(identifier):
    """
    Queries the ICGC data portal for the PanCancer variants based on Ensembl
    transcript identifier.
    :param str identifier: Ensembl transcript
    :return pandas.DataFrame: pandas table dataframe
    """
    transition_regex = '(?P<ref>[A-Z])(?P<position>[0-9]+)(?P<new>[A-Z\*])?'
    variation_endpoint = "protein/"
    url = config.api_icgc + variation_endpoint + identifier
    # fetches a nested json
    # normalise the data, making it flat
    data = fetch_from_url_or_retry(url, json=True)
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
    url = config.api_ebi_uniprot + endpoint + uniprot_idd

    data = fetch_from_url_or_retry(url, json=True)
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
    url = config.api_ensembl + ensembl_endpoint + ensembl_ptn_id

    rows = fetch_from_url_or_retry(url, json=True, **params)
    return pd.DataFrame(rows)


def _fetch_sequence_from_ensembl_protein(identifier, protein=True):
    """
    Gets the sequence for an Ensembl identifier.

    :param identifier: Ensembl ID
    :return: sequence
    """
    ensembl_endpoint = "sequence/id/"
    url = config.api_ensembl + ensembl_endpoint + str(identifier)
    header = {'content-type': 'text/plain'}
    if protein:
        params = {'type': 'protein'}
    else:
        params = {}

    return fetch_from_url_or_retry(url, header=header, **params)


def _fetch_uniprot_ensembl_mapping(identifier, species=None):
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
    url = config.api_ensembl + ensembl_endpoint + str(identifier)
    data = fetch_from_url_or_retry(url, json=True)
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