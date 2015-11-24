#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Created on 03/06/2015
Functions that handle the reading data files and extracting their information
as a pandas.DataFrame. Also include wrapper functions that select and index
the information. Prefers the use o the wrapper instead the private functions
for better error handling. Both levels are convered by test cases.
"""

import logging
from StringIO import StringIO
from os import path
import sys
from urlparse import parse_qs
from lxml import etree
import pandas as pd
import requests
from config import defaults
from fetcher import fetch_files
from utils import is_valid
from utils import compare_uniprot_ensembl_sequence
from utils import get_url_or_retry

log = logging.getLogger(__name__)


def to_unique(series):
    """Lambda-like expression for returning unique elements of a Series.
    :param series: pandas.Series
    :return: pandas.Series
    """""
    return series.unique()


def _dssp(filename):
    """Parses DSSP file output to a pandas DataFrame.

    :param filename: input SIFTS xml file
    :return: pandas table dataframe
    """
    # column width descriptors
    cols_widths = ((0, 5), (6, 11), (11, 12), (13, 14), (16, 17), (35, 38),
                   (103, 109), (109, 115))
    # simplified headers for the table
    dssp_header = ("dssp_index", "icode", "chain_id", "aa", "ss", "acc", "phi",
                   "psi")
    dssp_table = pd.read_fwf(filename, skiprows=28, names=dssp_header,
                             colspecs=cols_widths, index_col=0,
                             compression=None)
    if dssp_table.empty:
        log.error('DSSP file {} resulted in a empty Dataframe'.format(filename))
        raise ValueError('DSSP file {} resulted in a empty Dataframe'.format(
            filename))
    return dssp_table


def _mmcif_atom(filename, delimiter=None):
    """Parse mmCIF ATOM and HETEROATOM lines to a pandas DataFrame.

    :param filename: input CIF file path
    :return: pandas table dataframe
    """

    _header_mmcif = []
    lines = []
    with open(filename) as inlines:
        for line in inlines:
            if line.startswith("_atom_site."):
                _header_mmcif.append(line.split('.')[1].rstrip())
            elif line.startswith("ATOM"):
                lines.append(line)
            elif line.startswith("HETATM"):
                lines.append(line)
    lines = "".join(lines)

    if delimiter is None:
        return pd.read_table(StringIO(lines),
                             delim_whitespace=True,
                             low_memory=False,
                             names=_header_mmcif,
                             compression=None)
    else:
        return pd.read_table(StringIO(lines),
                             sep=str(delimiter),
                             low_memory=False,
                             names=_header_mmcif,
                             compression=None)


def _sifts_residues(filename, cols=None):
    """Parses the residue fields of a SIFTS XML file to a pandas DataFrame

    :param filename: input SIFTS xml file
    :return: pandas table dataframe
    """

    tree = etree.parse(filename)
    root = tree.getroot()
    namespace = root.nsmap[None]
    nsmap = {'ns': namespace}
    cross_reference = "{{{}}}crossRefDb".format(namespace)
    residue_detail = "{{{}}}residueDetail".format(namespace)
    rows = []
    reference = root.attrib['dbCoordSys']

    for segment in root.iterfind('.//ns:entity[@type="protein"]',
                                 namespaces=nsmap):
        for list_residue in segment.iterfind('.//ns:listResidue',
                                             namespaces=nsmap):
            for residue in list_residue:
                # get residue annotations
                residue_annotation = {}
                # key, value pairs
                for k, v in residue.attrib.items():
                    # skipping dbSource
                    if k == 'dbSource':
                        continue
                    # renaming all keys with dbSource prefix
                    try:
                        k = "{}_{}".format(residue.attrib["dbSource"], k)
                    except KeyError:
                        k = "{}_{}".format("REF", k)
                    # adding to the dictionary
                    residue_annotation[k] = v
                # parse extra annotations for each residue
                for annotation in residue:
                    for k, v in annotation.attrib.items():
                        # crossRefDb entries
                        if annotation.tag == cross_reference:
                            # skipping dbSource
                            if k == 'dbSource':
                                continue
                            # renaming all keys with dbSource prefix
                            try:
                                k = "{}_{}".format(
                                    annotation.attrib["dbSource"], k)
                            except KeyError:
                                k = "{}_{}".format("REF", k)

                        # residueDetail entries
                        elif annotation.tag == residue_detail:
                            # joining dbSource and property keys
                            k = "_".join([annotation.attrib["dbSource"],
                                          annotation.attrib["property"]])

                            # value is the text field in the XML
                            v = annotation.text

                        # adding to the dictionary
                        try:
                            if v in residue_annotation[k]:
                                continue
                            residue_annotation[k].append(v)
                        except KeyError:
                            residue_annotation[k] = v
                        except AttributeError:
                            residue_annotation[k] = [residue_annotation[k]]
                            residue_annotation[k].append(v)

                rows.append(residue_annotation)
    if cols:
        data = pd.DataFrame(rows, columns=cols)
    else:
        data = pd.DataFrame(rows)
    data.columns = [
        col if not col.startswith("PDBe")  # this should come from reference
        else col.replace(reference, "REF")
        for col in data.columns]
    return data


def _sifts_regions(filename):
    """Parse ther region field of the SIFTS XML file to a pandas dataframe

    :param filename: input SIFTS xml file path
    :return: pandas table dataframe
    """

    if not path.isfile(filename):
        raise IOError('File {} not found or unavailable.'.format(filename))

    tree = etree.parse(filename)
    root = tree.getroot()
    namespace = 'http://www.ebi.ac.uk/pdbe/docs/sifts/eFamily.xsd'
    namespace_map = {'ns': namespace}
    db_reference = "{{{}}}db".format(namespace)
    db_detail = "{{{}}}dbDetail".format(namespace)
    rows = []
    regions = {}

    for segment in root.find('.//ns:entity[@type="protein"]',
                             namespaces=namespace_map):
        for region in segment.find('.//ns:listMapRegion',
                                   namespaces=namespace_map):
            # get region annotations
            region_annotation = {}

            # parse extra annotations for each region
            for annotation in region:
                for k, v in annotation.attrib.items():
                    # db entries
                    if annotation.tag == db_reference:
                        # skipping dbSource
                        if k == 'dbSource':
                            continue

                        start = region.attrib['start']
                        end = region.attrib['end']
                        coord = annotation.attrib.get('dbCoordSys', '')
                        region_annotation['Start'] = start
                        region_annotation['End'] = end

                        # region id
                        r = (start, end, coord)

                        # renaming all keys with dbSource prefix
                        k = "{}_{}".format(annotation.attrib["dbSource"], k)

                    # dbDetail entries
                    elif annotation.tag == db_detail:
                        # joining dbSource and property keys
                        k = "_".join([annotation.attrib["dbSource"],
                                      annotation.attrib["property"]])
                        # value is the text field in the XML
                        v = annotation.text

                    # adding to the dictionary
                    try:
                        if v in region_annotation[k]:
                            continue
                        region_annotation[k].append(v)
                    except KeyError:
                        region_annotation[k] = v
                    except AttributeError:
                        region_annotation[k] = [region_annotation[k]]
                        region_annotation[k].append(v)

                    if r not in regions:
                        regions[r] = [region_annotation]
                    else:
                        regions[r].append(region_annotation)

        # group regions together
        for region in regions:
            region_annotation = {}
            for region_annot in regions[region]:
                for k in region_annot:
                    v = region_annot[k]
                    try:
                        if v in region_annotation[k]:
                            continue
                        region_annotation[k].append(v)
                    except KeyError:
                        region_annotation[k] = v
                    except AttributeError:
                        region_annotation[k] = [region_annotation[k]]
                        region_annotation[k].append(v)

            rows.append(region_annotation)
    return pd.DataFrame(rows)


def _pdb_uniprot_sifts_mapping(identifier):
    """Queries the PDBe API for SIFTS mapping between PDB - UniProt. One to many
     relationship expected.

    :param identifier: PDB id
    :return: pandas table dataframe
    """

    if not is_valid(identifier, 'pdbe'):
        raise ValueError(
            "{} is not a valid PDB identifier.".format(identifier))

    sifts_endpoint = "mappings/uniprot/"
    url = defaults.api_pdbe + sifts_endpoint + identifier
    information = get_url_or_retry(url, json=True)

    rows = []
    for uniprot in information[identifier]['UniProt']:
        uniprots = {'uniprot_id': uniprot}
        rows.append(uniprots)
    return pd.DataFrame(rows)


def _uniprot_pdb_sifts_mapping(identifier):
    """Queries the PDBe API for SIFTS mapping between UniProt - PDB entries.
    One to many relationship expected.

    :param identifier: UniProt ID
    :return: pandas table dataframe
    """
    sifts_endpoint = "mappings/best_structures/"
    url = defaults.api_pdbe + sifts_endpoint + str(identifier)
    information = get_url_or_retry(url, json=True)

    rows = []
    for entry in information[identifier]:
        rows.append(entry)
    return pd.DataFrame(rows)


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


def _transcript_variants_ensembl(identifier, missense=True, species=None):
    """Queries the Ensembl API for transcript variants (mostly dbSNP)
    based on Ensembl Protein identifiers (e.g. ENSP00000326864).

    :param identifier: Ensembl Protein ID
    :param missense: if True only fetches missense variants
    :return: pandas table dataframe
    """
    #TODO not using organism
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


def _pdb_validation_to_table(filename, global_parameters=False):
    """Parse the PDB's validation validation file to a pandas DataFrame. Private
     method, prefer its higher level wrapper.

    :type global_parameters: bool
    :param filename: path to file
    :return: table with validation information
    :rtype: pandas.DataFrame
    """
    if not path.isfile(filename):
        raise IOError('File {} not found or unavailable.'.format(filename))

    tree = etree.parse(filename)
    root = tree.getroot()
    if global_parameters:
        global_parameters = root.find('Entry').attrib
        log.info(global_parameters)
    rows = []
    header = set()
    for i, elem in enumerate(root.iterfind('ModelledSubgroup')):
        rows.append(dict(elem.attrib))
        header.update(rows[-1].keys())
    for row in rows:
        not_in = {k: None for k in header.difference(row.keys())}
        row.update(not_in)
    df = pd.DataFrame(rows, columns=header)
    return df


def select_cif(pdb_id, models='first', chains=None, lines=('ATOM',)):
    """Produce table read from mmCIF file

    :param pdb_id: PDB identifier
    :param models: protein structure entity
    :param chains: protein structure chain
    :param lines: choice of ATOM, HETEROATOMS or both.
    :return: Table read to be joined
    """

    cif_path = path.join(defaults.db_mmcif, pdb_id + '.cif')

    try:
        cif_table = _mmcif_atom(cif_path)
    except IOError:
        cif_path = fetch_files(pdb_id, sources='cif',
                               directory=defaults.db_mmcif)[0]
        cif_table = _mmcif_atom(cif_path)

    if models:
        if models == 'first':
            models = cif_table.pdbx_PDB_model_num.iloc[0]

        if isinstance(models, int):
            models = [models]
            try:
                cif_table = cif_table[(
                    cif_table.pdbx_PDB_model_num.isin(models))]
            except AttributeError:
                err = 'Structure {} has only one model, which was kept'.format
                log.info(err(pdb_id))

    if chains:
        if isinstance(chains, str):
            chains = [chains]
        cif_table = cif_table[cif_table.auth_asym_id.isin(chains)]
        if cif_table.empty:
            raise TypeError('Structure {} does not contain {} chain'.format(
                pdb_id, ' '.join(chains)))

    if lines:
        if isinstance(lines, str):
            lines = [lines]
        cif_table = cif_table[cif_table.group_PDB.isin(lines)]

    cif_table = cif_table[cif_table.label_atom_id == 'CA']  # TODO: Make this a parameter?

    if not cif_table['auth_seq_id'].duplicated().any():
        return cif_table.set_index(['auth_seq_id'])

    elif len(cif_table.label_alt_id.unique()) > 1:
        idx = cif_table.groupby(['auth_seq_id']).occupancy.idxmax()
        cif_table = cif_table.ix[idx]
        return cif_table.set_index(['auth_seq_id'])

    elif len(cif_table.pdbx_PDB_ins_code.unique()) > 1:
        cif_table.pdbx_PDB_ins_code.replace("?", "", inplace=True)
        cif_table.auth_seq_id = cif_table.auth_seq_id.astype(str) + \
                                cif_table.pdbx_PDB_ins_code
        return cif_table.set_index(['auth_seq_id'])

    elif not cif_table['pdbe_label_seq_id'].duplicated().any():
        return cif_table.set_index(['pdbe_label_seq_id'])

    log.error('Failed to find unique index for {}'.format(cif_path))
    return cif_table.set_index(['auth_seq_id'])


def select_sifts(pdb_id, chains=None):
    """Produce table ready from SIFTS XML file

    :param pdb_id: PDB identifier
    :param chains: Protein structure chain
    :return: table read to be merged
    """

    sift_path = path.join(defaults.db_sifts, pdb_id + '.xml')

    try:
        sift_table = _sifts_residues(sift_path)
    except IOError:
        sift_path = fetch_files(pdb_id, sources='sifts',
                                directory=defaults.db_sifts)[0]
        sift_table = _sifts_residues(sift_path)
    if chains:
        if isinstance(chains, str):
            chains = [chains]
        sift_table = sift_table[sift_table.PDB_dbChainId.isin(chains)]
    return sift_table


def select_dssp(pdb_id, chains=None):
    """Produce table from DSSP file output

    :param pdb_id: PDB identifier
    :param chains: PDB protein chain
    :return: pandas dataframe
    """

    dssp_path = path.join(defaults.db_dssp, pdb_id + '.dssp')
    try:
        dssp_table = _dssp(dssp_path)
    except IOError:
        dssp_path = fetch_files(pdb_id, sources='dssp',
                                directory=defaults.db_dssp)[0]
        dssp_table = _dssp(dssp_path)
    except StopIteration:
        raise IOError('{} is unreadable.'.format(dssp_path))
    if chains:
        if isinstance(chains, str):
            chains = [chains]
        sel = dssp_table.chain_id.isin(chains)
        if not dssp_table[sel].empty:
            dssp_table = dssp_table[sel]
        else:
            raise ValueError('{} structure DSSP file does not have chains {}'
                             ''.format(pdb_id, ' '.join(chains)))
    if chains and dssp_table.icode.duplicated().any():
        log.info('DSSP file for {} has not unique index'.format(pdb_id))

    return dssp_table.set_index(['icode'])


def select_validation(pdb_id, chains=None):
    """Produces table from PDB validation XML file

    :param pdb_id: PDB identifier
    :param chains: PDB protein chain
    :return: pandas dataframe
    """
    val_path = path.join(defaults.db_pdb,
                         pdb_id + defaults.validation_extension)
    try:
        val_table = _pdb_validation_to_table(val_path)
    except IOError:
        val_path = fetch_files(pdb_id, sources='validation',
                               directory=defaults.db_pdb)[0]
        val_table = _pdb_validation_to_table(val_path)

    if chains:
        if isinstance(chains, str):
            chains = [chains]
        val_table = val_table[val_table.chain.isin(chains)]
    if not val_table.empty:
        val_table.columns = ["val_" + name for name in val_table.columns]
        return val_table
    raise ValueError('Parsing {} resulted in a empty table'.format(val_path))


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


def select_uniprot_variants(identifier):
    """Sumarise variants for a protein in the UniProt

    :param identifier: UniProt ID
    :return: table with variants, rows are residues
    """
    # get organism and sequence for the provided ientifier

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


def sifts_best(identifier, first=False):
    """Retrives the best structures from the SIFTS endpoint in the PDBe api

    :param identifier: Uniprot ID
    :param first: gets the first entry
    :return: url content or url content in json data structure.
    """

    sifts_endpoint = "mappings/best_structures/"
    url = defaults.api_pdbe + sifts_endpoint + str(identifier)
    response = get_url_or_retry(url, json=True)
    return response if not first else response[identifier][0]


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


def _rcsb_description(pdb_id, tag, key):

    api = 'http://www.rcsb.org/pdb/rest/'
    endpoint = 'describeMol'
    query = '?structureId=' + pdb_id

    url = api + endpoint + query

    tree = etree.fromstring(get_url_or_retry(url))
    values = []
    for i in tree.iter(tag):
        values.append(i.attrib[key])

    return values


if __name__ == '__main__':
    X = _sifts_residues('tests/SIFTS/2pah.xml')
