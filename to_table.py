#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Created on 03/06/2015

"""

import json
import logging
from StringIO import StringIO
from os import path

from lxml import etree
import pandas as pd

from config import defaults
from utils import isvalid_ensembl, request_info_url, get_url_or_retry

log = logging.getLogger(__name__)


def _dssp_to_table(filename):
    """
    Loads and parses DSSP files generating a pandas dataframe.

    :param filename: input SIFTS xml file
    :return: pandas table dataframe
    """
    # column width descriptors
    cols_widths = ((0, 5), (6, 10), (11, 12), (13, 14), (16, 17), (35, 38),
                   (103, 109), (109, 115))
    # simplified headers for the table
    dssp_header = ("dssp_index", "icode", "chain_id", "aa", "ss", "acc", "phi",
                   "psi")
    try:
        return pd.read_fwf(filename, skiprows=28, names=dssp_header,
                           colspecs=cols_widths, index_col=0, compression=None)
    except IOError:
        raise IOError('File {} not found or unavailable.'.format(filename))


def _mmcif_atom_to_table(filename, delimiter=None):
    """
    Loader of mmCIF ATOM and HETEROATOM lines with pandas.

    :param filename: input CIF file
    :return: pandas table dataframe
    """

    if not path.isfile(filename):
        raise IOError('File {} not found or unavailable.'.format(filename))

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
        return pd.read_table(StringIO(lines), delim_whitespace=True,
                             names=_header_mmcif, compression=None)
    else:
        return pd.read_table(StringIO(lines), sep=str(delimiter),
                             names=_header_mmcif, compression=None)


def _sifts_residues_to_table(filename, cols=None):
    """
    Loads and parses SIFTS XML files generating a pandas dataframe.
    Parses the Residue entries.

    :param filename: input SIFTS xml file
    :return: pandas table dataframe
    """

    if not path.isfile(filename):
        raise IOError('File {} not found or unavailable.'.format(filename))

    tree = etree.parse(filename)
    root = tree.getroot()
    namespace = 'http://www.ebi.ac.uk/pdbe/docs/sifts/eFamily.xsd'
    namespace_map = {'ns': namespace}
    cross_reference = "{{{}}}crossRefDb".format(namespace)
    residue_detail = "{{{}}}residueDetail".format(namespace)
    rows = []

    for segment in root.iterfind('.//ns:entity[@type="protein"]', namespaces=namespace_map):
        for list_residue in segment.iterfind('.//ns:listResidue', namespaces=namespace_map):
            for residue in list_residue:
                # get residue annotations
                residue_annotation = {}
                # key, value pairs
                for k, v in residue.attrib.items():
                    # skipping dbSource
                    if k == 'dbSource':
                        continue
                    # renaming all keys with dbSource prefix
                    k = "{}_{}".format(residue.attrib["dbSource"], k)
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
                            k = "{}_{}".format(annotation.attrib["dbSource"], k)

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
    return data


def _sifts_regions_to_table(filename):
    """
    Loads and parses SIFTS XML files generating a pandas dataframe.
    Parses the Regions entries.

    :param filename: input SIFTS xml file
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


def _pdb_uniprot_sifts_mapping_to_table(identifier, verbose=False):
    """
    Queries the PDBe API for SIFTS mapping between PDB - UniProt.
    One to many relationship expected.

    :param identifier: PDB id
    :param verbose: boolean
    :return: pandas table dataframe
    """
    sifts_endpoint = "mappings/uniprot/"
    request = request_info_url(defaults.api_pdbe + sifts_endpoint + identifier,
                               verbose=verbose)
    information = json.loads(request.text)

    rows = []
    for uniprot in information[identifier]['UniProt']:
        uniprots = {'uniprot_id': uniprot}
        rows.append(uniprots)

    return pd.DataFrame(rows)


def _uniprot_pdb_sifts_mapping_to_table(identifier, verbose=False):
    """
    Queries the PDBe API for SIFTS mapping between UniProt - PDB entries.
    One to many relationship expected.

    :param identifier: UniProt ID
    :param verbose: boolean
    :return: pandas table dataframe
    """
    sifts_endpoint = "mappings/best_structures/"
    request = request_info_url(defaults.api_pdbe + sifts_endpoint + identifier,
                               verbose=verbose)
    information = json.loads(request.text)

    rows = []
    for entry in information[identifier]:
        rows.append(entry)
    return pd.DataFrame(rows)


def _uniprot_info_to_table(identifier, retry_in=(503, 500), cols=None):
    """
    Retrive uniprot information from the database.

    :param identifier: UniProt accession identifier
    :return: pandas table dataframe
    """

    if not cols:
        cols = ('entry name', 'reviewed', 'protein names', 'genes', 'organism',
                'sequence', 'length')
    elif isinstance(cols, str):
        cols = ('entry name', cols)

    params = {'query': 'accession:' + identifier,
              'columns': ",".join(cols),
              'format': 'tab',
              'contact': ""}
    url = "http://www.uniprot.org/uniprot/"
    response = get_url_or_retry(url=url, retry_in=retry_in, **params)
    try:
        data = pd.read_table(StringIO(response))
    except ValueError as e:
        log.errore(e)
        data = response
    return data


def _uniprot_ensembl_mapping_to_table(identifier, verbose=False):
    """
    Uses the UniProt mapping service to try and get Ensembl IDs for
    the UniProt accession identifier provided.

    :param identifier: UniProt accession identifier
    :param verbose: boolean
    :return: pandas table dataframe
    """

    information = {}
    rows = []

    # TODO: fix this assuming human variation
    ensembl_endpoint = 'xrefs/symbol/human/'
    params = {'content-type': 'application/json'}
    request = request_info_url("{}{}{}".format(defaults.api_ensembl, ensembl_endpoint,
                                               str(identifier)),
                               params=params,
                               verbose=verbose)

    data = json.loads(request.text)
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


def _transcript_variants_ensembl_to_table(identifier, verbose=False):
    """
    Queries the Ensembl API for transcript variants (mostly dbSNP)
    based on Ensembl Protein identifiers (e.g. ENSP00000326864).

    :param identifier: Ensembl Protein ID
    :param verbose: boolean
    :return: pandas table dataframe
    """
    #TODO delete when safe
    information = {}
    rows = []

    if not isvalid_ensembl(identifier):
        raise ValueError("{} is not a valid Ensembl Accession.".format(identifier))

    ensembl_endpoint = 'overlap/translation/'
    params = {'feature': 'transcript_variation',
              'content-type': 'application/json'}
    request = request_info_url("{}{}{}".format(defaults.api_ensembl,
                                               ensembl_endpoint,
                                               str(identifier)),
                               params=params,
                               verbose=verbose)
    rows = json.loads(request.text)
    return pd.DataFrame(rows)


def _somatic_variants_ensembl_to_table(identifier, verbose=False):
    """
    Queries the Ensembl API for somatic transcript variants (COSMIC)
    based on Ensembl Protein identifiers (e.g. ENSP00000326864).

    :param identifier: Ensembl Protein ID
    :param verbose: boolean
    :return: pandas table dataframe
    """

    if not isvalid_ensembl(identifier):
        raise ValueError("{} is not a valid Ensembl Accession.".format(identifier))

    ensembl_endpoint = 'overlap/translation/'
    params = {'feature': 'somatic_transcript_variation',
              'content-type': 'application/json'}
    request = request_info_url("{}{}{}".format(defaults.api_ensembl,
                                               ensembl_endpoint,
                                               str(identifier)),
                               params=params,
                               verbose=verbose)
    rows = json.loads(request.text)
    return pd.DataFrame(rows)


def _ensembl_variant_to_table(identifier, verbose=False):
    """
    Queries the Ensembl API for variant IDs (e.g rs376845802 or COSM302853).

    :param identifier: variant ID
    :param verbose: boolean
    :return: pandas table dataframe
    """

    if not isvalid_ensembl(identifier, variant=True):
        raise ValueError("{} is not a valid Variation Accession.".format(identifier))

    # TODO: fix this assuming human variation
    ensembl_endpoint = 'variation/human/'
    params = {'content-type': 'application/json'}
    # other params are {'pops': '1', 'phenotypes': '1', 'genotypes': '1'}
    request = request_info_url("{}{}{}".format(defaults.api_ensembl,
                                               ensembl_endpoint,
                                               str(identifier)),
                               params=params,
                               verbose=verbose)
    data = json.loads(request.text)

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

if __name__ == '__main__':
    X = _sifts_residues_to_table("tests/SIFTS/3edv.xml")
    pass
