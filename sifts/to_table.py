#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Created on 11/06/2015

"""

import sys
sys.path.insert(0, '../')
import os
import logging
import json
from lxml import etree
import pandas as pd

from utils.config import get_config
from utils.utils import request_info_url
from utils.utils import isvalid_uniprot
from utils.utils import isvalid_pdb

logger = logging.getLogger(__name__)


def _sifts_residues_to_table(filename):
    """
    Loads and parses SIFTS XML files generating a pandas dataframe.
    Parses the Residue entries.

    :param filename: input SIFTS xml file
    :return: pandas table dataframe
    """

    if not os.path.isfile(filename):
        raise IOError('File {} not found or unavailable.'.format(filename))

    tree = etree.parse(filename)
    root = tree.getroot()
    namespace = 'http://www.ebi.ac.uk/pdbe/docs/sifts/eFamily.xsd'
    namespace_map = {'ns': namespace}
    cross_reference = "{{{}}}crossRefDb".format(namespace)
    residue_detail = "{{{}}}residueDetail".format(namespace)
    rows = []

    for segment in root.find('.//ns:entity[@type="protein"]', namespaces=namespace_map):
        for residue in segment.find('.//ns:listResidue', namespaces=namespace_map):
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
    return pd.DataFrame(rows)


def _sifts_regions_to_table(filename):
    """
    Loads and parses SIFTS XML files generating a pandas dataframe.
    Parses the Regions entries.

    :param filename: input SIFTS xml file
    :return: pandas table dataframe
    """

    if not os.path.isfile(filename):
        raise IOError('File {} not found or unavailable.'.format(filename))

    tree = etree.parse(filename)
    root = tree.getroot()
    namespace = 'http://www.ebi.ac.uk/pdbe/docs/sifts/eFamily.xsd'
    namespace_map = {'ns': namespace}
    db_reference = "{{{}}}db".format(namespace)
    db_detail = "{{{}}}dbDetail".format(namespace)
    rows = []
    regions = {}

    for segment in root.find('.//ns:entity[@type="protein"]', namespaces=namespace_map):
        for region in segment.find('.//ns:listMapRegion', namespaces=namespace_map):
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

    if not isvalid_pdb(identifier):
        raise ValueError("{} is not a valid PDB Accession.".format(identifier))

    config = get_config('api_pdbe')
    sifts_endpoint = "mappings/uniprot/"
    request = request_info_url("{}{}{}".format(config.api_pdbe, sifts_endpoint,
                                               str(identifier), verbose=verbose))
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

    if not isvalid_uniprot(identifier):
        raise ValueError("{} is not a valid UniProt Accession.".format(identifier))

    config = get_config('api_pdbe')
    sifts_endpoint = "mappings/best_structures/"
    request = request_info_url("{}{}{}".format(config.api_pdbe, sifts_endpoint,
                                               str(identifier), verbose=verbose))
    information = json.loads(request.text)

    rows = []
    for entry in information[identifier]:
        rows.append(entry)
    return pd.DataFrame(rows)


if __name__ == '__main__':
    # testing routines
    pass
