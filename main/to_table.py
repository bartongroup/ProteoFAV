#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Created on 18:16 03/06/15 2015

"""
import sys
sys.path.insert(0, '../')
import os
import logging
from StringIO import StringIO
from lxml import etree

import pandas as pd

from utils.config import get_config
from utils.utils import request_info_url

logger = logging.getLogger(__name__)



def _dssp_to_table(filename):
    """
    Loads and parses DSSP files generating a pandas dataframe.

    :param filename: input SIFTS xml file
    :return: pandas table dataframe
    """

    if not os.path.isfile(filename):
        raise IOError('File {} not found or unavailable.'.format(filename))

    # column width descriptors
    cols_widths = ((0, 5), (6, 10), (11, 12), (13, 14), (16, 17), (35, 38),
                   (103, 109), (109, 115))
    # simplified headers for the table
    dssp_header = ("dssp_index", "icode", "chain_id", "aa", "ss", "acc", "phi",
                   "psi")
    return pd.read_fwf(filename, skiprows=28, names=dssp_header,
                       colspecs=cols_widths, index_col=0, compression=None)


def _mmcif_atom_to_table(filename, delimiter=None):
    """
    Testing a loader of mmCIF ATOM lines with pandas.

    :param filename: input CIF file
    :return: pandas table dataframe
    """

    if not os.path.isfile(filename):
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


def _uniprot_info_to_table(identifier, verbose=False):
    """
    Fetches some information including the sequence of a particular
    UniProt entry.

    :param identifier: UniProt accession identifier
    :param verbose: boolean
    :return: pandas table dataframe
    """

    information = {}
    rows = []

    params = {'query': 'accession:%s' % identifier,
              'columns': 'entry name,reviewed,protein names,genes,organism,sequence,length',
              'format': 'tab',
              'contact': 'fmmarquesmadeira@dundee.ac.uk'}
    config = get_config('http_uniprot')
    request = request_info_url(config.http_uniprot, params, verbose=verbose)

    data = request.text.split('\n')
    for i, line in enumerate(data):
        if i == 1 and line != '':
            line = line.split('\t')
            information['Name'] = line[0]
            information['Status'] = line[1]
            information['Protein'] = line[2]
            information['Genes'] = line[3]
            information['Organism'] = line[4]
            information['Sequence'] = line[5]
            information['Length'] = int(line[6])

    rows.append(information)
    return pd.DataFrame(rows)


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

    ensembl_mappings = ["ENSEMBL_ID", "ENSEMBL_PRO_ID", "ENSEMBL_TRS_ID"]
    for ensembl in ensembl_mappings:
        params = {'from': 'ACC',
                  'to': ensembl,
                  'format': 'tab',
                  'query': identifier,
                  'contact': 'fmmarquesmadeira@dundee.ac.uk'}

        config = get_config('http_uniprot_mapping')
        request = request_info_url(config.http_uniprot_mapping, params,
                                   verbose=verbose)

        data = request.text.split('\n')
        for i, line in enumerate(data):
            if i >= 1 and line != '':
                line = line.split('\t')
                try:
                    if line[1] in information[ensembl]:
                        continue
                    information[ensembl].append(line[1])
                except KeyError:
                    information[ensembl] = line[1]
                except AttributeError:
                    information[ensembl] = [information[ensembl]]
                    information[ensembl].append(line[1])

    rows.append(information)
    return pd.DataFrame(rows)


if __name__ == '__main__':
    # testing routines
    pass
    info = _uniprot_info_to_table('O96013')
    print(info)

    print(info.columns.values)

    info = _uniprot_ensembl_mapping_to_table('O96013')
    print(info)
