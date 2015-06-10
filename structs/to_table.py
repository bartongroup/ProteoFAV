#!/usr/bin/env python

"""
Created on 18:16 03/06/15 2015

"""
import logging
from StringIO import StringIO
from lxml import etree

import pandas as pd

logger = logging.getLogger(__name__)



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
    return pd.read_fwf(filename, skiprows=28, names=dssp_header,
                       colspecs=cols_widths, index_col=0, compression=None)


def _mmcif_to_table(lines):
    """
    Testing a loader of mmCIF ATOM lines with pandas.

    :param lines: mmCIF ATOM lines
    :return: pandas table dataframe
    """
    _header_mmcif = (
    "group_PDB", "id", "type_symbol", "label_atom_id", "label_alt_id",
    "label_comp_id", "label_asym_id", "label_entity_id", "label_seq_id",
    "pdbx_PDB_ins_code", "Cartn_x", "Cartn_y", "Cartn_z", "occupancy",
    "B_iso_or_equiv", "Cartn_x_esd", "Cartn_y_esd", "Cartn_z_esd",
    "occupancy_esd", "B_iso_or_equiv_esd", "pdbx_formal_charge",
    "auth_seq_id", "auth_comp_id", "auth_asym_id", "auth_atom_id",
    "pdbx_PDB_model_num", "pdbe_label_seq_id")

    for line in lines:
        if line.startswith("ATOM"):
            break
    rows = [line]
    for line in lines:
        if line.startswith("HETATM"):
            break
        rows.append(line)
    rows.append(line)
    for line in lines:
        if not line.startswith("HETATM"):
            break
        rows.append(line)

    rows = "".join(rows)
    return pd.read_table(StringIO(rows), sep=" ", names=_header_mmcif,
                         compression=None)


def _sifts_to_table_residues(filename):
    """
    Loads and parses SIFTS XML files generating a pandas dataframe.
    Parses the Residue entries.

    TODO: Test that the file exists and it is valid?

    :param filename: input SIFTS xml file
    :return: pandas table dataframe
    """

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


def _sifts_to_table_regions(filename):
    """
    Loads and parses SIFTS XML files generating a pandas dataframe.
    Parses the Regions entries.

    TODO: Test that the file exists and it is valid?

    :param filename: input SIFTS xml file
    :return: pandas table dataframe
    """

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
                        # skiping dbSource
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


if __name__ == '__main__':
    pass