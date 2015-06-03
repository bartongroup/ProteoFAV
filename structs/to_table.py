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

    :param lines:
    :return:
    """
    cols_widths = ((0, 5), (6, 10), (11, 12), (13, 14), (16, 17), (35, 38),
                   (103, 109), (109, 115))
    dssp_header = ("dssp_index", "icode", "chain_id", "aa", "ss", "acc", "phi",
                   "psi")
    return pd.read_fwf(filename, skiprows=28, names=dssp_header,
                       colspecs=cols_widths, index_col=0, compression=None)


def _mmcif_to_table(lines):
    """

    :param lines:
    :return:
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
    lines = [line]
    for line in lines:
        if line.startswith("HETATM"):
            break
        lines.append(line)
    lines.append(line)
    for line in lines:
        if not line.startswith("HETATM"):
            break
        lines.append(line)

    lines = "".join(lines)
    return pd.read_table(StringIO(lines), sep=" ", names=_header_mmcif,
                         compression=None)


def _sifts_to_table(filename):
    """

    :param filename:
    :return:
    """
    tree = etree.parse(filename)
    root = tree.getroot()
    ns = 'http://www.ebi.ac.uk/pdbe/docs/sifts/eFamily.xsd'
    ns_map = {'ns': ns}
    residue_detail = "{{{}}}residueDetail".format(ns)
    rows = []

    for segment in root.find('.//ns:entity[@type="protein"]', namespaces=ns_map):
        # seg_id = segment.attrib["segId"]
        for residue in segment.find('.//ns:listResidue', namespaces=ns_map):

            source = residue.attrib["dbSource"]
            residue_annot = {k.replace("db", "_" + source): v
                             for k, v in residue.attrib.items()
                             if k not in ["dbSource"]}

            for annotation in residue:
                if annotation.tag == residue_detail:
                    new_k = "_".join([annotation.attrib["dbSource"],
                                      annotation.attrib["property"]])
                    try:
                        residue_annot[new_k].append(annotation.text)
                    except KeyError:
                        residue_annot[new_k] = annotation.text
                    except AttributeError:
                        residue_annot[new_k] = [residue_annot[new_k]]
                        residue_annot[new_k].append(annotation.text)

                else:
                    source = annotation.attrib["dbSource"]
                    for k, v in annotation.attrib.items():
                        if k == "dbSource":
                            continue
                        new_k = k.replace("db", "{}_".format(source))
                        try:
                            if v in residue_annot[new_k]:
                                continue
                            residue_annot[new_k].append(v)
                        except KeyError:
                            residue_annot[new_k] = v
                        except AttributeError:
                            residue_annot[new_k] = [residue_annot[new_k]]
                            residue_annot[new_k].append(v)
                rows.append(residue_annot)
    return pd.DataFrame(rows)

