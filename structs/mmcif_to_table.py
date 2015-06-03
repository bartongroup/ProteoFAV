#!/usr/bin/env python
# coding=utf-8
"""
Created on 17:44 03/06/15 2015 

"""
import pandas as pd
from StringIO import StringIO

header_mmcif = (
    "group_PDB", "id", "type_symbol", "label_atom_id", "label_alt_id",
    "label_comp_id", "label_asym_id", "label_entity_id", "label_seq_id",
    "pdbx_PDB_ins_code", "Cartn_x", "Cartn_y", "Cartn_z", "occupancy",
    "B_iso_or_equiv", "Cartn_x_esd", "Cartn_y_esd", "Cartn_z_esd",
    "occupancy_esd", "B_iso_or_equiv_esd", "pdbx_formal_charge",
    "auth_seq_id", "auth_comp_id", "auth_asym_id", "auth_atom_id",
    "pdbx_PDB_model_num", "pdbe_label_seq_id")

def mmcif_to_table(lines):
#
# with open("/Users/tbrittoborges/Downloads/2pah_updated.cif.txt") as open_f:
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
    return pd.read_table(StringIO(lines), sep=" ",
                         names=header_mmcif, compression=None)