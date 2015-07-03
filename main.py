#!/usr/bin/env python
# -*- coding: utf-8

from os import path
import numpy as np

from to_table import _dssp_to_table, _sifts_residues_to_table, _mmcif_atom_to_table
from utils import _fetch_sifts_best


def merge_tables(uniprot_id, pdb_id=None, chain=None, groupby='CA',
                 defaults=None, cols=None):
    """
    Merges the output from multiple to table method.
    if no pdb_id uses sifts_best_structure
    if no chain uses first
    or if use_all_chains use all chains of the pdb_id
    """
    if cols:
        for col in ('PDB_dbResNum', 'PDB_dbChainId'):
            if col not in cols:
                if not isinstance(cols, list):
                    cols = list(cols)
                cols.append(col)

    def to_list(series):
        return series.tolist()

    def to_unique(series):
        return series.unique()

    if not defaults:
        from config import defaults

    groupby_opts = {'list': {'label_comp_id': 'unique',
                             'label_atom_id': to_list,
                             'label_asym_id': 'unique',
                             'Cartn_x': to_list,
                             'Cartn_y': to_list,
                             'Cartn_z': to_list,
                             'occupancy': 'unique',
                             'B_iso_or_equiv': np.mean},
                    'CA': {'label_comp_id': to_unique,
                           'label_atom_id': to_unique,
                           'label_asym_id': to_unique,
                           'Cartn_x': to_list,
                           'Cartn_y': to_list,
                           'Cartn_z': to_list,
                           'occupancy': to_unique,
                           'B_iso_or_equiv': to_unique}}  # TODO centroid?

    try:
        groupby_opts[groupby]
    except KeyError:
        raise TypeError('Groupby paramenter should be one of {}'.format(
            ', '.join(groupby_opts.keys())))

    if not pdb_id:
        best_pdb = _fetch_sifts_best(uniprot_id, first=True)
        pdb_id = best_pdb['pdb_id']
        chain = best_pdb['chain_id']
    if not chain:
        best_pdb = _fetch_sifts_best(uniprot_id)[pdb_id]
        chain = best_pdb['chain_id']

    cif_path = path.join(defaults.db_mmcif, pdb_id + ".cif")
    dssp_path = path.join(defaults.db_dssp, pdb_id + ".dssp")
    sifts_path = path.join(defaults.db_sifts, pdb_id + ".xml")

    cif_table = _mmcif_atom_to_table(cif_path)
    cif_table = cif_table.query("label_asym_id == @chain & group_PDB == 'ATOM'")
    # next line raises a SettingWithCopyWarning but seems to be correct
    # TODO review next line.
    cif_table.loc[:, "label_seq_id"] = cif_table.loc[:, "label_seq_id"].astype(np.int)

    # TODO keep heteroatoms?

    dssp_table = _dssp_to_table(dssp_path)
    dssp_table = dssp_table.query('chain_id == @chain')
    dssp_table.set_index(['icode'], inplace=True)

    sifts_table = _sifts_residues_to_table(sifts_path, cols=cols)
    sifts_table = sifts_table.query('PDB_dbChainId == @chain')
    sifts_table.set_index(['PDB_dbResNum'], inplace=True)

    if groupby == 'CA':
        cif_table = cif_table.query("label_atom_id == 'CA'")
    cif_table = cif_table.groupby('auth_seq_id').agg(groupby_opts[groupby])

    # This will remove all atoms annotated in sifts but not in mmcif seqres
    # To keep those, sift_table should be the left in join
    cif_sifts = cif_table.join(sifts_table)
    cif_sifts_dssp = cif_sifts.join(dssp_table)
    return cif_sifts_dssp