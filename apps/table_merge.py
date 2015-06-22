#!/usr/bin/env python
# -*- coding: utf-8
import logging
from os import path

import numpy as np

from pdbs.to_table import _mmcif_atom_to_table as mmcif_to_table
from pdbs.to_table import _dssp_to_table as dssp_to_table
from sifts.to_table import _sifts_residues_to_table as sifts_to_table

from utils.config import defaults
from utils.utils import get_url_or_retry

log = logging.getLogger(__name__)


def _fetch_sifts_best(uniprot_id, first=False):
    url = 'http://wwwdev.ebi.ac.uk/pdbe/api/mappings/best_structures/'
    url = url + uniprot_id
    response = get_url_or_retry(url, json=True)
    return response if not first else response[uniprot_id][0]


def merge_tables(uniprot_id, pdb_id=None, chain=None, groupby='CA'):
    """
    Merges the output from multiple to table method.
    if no pdb_id uses sifts_best_structure
    if no chain uses first
    or if use_all_chains use all chains of the pdb_id
    """

    def to_list(series):
        return series.tolist()

    groupby_opts = {'list': {'label_comp_id': 'unique',
                             'label_atom_id': to_list,
                             'label_asym_id': 'unique',
                             'Cartn_x': to_list,
                             'Cartn_y': to_list,
                             'Cartn_z': to_list,
                             'occupancy': 'unique',
                             'B_iso_or_equiv': np.mean},
                    'CA': {'label_comp_id': 'unique',
                           'label_atom_id': 'unique',
                           'label_asym_id': 'unique',
                           'Cartn_x': 'unique',
                           'Cartn_y': 'unique',
                           'Cartn_z': 'unique',
                           'occupancy': 'unique',
                           'B_iso_or_equiv': 'unique'}}  # TODO centroid?

    try:
        groupby_opts[groupby]
    except KeyError:
        raise TypeError('Groupby paramenter should be one of {}'.format(
            ', '.join(groupby_opts.keys())))

    if not pdb_id:
        best_pdb = _fetch_sifts_best(uniprot_id, first=True)
        pdb_id = best_pdb['pdb_id']
        chain = best_pdb['chain_id']
        # start stop
    if not chain:
        best_pdb = _fetch_sifts_best(uniprot_id)[pdb_id]
        chain = best_pdb['chain_id']

    cif_path = path.join(path.dirname(__file__), '../tests/CIF/2w4o.cif')
    # cif_path = path.join(defaults.db_mmcif, pdb_id, ".cif")
    dssp_path = path.join(path.dirname(__file__), '../tests/DSSP/2w4o.dssp')
    # dssp_path = path.join(defaults.db_dssp, pdb_id, ".dssp")
    sifts_path = path.join(path.dirname(__file__), '../tests/SIFTS/2w4o.xml')
    # sifts_path = path.join(defaults.db_sifts, pdb_id, ".xml")

    cif_table = mmcif_to_table(cif_path)
    cif_table = cif_table.query("label_asym_id == @chain & group_PDB == 'ATOM'")
    # TODO keep heteroatoms?

    dssp_table = dssp_to_table(dssp_path)
    dssp_table = dssp_table.query('chain_id == @chain')
    dssp_table.set_index(['icode'], inplace=True)

    sifts_table = sifts_to_table(sifts_path)
    sifts_table = sifts_table.query('PDB_dbChainId == @chain')
    sifts_table.set_index(['PDB_dbResNum'], inplace=True)

    if groupby == 'CA':
        cif_table = cif_table.query("label_atom_id == 'CA'")
    cif_lines = cif_table.groupby('label_seq_id').agg(groupby_opts[groupby])

    # This remove all atoms annotated in sifts but not in mmcif seqres
    # To keep those, sift_table should be the left in join
    cif_sifts = cif_table.join(sifts_table)
    cif_sifts_dssp = cif_sifts.join(dssp_table)
    return cif_sifts_dssp


if __name__ == '__main__':
    merge_tables('Q16566')
