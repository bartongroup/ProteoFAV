#!/usr/bin/env python
# -*- coding: utf-8
import logging

from os import path
import numpy as np

from to_table import _dssp_to_table, _sifts_residues_to_table, _mmcif_atom_to_table, pd
from utils import _fetch_sifts_best

log = logging.getLogger(__name__)

def merge_tables(uniprot_id=None, pdb_id=None, chain=None, groupby='CA', defaults=None, cols=None,
                 model=1):
    """
    Merges the output from multiple to table method.
    if no pdb_id uses sifts_best_structure
    if no chain uses first
    or if use_all_chains use all chains of the pdb_id
    :param uniprot_id:
    :param pdb_id:
    :param chain:
    :param groupby:
    :param defaults:
    :param cols:
    :param model:
    """

    # TODO: cols

    if not any((uniprot_id, pdb_id)):
        raise TypeError("One of the following arguments is expected:"
                        "uniprot_id or pdb_id")


    def to_list(series):
        return series.tolist()

    def to_unique(series):
        return series.unique()

    def to_first(series):
        return series.first()


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
                           'Cartn_x': to_unique,
                           'Cartn_y': to_unique,
                           'Cartn_z': to_unique,
                           'occupancy': to_unique,
                           'B_iso_or_equiv': to_unique,
                           'id': to_unique},
                    'SC': {'label_comp_id': to_unique,
                           'label_atom_id': to_unique,
                           'label_asym_id': to_unique,
                           'Cartn_x': 'mean',
                           'Cartn_y': 'mean',
                           'Cartn_z': 'mean',
                           'occupancy': 'mean',
                           'B_iso_or_equiv': 'mean',
                           'label_alt_id': to_unique,
                           'id': to_unique},
                    'centroid':{},
                    'all_atoms':{}} # TODO

    try:
        groupby_opts[groupby]
    except KeyError:
        raise TypeError('Groupby paramenter should be one of {}'.format(
            ', '.join(groupby_opts.keys())))

    if not pdb_id:
        best_pdb = _fetch_sifts_best(uniprot_id, first=True)
        pdb_id = best_pdb['pdb_id']
        chain = best_pdb['chain_id']
        log.info("Best structure, chain: {}|{} for {} ".format(pdb_id, chain, uniprot_id))

    if not chain:
        try:
            best_pdb = _fetch_sifts_best(uniprot_id)[pdb_id]
        except KeyError:
            err = "Structure {} not found in best structures".format
            log.error(err(pdb_id))
        chain = best_pdb['chain_id']
        log.info("Best structure, chain: {}|{} for {} ".format(pdb_id, chain, uniprot_id))

    cif_path = path.join(defaults.db_mmcif, pdb_id + '.cif')
    dssp_path = path.join(defaults.db_dssp, pdb_id + '.dssp')
    sifts_path = path.join(defaults.db_sifts, pdb_id + '.xml')

    cif_table = _mmcif_atom_to_table(cif_path)
    try:
        cif_table = cif_table[(cif_table.pdbx_PDB_model_num == model) &
                              (cif_table.label_asym_id == chain) &
                              (cif_table.group_PDB == 'ATOM')]
    except AttributeError:
        err = 'Structure {} has only one model'.format
        log.info(err(pdb_id))
        cif_table = cif_table[(cif_table.label_asym_id == chain) &
                              (cif_table.group_PDB == 'ATOM')]
    # next line raises a SettingWithCopyWarning but seems to be correct
    cif_table.loc[:, 'auth_seq_id'] = cif_table.loc[:, 'auth_seq_id'].astype(np.int)

    if groupby == 'CA':
        cif_table = cif_table[cif_table.label_atom_id == 'CA']
        # assert cif_table.label_seq_id is a sequence with number of aas
    elif groupby == 'SC':
        cif_table = cif_table[~cif_table.label_atom_id.isin(['CA', 'C', 'O', 'N'])]
        # assert cif_table.label_seq_id is a sequence with number of aas

    # Check the existence of alt locations
    if len(cif_table.label_alt_id.unique()) > 1:
        cif_table = cif_table.groupby(['auth_seq_id', 'label_alt_id']).agg(
            groupby_opts[groupby])
        idx = cif_table.groupby(level=0).apply(lambda x: x.occupancy.idxmax())
        cif_table = cif_table.ix[idx]
        cif_table.reset_index(level=1, drop=True, inplace=True)
        log.info('Has atoms in alternative location')
        # TODO what should be done with atoms in alt location?
    else:
        cif_table = cif_table.groupby('auth_seq_id').agg(groupby_opts[groupby])

    dssp_table = _dssp_to_table(dssp_path)
    dssp_table = dssp_table[dssp_table.chain_id == chain]
    dssp_table.loc[:, 'icode'] = dssp_table.loc[:, 'icode'].astype(np.int)
    dssp_table.set_index(['icode'], inplace=True)

    sifts_table = _sifts_residues_to_table(sifts_path)
    sifts_table = sifts_table[sifts_table.PDB_dbChainId == chain]
    sifts_table.loc[:, 'PDB_dbResNum'] = sifts_table.loc[:, 'PDB_dbResNum'].astype(np.int)
    sifts_table.set_index(['PDB_dbResNum'], inplace=True)

    # Sift data in used as left element int the join, because MMCIF atom lines and DSSP
    # Ignore not observed residues.
    sifts_cif = sifts_table.join(cif_table)
    cif_sifts_dssp = sifts_cif.join(dssp_table)
    return cif_sifts_dssp

if __name__ == '__main__':
    from config import Defaults
    defaults = Defaults('config.txt')
    defaults.db_mmcif = 'tests/CIF'
    defaults.db_dssp = 'tests/DSSP'
    defaults.db_sifts = 'tests/SIFTS'

    X = merge_tables(pdb_id='2w4o', chain='A', defaults=defaults)
    pass