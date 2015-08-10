#!/usr/bin/env python
# -*- coding: utf-8
import logging
from os import path

import numpy as np

from to_table import _dssp_to_table, _sifts_residues_to_table, \
    _mmcif_atom_to_table
from utils import _fetch_sifts_best
from library import three_to_single_aa

log = logging.getLogger(__name__)


def merge_tables(uniprot_id=None, pdb_id=None, chain=None, model=1,
                 groupby='CA', default=None, validate=True):
    """
    Merges the output from multiple to table method.
    if no pdb_id uses sifts_best_structure
    if no chain uses first
    or if use_all_chains use all chains of the pdb_id

    :param uniprot_id: UniProt ID
    :param pdb_id: PDB ID
    :param chain: CHAIN ID
    :param model: MODEL ID
    :param groupby: options on how to group residues
    :param default: default parameters
    :param validate: boolean for structure validation entries
    """
    # TODO: cols

    if not any((uniprot_id, pdb_id)):
        raise TypeError("One of the following arguments is expected:"
                        "uniprot_id or pdb_id")

    to_list = lambda series: series.tolist()
    to_unique = lambda series: series.unique()

    if not default:
        from config import defaults as default

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
                    'centroid': {},
                    'all_atoms': {}}  # TODO

    try:
        groupby_opts[groupby]
    except KeyError:
        raise TypeError('Groupby paramenter should be one of {}'.format(
            ', '.join(groupby_opts.keys())))

    if not pdb_id:
        best_pdb = _fetch_sifts_best(uniprot_id, first=True)
        pdb_id = best_pdb['pdb_id']
        chain = best_pdb['chain_id']
        log.info("Best structure, chain: {}|{} for {} ".format(pdb_id, chain,
                                                               uniprot_id))

    if not chain:
        try:
            best_pdb = _fetch_sifts_best(uniprot_id)[pdb_id]
        except KeyError:
            err = "Structure {} not found in best structures".format
            log.error(err(pdb_id))
        chain = best_pdb['chain_id']
        log.info("Best structure, chain: {}|{} for {} ".format(pdb_id, chain,
                                                               uniprot_id))

    cif_path = path.join(default.db_mmcif, pdb_id + '.cif')
    dssp_path = path.join(default.db_dssp, pdb_id + '.dssp')
    sifts_path = path.join(default.db_sifts, pdb_id + '.xml')

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
    if 'pdbe_label_seq_id' in cif_table.columns:
        groupby_opts[groupby].update({'pdbe_label_seq_id': to_unique})
        log.info('Column pdbe_label_seq_id present')

    # TODO here to line 139 should be a function, that handles the cif table
    if groupby == 'CA':
        cif_table = cif_table[cif_table.label_atom_id == 'CA']
    elif groupby == 'SC':
        cif_table = cif_table[
            ~cif_table.label_atom_id.isin(['CA', 'C', 'O', 'N'])]

    # Check the existence of alt locations
    if len(cif_table.label_alt_id.unique()) > 1:
        cif_table = cif_table.groupby(['auth_seq_id', 'label_alt_id']).agg(
            groupby_opts[groupby])
        # We get the atom with max occupancy
        idx = cif_table.groupby(level=0).apply(lambda x: x.occupancy.idxmax())
        cif_table = cif_table.ix[idx]
        cif_table.reset_index(level=1, drop=True, inplace=True)
        log.info('Has atoms in alternative location')
        # TODO what should be done with atoms in alt location?
    else:
        cif_table = cif_table.groupby('auth_seq_id').agg(groupby_opts[groupby])

    dssp_table = _dssp_to_table(dssp_path)
    dssp_table = dssp_table[dssp_table.chain_id == chain]
    dssp_table.set_index(['icode'], inplace=True)
    cif_dssp = cif_table.join(dssp_table)

    if validate:
        # Fill all missing values with gaps
        cif_dssp['dssp_aa'] = cif_dssp.aa.fillna('-')
        cif_dssp['cif_aa'] = cif_dssp.label_comp_id.fillna('-')
        # From three letter to sigle letters or X if not a standard aa
        cif_dssp['cif_aa'] = cif_dssp['cif_aa'].apply(three_to_single_aa.get,
                                                      args='X')
        # Check if the sequences are the same
        if not (cif_dssp.dssp_aa == cif_dssp.cif_aa).all():
            raise ValueError('{pdb_id}|{chain} Cif and DSSP files have diffent'
                             ' sequences.'.format(pdb_id=pdb_id, chain=chain))

    sifts_table = _sifts_residues_to_table(sifts_path)
    sifts_table = sifts_table[sifts_table.PDB_dbChainId == chain]
    try:
        sifts_table.loc[:, 'PDB_dbResNum'] = sifts_table.loc[:, 'PDB_dbResNum'].astype(np.int)
        sifts_table.set_index(['PDB_dbResNum'], inplace=True)
    except ValueError: # Means it has alpha numeric insertion code, use something else as index
        sifts_table.loc[:, 'REF_dbResNum'] = sifts_table.loc[:, 'REF_dbResNum'].astype(np.int)
        sifts_table.set_index(['REF_dbResNum'], inplace=True)
        cif_dssp.set_index(['pdbe_label_seq_id'], inplace=True)

    # Sift data in used as base to keep not observed residues info.
    sifts_cif_dssp = sifts_table.join(cif_dssp)
    if validate:
        # Sifts seq has more res than the other, so we compare the common res
        seq_idx = ((~sifts_cif_dssp.cif_aa.isnull()) &
                   (sifts_cif_dssp.cif_aa != '-'))
        val_aa = sifts_cif_dssp.cif_aa[seq_idx]
        sifts_cif_dssp['sifts_aa'] = sifts_cif_dssp['REF_dbResName'].fillna('-')
        sifts_cif_dssp['sifts_aa'] = sifts_cif_dssp['sifts_aa'].apply(
            three_to_single_aa.get, args='X')
        val_aa2 = sifts_cif_dssp.sifts_aa[seq_idx]
        # Check if the sequences are the same
        if not (val_aa == val_aa2).all():
            raise ValueError('{pdb_id}|{chain} Cif and Sifts files have '
                             'different sequences '.format(pdb_id=pdb_id,
                                                           chain=chain))
    return sifts_cif_dssp


if __name__ == '__main__':
    from config import Defaults

    defaults = Defaults('config.txt')
    defaults.db_mmcif = 'tests/CIF'
    defaults.db_dssp = 'tests/DSSP'
    defaults.db_sifts = 'tests/SIFTS'

    X = merge_tables(pdb_id='3mn5', chain='A', default=defaults)
    pass
