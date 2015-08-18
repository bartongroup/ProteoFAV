#!/usr/bin/env python
# -*- coding: utf-8
import logging

import numpy as np

from to_table import select_cif, select_dssp, select_sifts
from utils import _fetch_sifts_best
from library import three_to_single_aa

log = logging.getLogger(__name__)


def merge_tables(uniprot_id=None, pdb_id=None, chain=None, model=1, validate=True):
    """
    Merges the output from multiple to table method.
    if no pdb_id uses sifts_best_structure
    if no chain uses first
    or if use_all_chains use all chains of the pdb_id
    :param uniprot_id: (opt) given a finds the best structure for a uniprot
    identifier.
    :param pdb_id: (opt) Alternatively parser
    :param chain: (opt) protein chain to parse. If None uses Sifts
    best_structure api to select best structure.
    :param model: number of the model to retrive from mmcif files.
    :param validate: Validate protein sequence between differnet data files.
    """
    if not any((uniprot_id, pdb_id)):
        raise TypeError("One of the following arguments is expected:"
                        "uniprot_id or pdb_id")

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

    cif_table = select_cif(pdb_id, chains=chain, models=model)
    dssp_table = select_dssp(pdb_id, chains=chain)
    cif_dssp = cif_table.join(dssp_table)

    if validate:# and cif_dssp['aa', 'label_comp_id'].any():
        # Fill all missing values with gaps
        cif_dssp['dssp_aa'] = cif_dssp.aa.fillna('-')
        # DSSP treat disulfite bonding cyteines as lower-cased pairs
        lower_cased_aa = cif_dssp.dssp_aa.str.islower()
        if lower_cased_aa.any():
            cif_dssp.loc[lower_cased_aa, 'dssp_aa'] = "C"
        cif_dssp['cif_aa'] = cif_dssp.label_comp_id.fillna('-')
        # From three letter to sigle letters or X if not a standard aa
        cif_dssp['cif_aa'] = cif_dssp['cif_aa'].apply(three_to_single_aa.get,
                                                      args='X')
        # Check if the sequences are the same
        if not (cif_dssp.dssp_aa == cif_dssp.cif_aa).all():
            raise ValueError('{pdb_id}|{chain} Cif and DSSP files have diffent'
                             ' sequences.'.format(pdb_id=pdb_id, chain=chain))

    sifts_table = select_sifts(pdb_id, chains=chain)
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
    X = merge_tables(pdb_id='3mg7', chain='A')

    pass
