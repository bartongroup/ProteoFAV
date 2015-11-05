#!/usr/bin/env python
# -*- coding: utf-8
import logging

import numpy as np
import pandas as pd

from to_table import (select_cif, select_dssp, select_sifts, select_validation,
                      _fetch_sifts_best, _uniprot_variants_to_table)
from library import to_single_aa

log = logging.getLogger(__name__)
logging.captureWarnings(True)
logging.basicConfig(level=9,
                    format='%(asctime)s - %(levelname)s - %(message)s ')


def merge_tables(uniprot_id=None, pdb_id=None, chain=None, model='first',
                 validate=True, add_validation=False, add_variants=False):
    """
    Merges the output from multiple to table method.
    if no pdb_id uses sifts_best_structure
    if no chain uses first
    or if use_all_chains use all chains of the pdb_id
    :param add_validation:
    :param add_validation:
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
    cif_table['auth_seq_id'] = cif_table.index
    try:
        dssp_table = select_dssp(pdb_id, chains=chain)
    except ValueError:
        # This indicates we could find the corret PDB chain.
        # What happens in case in structures with dozen of chains, ie: 4v9d
        # So we parse all the PDB file
        dssp_table = select_dssp(pdb_id)
        cif_seq = cif_table.auth_comp_id.apply(to_single_aa.get)
        dssp_table.reset_index(inplace=True)
        dssp_seq = "".join(dssp_table.aa)
        i = dssp_seq.find("".join(cif_seq))
        dssp_table = dssp_table.iloc[i: i + len(cif_seq), :]
        dssp_table.set_index(['icode'], inplace=True)
    # Correction for some dssp index parsed as object instead int
    if dssp_table.index.dtype == 'O' and cif_table.index.dtype != 'O':
        dssp_table.index = dssp_table.index.astype(np.int)

    table = cif_table.join(dssp_table)

    if validate:
        if not table['aa'].any():
            raise ValueError('Empty DSSP sequence cannot be validated')
        if not table['label_comp_id'].any():
            raise ValueError('Empty Cif sequence cannot be validated')
        table['dssp_aa'] = table.aa
        # DSSP disulfite bonding cyteines as lower-cased pairs
        # mask nans since they are not comparable
        mask = table['dssp_aa'].isnull()
        lower_cased_aa = table[~mask].dssp_aa.str.islower()
        if lower_cased_aa.any():
            table.loc[lower_cased_aa, 'dssp_aa'] = "C"
        table['cif_aa'] = table['label_comp_id']
        mask = mask | table['cif_aa'].isnull()
        # mask the X in DSSP since you can't compare those
        mask = mask | (table.dssp_aa == 'X')
        # From three letter to sigle letters or X if not a standard aa
        table['cif_aa'] = table['cif_aa'].apply(to_single_aa.get,
                                                args='X')
        # Check if the sequences are the same
        if not (table['dssp_aa'][~mask] == table['cif_aa'][~mask]).all():
            raise ValueError('{pdb_id}|{chain} Cif and DSSP files have diffent'
                             ' sequences.'.format(pdb_id=pdb_id, chain=chain))

    sifts_table = select_sifts(pdb_id, chains=chain)
    try:
        sifts_table.loc[:, 'PDB_dbResNum'] = sifts_table.loc[
                                             :, 'PDB_dbResNum'].astype(np.int)
        sifts_table.set_index(['PDB_dbResNum'], inplace=True)
    except ValueError:
        # Means it has alpha numeric insertion code, use something else as index
        sifts_table.loc[:, 'REF_dbResNum'] = sifts_table.loc[
                                             :, 'REF_dbResNum'].astype(np.int)
        sifts_table.set_index(['REF_dbResNum'], inplace=True)
        table.set_index(['pdbe_label_seq_id'], inplace=True)

    # Sift data in used as base to keep not observed residues info.
    table = sifts_table.join(table)
    if validate:
        if not table['REF_dbResName'].any():
            raise ValueError('Empty Sifts sequence cannot be validated')
        # Mask here because Sifts conseve missing residues and other data don't
        mask = table.cif_aa.isnull()
        table['sifts_aa'] = table['REF_dbResName']
        table['sifts_aa'] = table['sifts_aa'].apply(
            to_single_aa.get, args='X')
        mask = mask | table['sifts_aa'].isnull()

        # Check if the sequences are the same
        if not (table['cif_aa'][~mask] == table['sifts_aa'][~mask]).all():
            raise ValueError('{pdb_id}|{chain} Cif and Sifts files have '
                             'different sequences '.format(pdb_id=pdb_id,
                                                           chain=chain))
    if add_validation:
        validation_table = select_validation(pdb_id, chains=chain)
        validation_table.loc[:, 'val_resnum'] = validation_table.loc[
                                                :, 'val_resnum'].astype(np.int)
        validation_table.set_index(['val_resnum'], inplace=True)
        table = table.join(validation_table)
        if not table['val_resname'].any():
            raise ValueError('The validation file has empty sequence and cannot'
                             ' be validated')
        mask = table['cif_aa'].isnull()
        table['val_aa'] = table['val_resname']
        table['val_aa'] = table['val_aa'].apply(to_single_aa.get, args='X')
        mask = mask | table['val_aa'].isnull()
        if not (table['cif_aa'][~mask] == table['val_aa'][~mask]).all():
            raise ValueError('{pdb_id}|{chain} Cif and validation files have '
                             'different sequences '.format(pdb_id=pdb_id,
                                                           chain=chain))

    if add_variants:
        structure_uniprots = table.UniProt_dbAccessionId
        structure_uniprots = structure_uniprots[structure_uniprots.notnull()]
        structure_uniprot = structure_uniprots.unique()[0]
        variants_table = _uniprot_variants_to_table(structure_uniprot)
        variants_table[["start"]] = variants_table[["start"]].astype(float)
        table[["UniProt_dbResNum"]] = table[["UniProt_dbResNum"]].astype(float)

        table = table.reset_index()
        table = pd.merge(table, variants_table, left_on = "UniProt_dbResNum", right_on = "start", how = "left")

    return table


if __name__ == '__main__':
    X = merge_tables(pdb_id='2pm7', chain='D')
    X.head()
