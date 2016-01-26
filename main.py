#!/usr/bin/env python
# -*- coding: utf-8

from __future__ import print_function

import logging

from library import to_single_aa
from pandas import DataFrame
from structures.to_table import (select_cif, select_dssp, select_sifts,
                                 select_validation, sifts_best, _rcsb_description)

from variants.to_table import (select_uniprot_gff, select_uniprot_variants,
                               _fetch_uniprot_variants)

log = logging.getLogger(__name__)
logging.captureWarnings(True)
logging.basicConfig(level=9,
                    format='%(asctime)s - %(levelname)s - %(message)s ')

__all__ = ["merge_tables"]


def merge_tables(uniprot_id=None, pdb_id=None, chain=None, model='first',
                 validate=True, add_validation=False, add_variants=None,
                 add_annotation=False, remove_redundant=False,
                 uniprot_variants=False):
    """
    Join multiple resource tables. If no pdb_id uses sifts_best_structure
    If no chain uses the first on.

    :type remove_redundant: bool
    :type add_annotation: bool
    :type add_variants: bool
    :type add_validation: bool
    :type validate: bool
    :type model: str or None
    :type chain: str or None
    :type pdb_id: str or None
    :type uniprot_id: str or None
    :rtype: pandas.DataFrame
    :param add_validation: join the PDB validation table?
    :param uniprot_id: gives sifts best representative for this UniProt entry
    :param pdb_id: Entry to be parsed
    :param chain: Protein chain to be parsed
    :param model: Which entity to use? Useful for multiple entities protein
    structures determined by NMR
    :param validate: Checks whether sequence is the same in all tables
    :type uniprot_variants: bool
    """
    if not any((uniprot_id, pdb_id)):
        raise TypeError("One of the following arguments is expected:"
                        "uniprot_id or pdb_id")

    if chain == 'all':
        # If we want to fetch all chains we can just find out what chains are
        # available in the specified PDB and then recursively call `merge_tables`
        # until we got them all. This should only work for a PDB based query and we
        # need to ensure that fields like 'chain_id' aren't dropped by 'remove_redundant'
        if not pdb_id:
            best_pdb = sifts_best(uniprot_id, first=True)
            pdb_id = best_pdb['pdb_id']
            if best_pdb is None:
                logging.error('Could not process {}'.format(uniprot_id))
                return None
            log.info("Best structure: {} for {} ".format(pdb_id, uniprot_id))
        if remove_redundant:
            remove_redundant = False
            log.warning("remove_redundant is ignored when chain='all' and is set to False")

        chain_ids = _rcsb_description(pdb_id, tag='chain', key='id')
        # TODO transition to chain = None strategy
        table = DataFrame()
        for current_chain in chain_ids:
            table = table.append(merge_tables(uniprot_id=uniprot_id, pdb_id=pdb_id,
                                              chain=current_chain, model=model,
                                              validate=validate,
                                              add_validation=add_validation,
                                              add_variants=add_variants,
                                              add_annotation=add_annotation,
                                              remove_redundant=remove_redundant,
                                              uniprot_variants=uniprot_variants))

        return table

    if not pdb_id:
        best_pdb = sifts_best(uniprot_id, first=True)
        if best_pdb is None:
            logging.error('Could not process {}'.format(uniprot_id))
            return None
        pdb_id = best_pdb['pdb_id']
        chain = best_pdb['chain_id']
        log.info("Best structure, chain: {}|{} for {} ".format(pdb_id, chain, uniprot_id))

    cif_table = select_cif(pdb_id, chains=chain, models=model)
    cif_table['auth_seq_id'] = cif_table.index

    dssp_table = select_dssp(pdb_id, chains=chain)
    if dssp_table.index.dtype == 'O' and cif_table.index.dtype != 'O':
        dssp_table.index = dssp_table.index.astype(int)

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
        dssp_aa = table.dssp_aa.fillna('X')
        lower_cased_aa = dssp_aa.str.islower()
        if lower_cased_aa.any():
            table.loc[lower_cased_aa, 'dssp_aa'] = "C"
        table['cif_aa'] = table['label_comp_id']
        mask |= table['cif_aa'].isnull()
        # mask the X in DSSP since you can't compare those
        mask |= table.dssp_aa == 'X'
        # From three letter to single letters or X if not a standard aa
        table['cif_aa'] = table['cif_aa'].apply(to_single_aa.get, args='X')
        # Check if the sequences are the same
        if not (table['dssp_aa'][~mask] == table['cif_aa'][~mask]).all():
            err = '{}|{} Cif and DSSP files have different sequences.'
            raise ValueError(err.format(pdb_id, chain))

    sifts_table = select_sifts(pdb_id, chains=chain)
    try:
        sifts_table.loc[:, 'PDB_dbResNum'] = sifts_table.loc[:, 'PDB_dbResNum'].astype(int)
        # TODO: Shouldn't the index be strictly defined no matter what?
        sifts_table.set_index(['PDB_dbResNum'], inplace=True)
    except ValueError:
        # Means it has alpha numeric insertion code, use something else as index
        sifts_table.set_index(['PDB_dbResNum'], inplace=True)
        table.pdbx_PDB_ins_code = table.pdbx_PDB_ins_code.replace('?', '')

        table['index'] = table.auth_seq_id.astype(str) + table.pdbx_PDB_ins_code
        # See Above
        table.set_index(['index'], inplace=True)

    # Sift data in used as base to keep not observed residues info.
    # TODO: this is inconsistent with other joins; should table be left?
    table = sifts_table.join(table)
    if validate:
        if not table['REF_dbResName'].any():
            raise ValueError('Empty Sifts sequence cannot be validated')
        # Mask here because Sifts conserve missing residues and other data don't
        mask = table.cif_aa.isnull()
        table['sifts_aa'] = table['REF_dbResName']
        table['sifts_aa'] = table['sifts_aa'].apply(to_single_aa.get, args='X')
        mask = mask | table['sifts_aa'].isnull()

        # Check if the sequences are the same
        if not (table['cif_aa'][~mask] == table['sifts_aa'][~mask]).all():
            err = '{}|{} Cif and Sifts files have different sequences '
            raise ValueError(err.format(pdb_id, chain))
    if add_validation:
        validation_table = select_validation(pdb_id, chains=chain)
        validation_table.loc[:, 'val_resnum'] = validation_table.loc[:, 'val_resnum'].astype(int)
        validation_table.set_index(['val_resnum'], inplace=True)
        table = table.join(validation_table)
        if not table['val_resname'].any():
            raise ValueError('The validation file has empty sequence and cannot  be validated')
        mask = table['cif_aa'].isnull()
        table['val_aa'] = table['val_resname']
        table['val_aa'] = table['val_aa'].apply(to_single_aa.get, args='X')
        mask = mask | table['val_aa'].isnull()
        if not (table['cif_aa'][~mask] == table['val_aa'][~mask]).all():
            err = '{}|{} Cif and validation files have different sequences '
            raise ValueError(err.format(pdb_id, chain))

    # Table preparations before adding variants
    table[["UniProt_dbResNum"]] = table[["UniProt_dbResNum"]].astype(float)
    no_mapped_uniprot = table[table.UniProt_dbAccessionId.isnull()]

    if add_variants:
        grouped_table = table.groupby('UniProt_dbAccessionId')
        new_table = DataFrame()
        for structure_uniprot, part_table in grouped_table:
            variants_table = select_uniprot_variants(structure_uniprot)
            variants_table[["start"]] = variants_table[["start"]].astype(float)
            part_table = part_table.reset_index()  # Gives access to UniProt_dbResNum
            merged_table = part_table.merge(variants_table, left_on="UniProt_dbResNum", right_on="start", how="left")
            if new_table.empty:
                col_order = merged_table.columns.tolist()  # Get the column order once
            new_table = new_table.append(merged_table)
        table = new_table.append(no_mapped_uniprot)[col_order]
        table.set_index(['PDB_dbResNum'], inplace=True)

    if uniprot_variants:
        grouped_table = table.groupby('UniProt_dbAccessionId')
        new_table = DataFrame()
        for structure_uniprot, part_table in grouped_table:
            variants = _fetch_uniprot_variants(structure_uniprot)
            part_table = part_table.reset_index()  # Gives access to UniProt_dbResNum
            merged_table = part_table.merge(variants, on='UniProt_dbResNum', how='left')
            if new_table.empty:
                col_order = merged_table.columns.tolist()  # Get the column order once
            new_table = new_table.append(merged_table)
        table = new_table.append(no_mapped_uniprot)[col_order]
        table.set_index(['PDB_dbResNum'], inplace=True)

    if add_annotation:
        for identifier in table['UniProt_dbAccessionId'].dropna().unique():
            uniprot_annotation = select_uniprot_gff(identifier)
            uniprot_annotation.index = uniprot_annotation.index.astype(float)
            table = table.reset_index().merge(uniprot_annotation, how='left',
                                              on="UniProt_dbResNum")
            table.set_index(['PDB_dbResNum'], inplace=True)

    # non-positional information goes to an attribute
    if remove_redundant:
        for col in table:
            try:
                value = table[col].dropna().unique()
            except TypeError:
                # break for list-like columns
                continue

            if value.shape[0] == 1:
                if value[0] == '?':
                    # if the only value is a `?` we don't need to keep
                    continue
                del table[col]
                setattr(table, col, value)
                log.info('Column {} is now an attribute.'.format(col))
    return table


if __name__ == '__main__':
    X = merge_tables(pdb_id='4v9d', chain='BD')
    print(X.head())
