#!/usr/bin/env python
# -*- coding: utf-8

from __future__ import print_function

import logging

from proteofav.library import to_single_aa
from proteofav.structures import (select_cif, select_dssp, select_sifts,
                                  select_validation, sifts_best)
from proteofav.variants import (map_gff_features_to_sequence, select_variants)

log = logging.getLogger(__name__)
logging.captureWarnings(True)
logging.basicConfig(level=9,
                    format='%(asctime)s - %(levelname)s - %(message)s ')


def merge_tables(uniprot_id=None, pdb_id=None, chain=None, atoms='CA', model='first',
                 sequence_check='raise', drop_empty_cols=False, add_validation=False,
                 add_annotation=False,
                 add_ensembl_variants=None, add_uniprot_variants=False):
    """
    Join multiple resource tables. If no pdb_id uses sifts_best_structure
    If no chain uses the first on.

    :param str or None uniprot_id: Fetch Sifts best (higher coverage) protein structure for
    this UniProt entry
    :param str or None pdb_id: Entry to be loaded
    :param str or None chain: Protein chain to loaded
    :param str or None atoms: Atom to be selected in the
    :param str or None model: Select the PDB entity, like in structures determined by NMR
    :param bool add_validation: Attach the PDB validation table
    :param str sequence_check: Whether to compare sequence from different sources.
    Choose from raise, warn or ignore.
    :param bool drop_empty_cols: Whether to drop columns without positional information
    :param bool add_ensembl_variants: Whether to add variant table from Ensembl
    :param bool add_annotation: Whether to add variant table from Ensembl
    :param bool add_uniprot_variants: Whether to add  variant table from UniProt
    :rtype: pandas.DataFrame

    """
    if not any((uniprot_id, pdb_id)):
        raise TypeError('One of the following arguments is expected:'
                        'uniprot_id or pdb_id')

    if model != 'first':
        raise NotImplementedError('Proteofav current implementation ignore alternative models.')

    if sequence_check not in ['raise', 'warn ' or 'ignore']:
        sequence_check = 'raise'

    if not pdb_id:
        best_pdb = sifts_best(uniprot_id, first=True)
        if best_pdb is None:
            logging.error('Could not process {}'.format(uniprot_id))
            return None
        pdb_id = best_pdb['pdb_id']
        chain = best_pdb['chain_id']
        log.info('Best structure, chain: {}|{} for {} '.format(pdb_id, chain, uniprot_id))

    cif_table = select_cif(pdb_id, chains=chain, models=model, atoms=atoms)

    dssp_table = select_dssp(pdb_id, chains=chain)

    cif_table.loc[:, 'label_seq_id'] = cif_table.loc[:, 'label_seq_id'].astype(int)
    table = cif_table.merge(dssp_table, how='left',
                            left_on=['auth_seq_id', 'auth_asym_id'],
                            right_on=['icode', 'chain_id'])

    if sequence_check == 'ignore' or atoms is None:
        # sequence check not support for multiple atoms
        pass
    else:
        # exchange lower cased aa's for  for cysteines
        lower_cased_aa = table['aa'].str.islower()
        if lower_cased_aa.any():
            table.loc[lower_cased_aa, 'aa'] = 'C'

        mask = table['label_comp_id'].isnull() | table['aa'].isnull()
        mask |= table['aa'] == 'X'
        cif_seq = table['label_comp_id'].apply(to_single_aa.get, args='X')
        dssp_seq = table['aa']

        # Check if the sequences are the same
        if not (dssp_seq[~mask] == cif_seq[~mask]).all():
            log_msg = '{}|{} Cif and DSSP files have different sequences.'.format(pdb_id, chain)
            if sequence_check == 'raise':
                raise ValueError(log_msg)
            else:
                log.warn(log_msg)

    sifts_table = select_sifts(pdb_id, chains=chain)
    try:
        sifts_table.loc[:, 'PDB_dbResNum'] = sifts_table.loc[:, 'PDB_dbResNum'].astype(int)
        table = sifts_table.merge(table, how='left',
                                  left_on=['PDB_dbResNum', 'PDB_dbChainId'],
                                  right_on=['auth_seq_id', 'auth_asym_id'])

    except ValueError:
        # PDB resnumber has a insertion code
        table.pdbx_PDB_ins_code = table.pdbx_PDB_ins_code.replace('?', '')
        table['index'] = table.auth_seq_id.astype(str) + table.pdbx_PDB_ins_code
        table = table.merge(sifts_table, how='left',
                            left_on=['index', 'auth_asym_id'],
                            right_on=['PDB_dbResNum', 'PDB_dbChainId'])

    if sequence_check == 'ignore' or atoms is None:
        pass
    else:
        # Update reference sequence with the new table
        mask = table['auth_comp_id'].isnull() | table['PDB_dbResName'].isnull()
        cif_seq = table['auth_comp_id'].apply(to_single_aa.get, args='X')
        sifts_seq = table['PDB_dbResName'].apply(to_single_aa.get, args='X')

        # Check if the sequences are the same
        if not (sifts_seq[~mask] == cif_seq[~mask]).all():
            log_msg = '{}|{} Cif and Sifts files have different sequences.'.format(pdb_id, chain)
            if sequence_check == 'raise':
                raise ValueError(log_msg)
            else:
                log.warn(log_msg)

    if add_validation:
        validation_table = select_validation(pdb_id, chains=chain)
        validation_table.loc[:, 'validation_resnum'] = validation_table.loc[
                                                       :, 'validation_resnum'].astype(int)
        table = table.merge(table, how='left',
                            left_on=['PDB_dbResNum', 'PDB_dbChainId'],
                            right_on=['validation_resnum', 'validation_chain'])

    variant_features = []
    if add_ensembl_variants:
        variant_features.extend(['ensembl_somatic', 'ensembl_germline'])
    if add_uniprot_variants:
        variant_features.append(['uniprot'])
    if variant_features:
        for identifier in table['UniProt_dbAccessionId'].dropna().unique():
            variants_table = select_variants(identifier, features=variant_features)
            variants_table.reset_index(inplace=True)
            variants_table['UniProt_dbAccessionId'] = identifier
            variants_table.rename(columns={'start': 'UniProt_dbResNum'}, inplace=True)
            table = table.merge(variants_table, how='left',
                                on=['UniProt_dbResNum', 'UniProt_dbAccessionId'])

    if add_annotation:
        for identifier in table['UniProt_dbAccessionId'].dropna().unique():
            uniprot_annotation = map_gff_features_to_sequence(identifier)
            uniprot_annotation.reset_index(inplace=True)
            uniprot_annotation['UniProt_dbAccessionId'] = identifier
            uniprot_annotation.rename(columns={'index': 'UniProt_dbResNum'}, inplace=True)
            table = table.merge(uniprot_annotation, how='left',
                                on=['UniProt_dbAccessionId', 'UniProt_dbResNum'])

    # non-positional information goes to an attribute
    if drop_empty_cols:
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
    pass
