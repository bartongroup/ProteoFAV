# -*- coding: utf-8

import logging
import numpy as np
import pandas as pd

from proteofav.structures import select_structures
from proteofav.validation import select_validation
from proteofav.sifts import select_sifts
from proteofav.sifts import sifts_best
from proteofav.dssp import select_dssp
from proteofav.variants import select_variants
from proteofav.annotation import select_annotation
from proteofav.utils import merging_down_by_key
from proteofav.library import to_single_aa

__all__ = ['uniprot_vars_ensembl_vars_merger',
           'merge_tables']

log = logging.getLogger('proteofav.config')


class TableMergerError(Exception):
    pass


def uniprot_vars_ensembl_vars_merger(uniprot_vars_table, ensembl_vars_table):
    """
    Merges the tables provided using appropriate columns.

    :param uniprot_vars_table: UniProt Variants pandas DataFrame
    :param ensembl_vars_table: Ensembl Variants pandas DataFrame
    :return: merged pandas DataFrame
    """

    # bare minimal columns needed
    merge_on = ['begin', 'end', 'xrefs_id', 'frequency',
                'consequenceType', 'siftScore', 'polyphenScore']

    if (set(merge_on).issubset(uniprot_vars_table.columns) and
            set(merge_on).issubset(ensembl_vars_table.columns)):

        table = uniprot_vars_table.merge(ensembl_vars_table, how='outer',
                                         on=merge_on).reset_index(drop=True)

        table = merging_down_by_key(table, key='xrefs_id')
        table.fillna(np.nan, inplace=True)
    else:
        raise TableMergerError('Not possible to merge UniProt and Ensembl Vars table! '
                               'Some of the necessary columns are missing...')

    log.info("Merged UniProt and Ensembl Vars tables...")
    return table


def merge_tables(uniprot_id=None,
                 pdb_id=None,
                 chain=None,
                 atoms='CA',
                 model='first',
                 sequence_check='raise',
                 drop_empty_cols=False,
                 add_validation=False,
                 add_annotation=False,
                 add_ensembl_variants=False,
                 add_uniprot_variants=False):
    """
    Automatically merge data tables from various resources. If no pdb_id set
        sifts_best_structure, which is sorted by sequence coverage.
        If no chain is set, uses the first one.

    :param bool add_ensembl_variants: Whether to add variant table from Ensembl
    :param bool add_validation: whether to merge PDB validation information
    :param bool add_annotation: whether to merge UniProt GFF information
    :param bool drop_empty_cols: whether to drop columns without useful information
    :param str sequence_check: how to handle sequence inconsistencies. Choose from raise,
    warn or ignore
    :param pdb_id: Entry to be loaded
    :type pdb_id: str or None
    :param uniprot_id: Select PDB entry with highest sequence
        coverage for the selected UniProt protein sequence
    :type uniprot_id: str or None
    :param chain: Protein chain to loaded
    :type chain: str or None
    :param atoms: Atom to be selected in the
    :type atoms: str or None
    :param model: Select the PDB entity, like in structures determined by NMR
    :type model: str or None
    """
    if not any((uniprot_id, pdb_id)):
        raise TypeError('One of the following arguments is expected:'
                        'uniprot_id or pdb_id')

    if model != 'first':
        raise NotImplementedError('Proteofav current implementation ignore alternative models.')

    if sequence_check not in ['raise', 'warn', 'ignore']:
        sequence_check = 'raise'

    if not pdb_id:
        best_pdb = sifts_best(uniprot_id, first=True)
        if best_pdb is None:
            logging.error('Could not process {}'.format(uniprot_id))
            return None
        pdb_id = best_pdb['pdb_id']
        chain = best_pdb['chain_id']
        log.info('Best structure, chain: {}|{} for {} '.format(pdb_id, chain, uniprot_id))

    cif_table = select_structures(pdb_id, chains=chain, models=model, atoms=atoms)

    dssp_table = select_dssp(pdb_id, chains=chain)

    cif_table.loc[:, 'auth_seq_id'] = cif_table.loc[:, 'auth_seq_id'].astype(str)
    dssp_table.loc[:, 'RES_FULL'] = dssp_table.loc[:, 'RES_FULL'].astype(str)
    table = cif_table.merge(dssp_table, how='left',
                            left_on=['auth_seq_id', 'auth_asym_id'],
                            right_on=['RES_FULL', 'CHAIN'])

    if sequence_check == 'ignore' or atoms is None:
        # sequence check not support for multiple atoms
        pass
    else:
        # exchange lower cased aa's for  for cysteines
        lower_cased_aa = table['AA'].str.islower()
        if lower_cased_aa.any():
            table.loc[lower_cased_aa, 'AA'] = 'C'

        mask = table['label_comp_id'].isnull() | table['AA'].isnull()
        mask |= table['AA'] == 'X'
        cif_seq = table['label_comp_id'].apply(to_single_aa.get, args='X')
        dssp_seq = table['AA']

        # Check if the sequences are the same
        if not (dssp_seq[~mask] == cif_seq[~mask]).all():
            log_msg = '{}|{} Cif and DSSP files have different sequences.'.format(pdb_id, chain)
            if sequence_check == 'raise':
                raise ValueError(log_msg)
            else:
                log.warning(log_msg)

    sifts_table = select_sifts(pdb_id, chains=chain)
    try:
        sifts_table.loc[:, 'PDB_dbResNum'] = sifts_table.loc[:, 'PDB_dbResNum'].astype(str)
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
                log.warning(log_msg)

    if add_validation:
        validation_table = select_validation(pdb_id, chains=chain)
        validation_table.loc[:, 'validation_resnum'] = validation_table.loc[:,
                                                       'validation_resnum'].astype(int)
        table = table.merge(table, how='left',
                            left_on=['PDB_dbResNum', 'PDB_dbChainId'],
                            right_on=['validation_resnum', 'validation_chain'])

    if add_uniprot_variants or add_ensembl_variants:
        for identifier in table['UniProt_dbAccessionId'].dropna().unique():
            uniprot_vars, ensembl_vars = select_variants(identifier,
                                                         uniprot_vars=add_uniprot_variants,
                                                         ensembl_germline_vars=add_ensembl_variants,
                                                         ensembl_somatic_vars=add_ensembl_variants)

            if isinstance(uniprot_vars, pd.DataFrame) and isinstance(ensembl_vars, pd.DataFrame):
                variants_table = uniprot_vars_ensembl_vars_merger(uniprot_vars, ensembl_vars)
            elif isinstance(uniprot_vars, pd.DataFrame):
                variants_table = uniprot_vars
            elif isinstance(ensembl_vars, pd.DataFrame):
                variants_table = ensembl_vars
            else:
                log.info('No variants found...')
                variants_table = pd.DataFrame()

            if not variants_table.empty:
                variants_table.reset_index(inplace=True)
                variants_table['UniProt_dbAccessionId'] = identifier
                variants_table.rename(columns={'begin': 'UniProt_dbResNum'}, inplace=True)
                table = table.merge(variants_table,
                                    how='left',
                                    on=['UniProt_dbResNum', 'UniProt_dbAccessionId'])

    if add_annotation:
        for identifier in table['UniProt_dbAccessionId'].dropna().unique():
            uniprot_annotation = select_annotation(identifier, annotation_agg=True)
            uniprot_annotation.reset_index(inplace=True)
            uniprot_annotation['UniProt_dbAccessionId'] = identifier
            uniprot_annotation.rename(columns={'index': 'UniProt_dbResNum'}, inplace=True)
            table = table.merge(uniprot_annotation,
                                how='left',
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
