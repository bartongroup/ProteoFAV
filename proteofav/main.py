#!/usr/bin/env python
# -*- coding: utf-8
import argparse
import click
import logging

import sys

from .library import to_single_aa
from .structures import (select_cif,
                         select_dssp,
                         select_sifts,
                         select_validation,
                         sifts_best)

from .variants import (map_gff_features_to_sequence,
                       select_variants)

__all__ = ['merge_tables',
           'parse_args',
           'main']
log = logging.getLogger('proteofav.config')



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

    :param pdb_id: Entry to be loaded
    :type pdb_id: str or None
    :param uniprot_id: Select PDB entry with highest sequence
        coverage for a UniProt protein sequence
    :type uniprot_id: str or None
    :param chain: Protein chain to loaded
    :type chain: str or None
    :param atoms: Atom to be selected in the
    :type atoms: str or None
    :param model: Select the PDB entity, like in structures determined by NMR
    :type model: str or None
    :param bool: add_validation: Attach the PDB validation table
    :param bool: sequence_check: Whether to compare sequence from different
        sources. Choose from raise, warn or ignore
    :param bool: drop_empty_cols: Whether to drop columns without positional
        information
    :param bool: add_ensembl_variants: Whether to add variant table from Ensembl
    :param bool: add_annotation: Whether to add variant table from Ensembl
    :param bool: add_uniprot_variants: Whether to add  variant table from
        UniProt
    :raises: ValueError
    :rtype: pd.DataFrame

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
        validation_table.loc[:, 'validation_resnum'] = validation_table.loc[:,
                                                       'validation_resnum'].astype(int)
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
            table = table.merge(variants_table,
                                how='left',
                                on=['UniProt_dbResNum', 'UniProt_dbAccessionId'])

    if add_annotation:
        for identifier in table['UniProt_dbAccessionId'].dropna().unique():
            uniprot_annotation = map_gff_features_to_sequence(identifier)
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


@click.command()
@click.option('--pdb', type=str, default=None, help='Protein data bank identifier.')
@click.option('--chain', type=str, default=None, help='Protein structure chain identifier.')
@click.option('--uniprot', type=str, default=None, help='UniProt knowledgebase accession.')
@click.option('--output_type', default='csv',
              type=click.Choice(['csv', 'json', 'tab']),  # TODO add JALVIEW and chimera
              help='File format for the output.')
@click.option('-v', '--verbose', is_flag=True, help="Show more verbose logging.")
@click.option('-l', '--log', default=sys.stderr, help="Path to the logfile.")
@click.option('--add_annotation', is_flag=True,
              help="Whether to merge annotation information to the output.")
@click.option('--add_validation', is_flag=True,
              help="Whether to merge protein structure validation information to the output.")
@click.option('--add_variants', is_flag=True,
              help="Whether to merge genetic variant information to the output.")
@click.option('--remove_redundants', is_flag=True,
              help="Whether to remove columns with redundant information from the output.")
@click.argument ('output', type=click.File('wb'))
def main(pdb, chain, uniprot, output_type, verbose, log, add_annotation, add_validation,
         add_variants, remove_redundants, output):
    """
    ProteFAV: a Python framework to process and integrate protein structure and
    features to genetic variants. 
    OUTPUT: Path to the output file. Use `-` for stdout
    """

    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(stream=log, level=level,
                        format='%(asctime)s - %(levelname)s - %(message)s ')

    table = merge_tables(pdb_id=pdb,
                         chain=chain,
                         uniprot_id=uniprot,
                         add_annotation=add_annotation,
                         add_validation=add_validation,
                         add_ensembl_variants=add_variants,
                         drop_empty_cols=remove_redundants)

    if output_type == 'csv':
        table.to_csv(output)
    elif output_type in {'jalview', 'chimera'}:
        raise NotImplementedError
    else:
        if hasattr(table, "to_" + output_type):
            fun = getattr(table, "to_" + output_type)
            fun(output)
    return 0


if __name__ == '__main__':
    sys.exit(main())
