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
                 validate=True, drop_empty_cols=False, add_validation=False, add_annotation=False,
                 add_ensembl_variants=None, add_uniprot_variants=False):
    """
    Join multiple resource tables. If no pdb_id uses sifts_best_structure
    If no chain uses the first on.

    :param str or None uniprot_id: gives sifts best representative for this UniProt entry
    :param str or None pdb_id: Entry to be loaded
    :param str or None chain: Protein chain to loaded
    :param str or None atoms: Atom to be selected in the
    :param str or None model: Select the PDB entity, like in structures determined by NMR
    :param bool add_validation: Attach the PDB validation table
    :param bool validate: Whether to validate the protein sequence in different tables
    :param bool drop_empty_cols: Whether to drop columns without positional information
    :param bool add_ensembl_variants: Whether to add variant table from Ensembl
    :param bool add_annotation: Whether to add variant table from Ensembl
    :param bool add_uniprot_variants: Whether to add  variant table from UniProt
    :rtype: pandas.DataFrame

    """
    if not any((uniprot_id, pdb_id)):
        raise TypeError("One of the following arguments is expected:"
                        "uniprot_id or pdb_id")

    if not pdb_id:
        best_pdb = sifts_best(uniprot_id, first=True)
        if best_pdb is None:
            logging.error('Could not process {}'.format(uniprot_id))
            return None
        pdb_id = best_pdb['pdb_id']
        chain = best_pdb['chain_id']
        log.info("Best structure, chain: {}|{} for {} ".format(pdb_id, chain, uniprot_id))

    cif_table = select_cif(pdb_id, chains=chain, models=model, atoms=atoms)
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

    # pragma
    table.loc[:, 'UniProt_dbResNum'] = table.loc[:, 'UniProt_dbResNum'].fillna(-9999)
    table.loc[:, 'UniProt_dbResNum'] = table.loc[:, 'UniProt_dbResNum'].astype(int)

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
            variants_table.rename(columns={'index': 'UniProt_dbResNum'}, inplace=True)
            table = table.merge(uniprot_annotation,
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
    X2 = merge_tables(pdb_id='2pah', atoms=None)
    print(X.head())
