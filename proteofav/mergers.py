# -*- coding: utf-8

import logging
import numpy as np
import pandas as pd

from proteofav.structures import select_structures, mmCIF
from proteofav.sifts import select_sifts, sifts_best, SIFTS
from proteofav.dssp import select_dssp, DSSP
from proteofav.variants import select_variants, Variants
from proteofav.annotation import select_annotation, Annotation
from proteofav.validation import select_validation, Validation
from proteofav.utils import merging_down_by_key
from proteofav.library import to_single_aa

__all__ = ['uniprot_vars_ensembl_vars_merger',
           'merge_tables']

log = logging.getLogger('proteofav.config')


class TableMergerError(Exception):
    pass


def mmcif_sifts_table_merger(mmcif_table, sifts_table):
    """
    Merge the mmCIF and SIFTS Tables.

    :param mmcif_table: mmCIF pandas DataFrame
    :param sifts_table: SIFTS pandas DataFrame
    :return: merged pandas DataFrame
    """

    # bare minimal columns needed
    if ('auth_seq_id_full' in mmcif_table and 'label_asym_id' in mmcif_table and
                'PDB_dbResNum' in sifts_table and 'PDB_entityId' in sifts_table):

        # workaround for BioUnits
        if 'orig_label_asym_id' in mmcif_table:
            table = mmcif_table.merge(sifts_table, how='left',
                                      left_on=['auth_seq_id_full', 'orig_label_asym_id'],
                                      right_on=['PDB_dbResNum', 'PDB_entityId'])
        else:
            table = mmcif_table.merge(sifts_table, how='left',
                                      left_on=['auth_seq_id_full', 'label_asym_id'],
                                      right_on=['PDB_dbResNum', 'PDB_entityId'])

    else:
        raise TableMergerError('Not possible to merge mmCIF and SIFTS table! '
                               'Some of the necessary columns are missing...')

    log.info("Merged mmCIF and SIFTS Tables...")
    return table


def mmcif_dssp_table_merger(mmcif_table, dssp_table):
    """
    Merge the mmCIF and DSSP Tables.

    :param mmcif_table: mmCIF pandas DataFrame
    :param dssp_table: DSSP pandas DataFrame
    :return: merged pandas DataFrame
    """

    # bare minimal columns needed
    if ('auth_seq_id_full' in mmcif_table and 'label_asym_id' in mmcif_table and
                'RES' in dssp_table and 'CHAIN_FULL' in dssp_table):

        table = mmcif_table.merge(dssp_table, how='left',
                                  left_on=['auth_seq_id_full', 'label_asym_id'],
                                  right_on=['RES', 'CHAIN_FULL'])
    else:
        print()
        print(list(mmcif_table))
        print(list(dssp_table))
        print()
        raise TableMergerError('Not possible to merge mmCIF and DSSP table! '
                               'Some of the necessary columns are missing...')

    log.info("Merged mmCIF and DSSP Tables...")
    return table


def dssp_sifts_table_merger(dssp_table, sifts_table):
    """
    Merge the DSSP and SIFTS Tables.

    :param dssp_table: DSSP pandas DataFrame
    :param sifts_table: SIFTS pandas DataFrame
    :return: merged pandas DataFrame
    """

    # bare minimal columns needed
    if ('RES' in dssp_table and 'CHAIN' in dssp_table and
                'PDB_dbResNum' in sifts_table and 'PDB_entityId' in sifts_table):

        table = dssp_table.merge(sifts_table, how='left',
                                 left_on=['RES', 'CHAIN'],
                                 right_on=['PDB_dbResNum', 'PDB_entityId'])

    else:
        raise TableMergerError('Not possible to merge DSSP and SIFTS table! '
                               'Some of the necessary columns are missing...')

    log.info("Merged DSSP and SIFTS Tables...")
    return table


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
    Automatically merge data Tables from various resources. If no pdb_id set
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


# FIXME - add missing
def table_merger(mmcif_table=None, dssp_table=None, sifts_table=None,
                 validation_table=None, annotation_table=None,
                 variants_table=None):
    """
    Merges the Tables provided using appropriate columns.

    :param mmcif_table: (optional) mmCIF pandas DataFrame
    :param dssp_table: (optional) DSSP pandas DataFrame
    :param sifts_table: (optional) SIFTS pandas DataFrame
    :param validation_table: (optional) Validation pandas DataFrame
    :param annotation_table: (optional) Annotation pandas DataFrame
    :param variants_table: (optional) Variants pandas DataFrame
    :return: merged pandas DataFrame
    """

    available = [mmcif_table, dssp_table, sifts_table]
    available = [k for k in available if k is not None]
    if len(available) < 2:
        raise TableMergerError(
            "At least two Tables are needed in order to merge...")

    table = None
    if mmcif_table is not None:
        if dssp_table is not None:
            mmcif_table = mmcif_dssp_table_merger(mmcif_table, dssp_table)
        if sifts_table is not None:
            mmcif_table = mmcif_sifts_table_merger(mmcif_table, sifts_table)
        table = mmcif_table

    elif dssp_table is not None:
        if sifts_table is not None:
            dssp_table = dssp_sifts_table_merger(dssp_table, sifts_table)
        table = dssp_table

    return table


def table_generator(uniprot_id=None, pdb_id=None, bio_unit=False,
                    sifts=True, dssp=False, variants=False, annotations=False,
                    chains=None, res=None, sites=None, atoms=None, lines=None,
                    residue_agg=False, overwrite=False):
    """
    Simplifies the process of generating Tables and merging them.

    :param uniprot_id: (str) UniProt ID
    :param pdb_id: (str) PDB ID
    :param bio_unit: boolean for using AsymUnits or BioUnits
    :param sifts: boolean
    :param dssp: boolean
    :param variants: boolean
    :param annotations: boolean
    :param chains: (tuple) Chain ID ('*_asym_id')
    :param res: (tuple) PDB ResNum ('*_seq_id_full')
    :param sites: (tuple) UniProt ResNum (positional index)
    :param atoms: (tuple) atom IDs or None
    :param lines: 'ATOM' or 'HETATM' or None (both)
    :param residue_agg: boolean
    :param overwrite: boolean
    :returns: mmcif, sifts, dssp, variants and annotations Tables
    """

    # generates the Tables for the uniprot_id, pdb_id
    # chains, res, sites are only processed after
    if uniprot_id or pdb_id:
        # get pdb_id if only uniprot_id is provided
        if uniprot_id and not pdb_id:
            # fetch the best structures from the PDBe API
            data = sifts_best(uniprot_id)
            if data is not None:
                # uses the first structure
                pdb_id = data[0]['pdb_id']
                chains = (data[0]['chain_id'],)
            else:
                log.info("Best structures not available from the PDBe API for %s",
                         uniprot_id)
                raise TableMergerError('Nothing to merge...')

        # mmCIF table
        mmcif_table = mmCIF.select(identifier=pdb_id, bio_unit=bio_unit,
                                   bio_unit_preferred=True, overwrite=overwrite,
                                   add_atom_altloc=True, add_res_full=True,
                                   category='auth', residue_agg=residue_agg,
                                   chains=chains, res_full=res, atoms=atoms, lines=lines)

        # SIFTS table
        if sifts:
            sifts_table = SIFTS.select(identifier=pdb_id, add_regions=True, add_dbs=False,
                                       chains=chains, res=res, uniprot=uniprot_id, site=sites)
        else:
            sifts_table = None

        # DSSP table
        if dssp:
            dssp_table = DSSP.select(identifier=pdb_id,
                                     add_ss_reduced=True, add_rsa_class=True,
                                     chains_full=chains, res=res)
        else:
            dssp_table = None

        # Validation table
        if annotations:
            validation_table = Validation.select(identifier=uniprot_id)
        else:
            validation_table = None

        # Variants table
        if variants:
            uni_vars, ens_vars = Variants.select(uniprot_id, id_source='uniprot',
                                                 synonymous=True,
                                                 uniprot_vars=True,
                                                 ensembl_germline_vars=True,
                                                 ensembl_somatic_vars=True)
            variants_table = uniprot_vars_ensembl_vars_merger(uni_vars, ens_vars)
        else:
            variants_table = None

        # Annotations table
        if annotations:
            annotations_table = Annotation.select(identifier=uniprot_id)
        else:
            annotations_table = None

        return (mmcif_table, dssp_table, sifts_table,
                validation_table, annotations_table, variants_table)

    else:
        raise ValueError('No UniProt ID or PDB ID provided...')


class Tables(object):
    def __init__(self):
        self.mmcif = None
        self.sifts = None
        self.dssp = None
        self.validation = None
        self.annotation = None
        self.variants = None
        self.table = None

    def merge(self, mmcif=None, dssp=None, sifts=None,
              validation=None, annotation=None, variants=None):
        """Merges the provided Tables and stores"""

        if mmcif is not None:
            self.mmcif = mmcif
        if dssp is not None:
            self.dssp = dssp
        if sifts is not None:
            self.sifts = sifts
        if validation is not None:
            self.validation = validation
        if annotation is not None:
            self.annotation = annotation
        if variants is not None:
            self.variants = variants
        self.table = table_merger(self.mmcif, self.dssp, self.sifts,
                                  self.validation, self.annotation, self.variants)
        return self.table

    def generate(self, **kwargs):
        """Generates the Tables, merges and stores"""
        self.mmcif, self.dssp, self.sifts, self.validation, self.annotation, self.variants = \
            table_generator(**kwargs)
        return (self.mmcif, self.dssp, self.sifts,
                self.validation, self.annotation, self.variants)


Tables = Tables()