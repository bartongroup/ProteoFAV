# -*- coding: utf-8

import logging
import numpy as np
import pandas as pd

from proteofav.structures import MMCIF
from proteofav.sifts import sifts_best, SIFTS
from proteofav.dssp import DSSP
from proteofav.variants import Variants, fetch_pdb_uniprot_mapping
from proteofav.annotation import Annotation
from proteofav.validation import Validation
from proteofav.utils import merging_down_by_key

__all__ = ['mmcif_sifts_table_merger', 'mmcif_dssp_table_merger',
           'mmcif_validation_table_merger', 'sifts_annotation_table_merger',
           'sifts_variants_table_merger', 'uniprot_vars_ensembl_vars_merger',
           '_Tables', 'Tables']
           'table_merger', 'table_generator',

log = logging.getLogger('proteofav.config')


class TableMergerError(Exception):
    pass


def mmcif_sifts_table_merger(mmcif_table, sifts_table, category='auth'):
    """
    Merge the mmCIF and SIFTS Tables.

    :param mmcif_table: mmCIF pandas DataFrame
    :param sifts_table: SIFTS pandas DataFrame
    :param category: data category to be used as precedence in _atom_site.*_*
        asym_id, seq_id and atom_id
    :return: merged pandas DataFrame
    """

    # bare minimal columns needed
    if ('{}_seq_id_full'.format(category) in mmcif_table and
                '{}_asym_id'.format(category) in mmcif_table and
                'PDB_dbResNum' in sifts_table and 'PDB_dbChainId' in sifts_table):

        # workaround for BioUnits
        if 'orig_{}_asym_id'.format(category) in mmcif_table:
            table = mmcif_table.merge(sifts_table, how='left',
                                      left_on=['{}_seq_id_full'.format(category),
                                               'orig_{}_asym_id'.format(category)],
                                      right_on=['PDB_dbResNum', 'PDB_dbChainId'])
        else:
            table = mmcif_table.merge(sifts_table, how='left',
                                      left_on=['{}_seq_id_full'.format(category),
                                               '{}_asym_id'.format(category)],
                                      right_on=['PDB_dbResNum', 'PDB_dbChainId'])

    else:
        raise TableMergerError('Not possible to merge mmCIF and SIFTS table! '
                               'Some of the necessary columns are missing...')

    log.info("Merged mmCIF and SIFTS Tables...")
    return table


def mmcif_dssp_table_merger(mmcif_table, dssp_table, category='auth'):
    """
    Merge the mmCIF and DSSP Tables.

    :param mmcif_table: mmCIF pandas DataFrame
    :param dssp_table: DSSP pandas DataFrame
    :param category: data category to be used as precedence in _atom_site.*_*
        asym_id, seq_id and atom_id
    :return: merged pandas DataFrame
    """

    # bare minimal columns needed
    if ('{}_seq_id_full'.format(category) in mmcif_table and
                '{}_asym_id'.format(category) in mmcif_table and
                'RES_FULL' in dssp_table and 'CHAIN_FULL' in dssp_table):

        # workaround for BioUnits
        if ('orig_{}_asym_id'.format(category) in mmcif_table and not
        (set(mmcif_table['{}_asym_id'.format(category)].unique()) ==
             set(dssp_table['CHAIN_FULL'].unique()))):
            table = mmcif_table.merge(dssp_table, how='left',
                                      left_on=['{}_seq_id_full'.format(category),
                                               'orig_{}_asym_id'.format(category)],
                                      right_on=['RES_FULL', 'CHAIN_FULL'])
        else:
            table = mmcif_table.merge(dssp_table, how='left',
                                      left_on=['{}_seq_id_full'.format(category),
                                               '{}_asym_id'.format(category)],
                                      right_on=['RES_FULL', 'CHAIN_FULL'])

    else:
        raise TableMergerError('Not possible to merge mmCIF and DSSP table! '
                               'Some of the necessary columns are missing...')

    log.info("Merged mmCIF and DSSP Tables...")
    return table


def mmcif_validation_table_merger(mmcif_table, validation_table, category='auth'):
    """
    Merge the mmCIF/PDB and Validation Tables.

    :param mmcif_table: mmCIF pandas DataFrame
    :param validation_table: Validation pandas DataFrame
    :param category: data category to be used as precedence in _atom_site.*_*
        asym_id, seq_id and atom_id
    :return: merged pandas DataFrame
    """

    # bare minimal columns needed
    if ('{}_seq_id_full'.format(category) in mmcif_table and
                '{}_asym_id'.format(category) in mmcif_table and
                'validation_resnum_full' in validation_table and
                'validation_chain' in validation_table):

        # workaround for BioUnits
        if 'orig_{}_asym_id'.format(category) in mmcif_table:
            table = mmcif_table.merge(validation_table, how='left',
                                      left_on=['{}_seq_id_full'.format(category),
                                               'orig_{}_asym_id'.format(category)],
                                      right_on=['validation_resnum_full',
                                                'validation_chain'])
        else:
            table = mmcif_table.merge(validation_table, how='left',
                                      left_on=['{}_seq_id_full'.format(category),
                                               '{}_asym_id'.format(category)],
                                      right_on=['validation_resnum_full',
                                                'validation_chain'])
    else:
        raise TableMergerError('Not possible to merge mmCIF and Validation table! '
                               'Some of the necessary columns are missing...')

    log.info("Merged mmCIF and Validation Tables...")
    return table


def sifts_annotation_table_merger(sifts_table, annotation_table):
    """
    Merge the SIFTS and Annotation Tables (annotation needs to be aggregated
    by residue - via `filter_annotation`).

    :param sifts_table: SIFTS pandas DataFrame
    :param annotation_table: Annotation pandas DataFrame
    :return: merged pandas DataFrame
    """

    # bare minimal columns needed
    if ('UniProt_dbAccessionId' in sifts_table and 'UniProt_dbResNum' in sifts_table and
                'accession' in annotation_table and 'site' in annotation_table):

        table = sifts_table.merge(annotation_table, how='left',
                                  left_on=['UniProt_dbAccessionId', 'UniProt_dbResNum'],
                                  right_on=['accession', 'site'])

    else:
        raise TableMergerError('Not possible to merge SIFTS and Annotation table! '
                               'Some of the necessary columns are missing...')

    log.info("Merged SIFTS and Annotation Tables...")
    return table


def sifts_variants_table_merger(sifts_table, variants_table):
    """
    Merge the mmCIF/PDB and Variants Tables.

    :param sifts_table: SIFTS pandas DataFrame
    :param variants_table: Variants pandas DataFrame
    :return: merged pandas DataFrame
    """

    # bare minimal columns needed
    if ('UniProt_dbAccessionId' in sifts_table and 'UniProt_dbResNum' in sifts_table and
                'accession' in variants_table and 'begin' in variants_table):
        variants_table['begin'] = variants_table['begin'].astype(str)

        table = sifts_table.merge(variants_table, how='left',
                                  left_on=['UniProt_dbAccessionId', 'UniProt_dbResNum'],
                                  right_on=['accession', 'begin'])

    else:
        raise TableMergerError('Not possible to merge SIFTS and Variants table! '
                               'Some of the necessary columns are missing...')

    log.info("Merged SIFTS and Variants Tables...")
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

    available = [mmcif_table, dssp_table, sifts_table,
                 validation_table, annotation_table, variants_table]
    available = [k for k in available if k is not None]
    if len(available) < 2 and (mmcif_table or sifts_table):
        raise TableMergerError(
            "At least two Tables are needed in order to merge...")

    table = None
    if mmcif_table is not None:
        if dssp_table is not None:
            mmcif_table = mmcif_dssp_table_merger(mmcif_table, dssp_table)
        if validation_table is not None:
            mmcif_table = mmcif_validation_table_merger(mmcif_table,
                                                        validation_table)
        table = mmcif_table

    if sifts_table is not None:
        if annotation_table is not None:
            sifts_table = sifts_annotation_table_merger(sifts_table,
                                                        annotation_table)
        if variants_table is not None:
            sifts_table = sifts_variants_table_merger(sifts_table,
                                                      variants_table)
        table = sifts_table

    if mmcif_table is not None and sifts_table is not None:
        table = mmcif_sifts_table_merger(mmcif_table, sifts_table)

    return table


def table_generator(uniprot_id=None, pdb_id=None, bio_unit=False,
                    sifts=True, dssp=False, validation=False, annotations=False, variants=False,
                    chains=None, res=None, sites=None, atoms=None, lines=None,
                    residue_agg=False, overwrite=False):
    """
    Simplifies the process of generating Tables and merging them.

    Automatically merge data Tables from various resources. If no pdb_id set
    `sifts_best`, which is sorted by sequence coverage. If no chain is set, uses all.

    :param uniprot_id: (str) UniProt ID
    :param pdb_id: (str) PDB ID
    :param bio_unit: boolean for using AsymUnits or BioUnits
    :param sifts: boolean
    :param dssp: boolean
    :param validation: boolean
    :param annotations: boolean
    :param variants: boolean
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
            data = sifts_best(uniprot_id)[uniprot_id]
            if data is not None:
                # uses the first structure
                pdb_id = data[0]['pdb_id']
                chains = (data[0]['chain_id'],)
            else:
                log.warning("Best structures not available from the PDBe API for %s",
                            uniprot_id)
                raise TableMergerError('Nothing to merge...')

        # mmCIF table
        mmcif_table = MMCIF.load(identifier=pdb_id, bio_unit=bio_unit,
                                 bio_unit_preferred=True, overwrite=overwrite,
                                 add_atom_altloc=True, add_res_full=True,
                                 category='auth', residue_agg=residue_agg,
                                 chains=chains, res_full=res, atoms=atoms, lines=lines)

        # SIFTS table
        if sifts:
            sifts_table = SIFTS.load(identifier=pdb_id, add_regions=True, add_dbs=False,
                                     chains=chains, res=res, uniprot=uniprot_id, site=sites)
        else:
            sifts_table = None

        # DSSP table
        if dssp:
            dssp_table = DSSP.load(identifier=pdb_id,
                                   add_ss_reduced=True, add_rsa_class=True,
                                   chains_full=chains, res=res)
        else:
            dssp_table = None

        # Validation table
        if validation:
            validation_table = Validation.load(identifier=pdb_id)
        else:
            validation_table = None

        # Annotations table
        if annotations:
            if uniprot_id:
                annotations_table = Annotation.load(identifier=uniprot_id,
                                                    annotation_agg=True)
            else:
                annotations_table = None
                annot_tables = []
                data = fetch_pdb_uniprot_mapping(identifier=pdb_id).json()[pdb_id]['UniProt']
                for uniprot_id in data:
                    annot_tables.append(Annotation.load(identifier=uniprot_id,
                                                        annotation_agg=True))
                if len(annot_tables) > 1:
                    annotations_table = pd.concat(annot_tables).reset_index(drop=True)
                elif len(annot_tables):
                    annotations_table = annot_tables[0]
        else:
            annotations_table = None

        # Variants table
        if variants:
            if uniprot_id:
                uni_vars, ens_vars = Variants.load(identifier=uniprot_id,
                                                   id_source='uniprot',
                                                   synonymous=True,
                                                   uniprot_vars=True,
                                                   ensembl_germline_vars=True,
                                                   ensembl_somatic_vars=True)
                variants_table = uniprot_vars_ensembl_vars_merger(uni_vars, ens_vars)
            else:
                variants_table = None
                vars_tables = []
                data = fetch_pdb_uniprot_mapping(identifier=pdb_id).json()[pdb_id]['UniProt']
                for uniprot_id in data:
                    uni_vars, ens_vars = Variants.load(identifier=uniprot_id,
                                                       id_source='uniprot',
                                                       synonymous=True,
                                                       uniprot_vars=True,
                                                       ensembl_germline_vars=True,
                                                       ensembl_somatic_vars=True)
                    vars_tables.append(uniprot_vars_ensembl_vars_merger(uni_vars, ens_vars))
                if len(vars_tables) > 1:
                    variants_table = pd.concat(vars_tables).reset_index(drop=True)
                elif len(vars_tables):
                    variants_table = vars_tables[0]
        else:
            variants_table = None

        return (mmcif_table, dssp_table, sifts_table,
                validation_table, annotations_table, variants_table)

    else:
        raise ValueError('No UniProt ID or PDB ID provided...')


class _Tables(object):
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

    def load(self, merge_tables=False, sequence_check='ignore', **kwargs):
        """Generates the Tables, merges and stores

        :param merge_tables: (bool) merge the Tables after generating them
        :param str sequence_check: how to handle sequence inconsistencies.
            Choose from 'raise', 'warn' or 'ignore'
        """
        self.mmcif, self.dssp, self.sifts, self.validation, self.annotation, self.variants = \
            table_generator(**kwargs)
        if not merge_tables:
            return (self.mmcif, self.dssp, self.sifts,
                    self.validation, self.annotation, self.variants)
        else:
            table = table_merger(self.mmcif, self.dssp, self.sifts,
                                 self.validation, self.annotation, self.variants)

            if sequence_check == 'raise' or sequence_check == 'warn':
                # TODO
                pass
            elif sequence_check != 'ignore':
                raise ValueError("Sequence check method '{}' not implemented."
                                 "".format(sequence_check))

            return table


Tables = _Tables()
