# -*- coding: utf-8 -*-

import os
import logging

from proteofav.annotations import Annotations
from proteofav.config import defaults as config
from proteofav.dssp import DSSP
from proteofav.fetchers import fetch_best_structures_pdbe
from proteofav.sifts import SIFTS
from proteofav.structures import mmCIF
from proteofav.variants import Variants
from proteofav.utils import TableMergerError

log = logging.getLogger("proteofav")


def mmcif_sifts_table_merger(mmcif_table, sifts_table):
    """
    Merge the mmCIF and SIFTS tables.

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

    log.info("Merged mmCIF and SIFTS tables...")
    return table


def mmcif_dssp_table_merger(mmcif_table, dssp_table, pro_format=False):
    """
    Merge the mmCIF and DSSP tables.

    :param mmcif_table: mmCIF pandas DataFrame
    :param dssp_table: DSSP pandas DataFrame
    :param pro_format: boolean
    :return: merged pandas DataFrame
    """

    # bare minimal columns needed
    if ('new_seq_id' in mmcif_table and 'new_asym_id' in mmcif_table and
                'RES' in dssp_table and 'CHAIN' in dssp_table and pro_format):

        table = mmcif_table.merge(dssp_table, how='left',
                                  left_on=['new_seq_id', 'new_asym_id'],
                                  right_on=['RES', 'CHAIN_FULL'])
    elif ('auth_seq_id_full' in mmcif_table and 'label_asym_id' in mmcif_table and
                  'RES' in dssp_table and 'CHAIN_FULL' in dssp_table):

        table = mmcif_table.merge(dssp_table, how='left',
                                  left_on=['auth_seq_id_full', 'label_asym_id'],
                                  right_on=['RES', 'CHAIN_FULL'])
    else:
        raise TableMergerError('Not possible to merge mmCIF and DSSP table! '
                               'Some of the necessary columns are missing...')

    log.info("Merged mmCIF and DSSP tables...")
    return table


def dssp_sifts_table_merger(dssp_table, sifts_table):
    """
    Merge the DSSP and SIFTS tables.

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

    log.info("Merged DSSP and SIFTS tables...")
    return table


def table_merger(mmcif_table=None, dssp_table=None, sifts_table=None):
    """
    Merges the tables provided using appropriate columns.

    :param mmcif_table: (optional) mmCIF pandas DataFrame
    :param dssp_table: (optional) DSSP pandas DataFrame
    :param sifts_table: (optional) SIFTS pandas DataFrame
    :return: merged pandas DataFrame
    """

    available = [mmcif_table, dssp_table, sifts_table]
    available = [k for k in available if k is not None]
    if len(available) < 2:
        raise TableMergerError(
            "At least two tables are needed in order to merge...")

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
    Simplifies the process of generating tables and merging them.

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
    :returns: mmcif, sifts, dssp, variants and annotations tables
    """

    # generates the tables for the uniprot_id, pdb_id
    # chains, res, sites are only processed after
    if uniprot_id or pdb_id:
        # get pdb_id if only uniprot_id is provided
        if uniprot_id and not pdb_id:
            # fetch the best structures from the PDBe API
            data = fetch_best_structures_pdbe(uniprot_id)
            if data is not None:
                # uses the first structure
                pdb_id = data[0]['pdb_id']
                chains = (data[0]['chain_id'],)
            else:
                log.info("Best structures not available from the PDBe API for %s",
                         uniprot_id)
                raise TableMergerError('Nothing to merge...')

        # mmCIF table
        if bio_unit:
            inputcif = os.path.join(config.db_root, config.db_mmcif,
                                    "{}_bio.cif".format(pdb_id))
        else:
            inputcif = os.path.join(config.db_root, config.db_mmcif,
                                    "{}.cif".format(pdb_id))

        r = mmCIF(identifier=pdb_id, filename=inputcif)

        if not os.path.exists(inputcif) or overwrite:
            r.download(bio_unit=bio_unit, bio_unit_preferred=True, overwrite=overwrite)

        mmcif_table = r.read(add_atom_altloc=True,
                             category='auth', residue_agg=residue_agg,
                             chains=chains, res_full=res, atoms=atoms, lines=lines)

        # SIFTS table
        if sifts:
            inputsifts = os.path.join(config.db_root, config.db_sifts,
                                      "{}.xml".format(pdb_id))

            r = SIFTS(identifier=pdb_id, filename=inputsifts)

            if not os.path.exists(inputsifts) or overwrite:
                r.download(overwrite=overwrite)

            sifts_table = r.read(add_regions=True, add_dbs=False, chains=chains,
                                 res=res, uniprot=uniprot_id, site=sites)
        else:
            sifts_table = None

        # DSSP table
        if dssp:
            if bio_unit:
                outputdssp = os.path.join(config.db_root, config.db_dssp,
                                          "{}_bio.dssp".format(pdb_id))
            else:
                outputdssp = os.path.join(config.db_root, config.db_dssp,
                                          "{}.dssp".format(pdb_id))

            r = DSSP(identifier=pdb_id, filename=outputdssp)
            if not os.path.isfile(outputdssp) or overwrite:
                r.download(overwrite=overwrite)

            dssp_table = r.read(add_ss_reduced=True, add_rsa_class=True,
                                chains_full=chains, res=res)
        else:
            dssp_table = None

        # Variants table
        if variants:
            r = Variants(identifier=uniprot_id, uniprot=True)
            variants_table = r.fetch(synonymous=False,
                                     uniprot_vars=True,
                                     ensembl_transcript_vars=True,
                                     ensembl_somatic_vars=True)
        else:
            variants_table = None

        # Annotations table
        if annotations:
            r = Annotations(identifier=uniprot_id)
            annotations_table = r.fetch()
        else:
            annotations_table = None

        return mmcif_table, dssp_table, sifts_table, variants_table, annotations_table

    else:
        raise ValueError('No UniProt ID or PDB ID provided...')


class Tables(object):
    def __init__(self, pdbx_table=None, dssp_table=None, sifts_table=None,
                 variants_table=None, annotations_table=None):
        self.table = None
        self.mmcif_table = pdbx_table
        self.dssp_table = dssp_table
        self.sifts_table = sifts_table
        self.variants_table = variants_table
        self.annotations_table = annotations_table

    def merge(self):
        """Merges the provided tables and stores"""

        self.table = table_merger(self.mmcif_table,
                                  self.dssp_table,
                                  self.sifts_table)
        return self.table

    def generate(self, **kwargs):
        """Generates the tables, merges and stores"""

        mmcif, dssp, sifts, variants, annotations = table_generator(**kwargs)
        self.mmcif_table = mmcif
        self.dssp_table = dssp
        self.sifts_table = sifts
        self.variants_table = variants
        self.annotations_table = annotations

        return (self.mmcif_table, self.dssp_table,
                self.sifts_table, self.variants_table, self.annotations_table)
