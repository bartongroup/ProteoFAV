# -*- coding: utf-8 -*-


import os
import sys
import logging
import unittest
from unittest.mock import patch

from proteofav.dssp import DSSP
from proteofav.sifts import SIFTS
from proteofav.structures import mmCIF

from proteofav.fetchers import (fetch_uniprot_variants_ebi,
                                fetch_ensembl_transcript_variants)

from proteofav.variants import (Variants, flatten_uniprot_variants_ebi,
                                flatten_ensembl_variants, uniprot_vars_ensembl_vars_merger)

from proteofav.mergers import (Tables, table_merger,
                               mmcif_dssp_table_merger, mmcif_sifts_table_merger,
                               dssp_sifts_table_merger, table_generator)

from proteofav.config import defaults as config

root = os.path.abspath(os.path.dirname(__file__))
config.db_root = os.path.join(root, "testdata")


@patch("proteofav.config.defaults.db_root", config.db_root)
class TestMerger(unittest.TestCase):
    """Test the Merger methods."""

    def setUp(self):
        """Initialize the framework for testing."""

        self.pdbid = TestMerger.pdbid
        self.inputcif = TestMerger.inputcif
        self.inputdssp = TestMerger.inputdssp
        self.inputbiocif = TestMerger.inputbiocif
        self.inputbiodssp = TestMerger.inputbiodssp
        self.inputsifts = TestMerger.inputsifts

        self.mmcif = TestMerger.mmcif
        self.mmcif_bio = TestMerger.mmcif_bio
        self.dssp = TestMerger.dssp
        self.dssp_bio = TestMerger.dssp_bio
        self.sifts = TestMerger.sifts

        self.mmcif_sifts = mmcif_sifts_table_merger
        self.mmcif_dssp = mmcif_dssp_table_merger
        self.dssp_sifts = dssp_sifts_table_merger
        self.tables = Tables
        self.table_merger = table_merger
        self.generator = table_generator

        self.uni_ens_vars = uniprot_vars_ensembl_vars_merger

        logging.disable(logging.DEBUG)

    def tearDown(self):
        """Remove testing framework."""

        self.pdbid = None
        self.inputcif = None
        self.inputdssp = None
        self.inputbiocif = None
        self.inputbiodssp = None
        self.inputsifts = None

        self.mmcif = None
        self.mmcif_bio = None
        self.dssp = None
        self.dssp_bio = None
        self.sifts = None

        self.mmcif_sifts = None
        self.mmcif_dssp = None
        self.dssp_sifts = None
        self.dssp_dssp = None
        self.contacts_mmcif = None
        self.tables = None
        self.table_merger = None
        self.generator = None

        self.uni_ens_vars = None

        logging.disable(logging.NOTSET)

    @classmethod
    def setUpClass(cls):
        # to be run only once
        super(TestMerger, cls).setUpClass()

        cls.pdbid = '2pah'
        cls.inputcif = os.path.join(root, "testdata", config.db_mmcif,
                                    "{}.cif".format(cls.pdbid))
        cls.inputdssp = os.path.join(root, "testdata", config.db_dssp,
                                     "{}.dssp".format(cls.pdbid))

        cls.inputbiocif = os.path.join(root, "testdata", config.db_mmcif,
                                       "{}_bio.cif".format(cls.pdbid))
        cls.inputbiodssp = os.path.join(root, "testdata", config.db_dssp,
                                        "{}_bio.dssp".format(cls.pdbid))

        cls.inputsifts = os.path.join(root, "testdata", config.db_sifts,
                                      "{}.xml".format(cls.pdbid))

        d = mmCIF(filename=cls.inputcif)
        cls.mmcif = d.read(add_res_full=True, atoms=('CA',))

        d = mmCIF(filename=cls.inputbiocif)
        cls.mmcif_bio = d.read(add_res_full=True, atoms=('CA',))

        d = DSSP(filename=cls.inputdssp)
        cls.dssp = d.read(add_rsa_class=True, add_ss_reduced=True)
        d = DSSP(filename=cls.inputbiodssp)
        cls.dssp_bio = d.read(add_rsa_class=True, add_ss_reduced=True)

        d = SIFTS(filename=cls.inputsifts)
        cls.sifts = d.read(add_regions=True, add_dbs=False)

        cls.uniprotid = 'P40227'
        v = Variants(cls.uniprotid, uniprot=True)
        r = fetch_uniprot_variants_ebi(v.uniprot_id)
        if r is not None:
            cls.uni_vars = flatten_uniprot_variants_ebi(r)

        r = fetch_ensembl_transcript_variants(v.ensembl_id)
        if r is not None:
            cls.ens_vars = flatten_ensembl_variants(r)

    @classmethod
    def tearDownClass(cls):

        cls.pdbid = None
        cls.inputcif = None
        cls.inputdssp = None
        cls.inputbiocif = None
        cls.inputbiodssp = None
        cls.inputsifts = None

        cls.mmcif = None
        cls.mmcif_bio = None
        cls.dssp = None
        cls.dssp_bio = None
        cls.sifts = None

        cls.uniprotid = None
        cls.uni_vars = None
        cls.ens_vars = None

    def test_mmcif_dssp_merger(self):
        table = self.mmcif_dssp(self.mmcif, self.dssp)
        # Chain level
        self.assertIn('label_asym_id', table)
        self.assertIn('CHAIN_FULL', table)
        self.assertNotIn('PDB_entityId', table)
        # Res level
        self.assertIn('label_seq_id_full', table)
        self.assertIn('RES', table)
        self.assertNotIn('PDB_dbResNum', table)
        # values
        self.assertEqual('CA', table.loc[0, 'label_atom_id'])
        self.assertEqual('A', table.loc[0, 'label_asym_id'])
        self.assertEqual('118', table.loc[0, 'RES'])
        self.assertEqual('V', table.loc[0, 'AA'])

    def test_mmcif_dssp_bio_merger(self):
        table = self.mmcif_dssp(self.mmcif_bio, self.dssp_bio)
        # Chain level
        self.assertIn('label_asym_id', table)
        self.assertIn('CHAIN_FULL', table)
        self.assertNotIn('PDB_entityId', table)
        # Res level
        self.assertIn('label_seq_id_full', table)
        self.assertIn('RES', table)
        self.assertNotIn('PDB_dbResNum', table)
        # values
        self.assertEqual('CA', table.loc[329, 'label_atom_id'])
        self.assertEqual('AA', table.loc[329, 'label_asym_id'])
        self.assertEqual('118', table.loc[329, 'RES'])
        self.assertEqual('V', table.loc[329, 'AA'])

    def test_mmcif_sifts_merger(self):
        table = self.mmcif_sifts(self.mmcif, self.sifts)
        # Chain level
        self.assertIn('label_asym_id', table)
        self.assertNotIn('CHAIN_FULL', table)
        self.assertIn('PDB_entityId', table)
        # Res level
        self.assertIn('label_seq_id_full', table)
        self.assertNotIn('RES', table)
        self.assertIn('PDB_dbResNum', table)
        # values
        self.assertEqual('CA', table.loc[0, 'label_atom_id'])
        self.assertEqual('A', table.loc[0, 'label_asym_id'])
        self.assertEqual('118', table.loc[0, 'PDB_dbResNum'])
        self.assertEqual('VAL', table.loc[0, 'PDB_dbResName'])

    def test_mmcif_sifts_bio_merger(self):
        table = self.mmcif_sifts(self.mmcif_bio, self.sifts)
        # Chain level
        self.assertIn('label_asym_id', table)
        self.assertNotIn('CHAIN_FULL', table)
        self.assertIn('PDB_entityId', table)
        # Res level
        self.assertIn('label_seq_id_full', table)
        self.assertNotIn('RES', table)
        self.assertIn('PDB_dbResNum', table)
        # values
        self.assertEqual('CA', table.loc[329, 'label_atom_id'])
        self.assertEqual('AA', table.loc[329, 'label_asym_id'])
        self.assertEqual('118', table.loc[329, 'PDB_dbResNum'])
        self.assertEqual('VAL', table.loc[329, 'PDB_dbResName'])

    def test_dssp_sifts_merger(self):
        table = self.dssp_sifts(self.dssp, self.sifts)
        # Chain level
        self.assertNotIn('label_asym_id', table)
        self.assertIn('CHAIN_FULL', table)
        self.assertIn('PDB_entityId', table)
        # Res level
        self.assertNotIn('label_seq_id_full', table)
        self.assertIn('RES', table)
        self.assertIn('PDB_dbResNum', table)
        # values
        self.assertEqual('A', table.loc[0, 'PDB_entityId'])
        self.assertEqual('118', table.loc[0, 'RES'])
        self.assertEqual('VAL', table.loc[0, 'PDB_dbResName'])

    def test_dssp_sifts_bio_merger(self):
        table = self.dssp_sifts(self.dssp, self.sifts)
        # Chain level
        self.assertNotIn('label_asym_id', table)
        self.assertIn('CHAIN_FULL', table)
        self.assertIn('PDB_entityId', table)
        # Res level
        self.assertNotIn('label_seq_id_full', table)
        self.assertIn('RES', table)
        self.assertIn('PDB_dbResNum', table)
        # values
        self.assertEqual('B', table.loc[329, 'PDB_entityId'])
        self.assertEqual('118', table.loc[329, 'RES'])
        self.assertEqual('VAL', table.loc[329, 'PDB_dbResName'])

    def test_table_merger(self):
        table = self.tables(self.mmcif, self.dssp, self.sifts).merge()
        # Chain level
        self.assertIn('label_asym_id', table)
        self.assertIn('CHAIN_FULL', table)
        self.assertIn('PDB_entityId', table)
        # Res level
        self.assertIn('label_seq_id_full', table)
        self.assertIn('RES', table)
        self.assertIn('PDB_dbResNum', table)
        # values
        self.assertEqual('CA', table.loc[0, 'label_atom_id'])
        self.assertEqual('A', table.loc[0, 'label_asym_id'])
        self.assertEqual('118', table.loc[0, 'RES'])
        self.assertEqual('VAL', table.loc[0, 'PDB_dbResName'])

    def test_table_merger_method(self):
        table = self.table_merger(self.mmcif, self.dssp, self.sifts)
        # Chain level
        self.assertIn('label_asym_id', table)
        self.assertIn('CHAIN_FULL', table)
        self.assertIn('PDB_entityId', table)
        # Res level
        self.assertIn('label_seq_id_full', table)
        self.assertIn('RES', table)
        self.assertIn('PDB_dbResNum', table)
        # values
        self.assertEqual('CA', table.loc[0, 'label_atom_id'])
        self.assertEqual('A', table.loc[0, 'label_asym_id'])
        self.assertEqual('118', table.loc[0, 'RES'])
        self.assertEqual('VAL', table.loc[0, 'PDB_dbResName'])
        self.assertEqual('V', table.loc[0, 'UniProt_dbResName'])

    def test_table_merger_method_bio(self):
        table = self.table_merger(self.mmcif_bio, self.dssp_bio, self.sifts)
        # Chain level
        self.assertIn('label_asym_id', table)
        self.assertIn('CHAIN_FULL', table)
        self.assertIn('PDB_entityId', table)
        # Res level
        self.assertIn('label_seq_id_full', table)
        self.assertIn('RES', table)
        self.assertIn('PDB_dbResNum', table)
        # values
        self.assertEqual('CA', table.loc[329, 'label_atom_id'])
        self.assertEqual('AA', table.loc[329, 'label_asym_id'])
        self.assertEqual('118', table.loc[329, 'RES'])
        self.assertEqual('VAL', table.loc[329, 'PDB_dbResName'])
        self.assertEqual('V', table.loc[329, 'UniProt_dbResName'])

    def test_table_generator(self):
        mmcif_table, dssp_table, sifts_table, variants_table, annot_table = \
            self.generator(uniprot_id=None, pdb_id=self.pdbid, bio_unit=False,
                           sifts=True, dssp=True, variants=False, annotations=False,
                           chains=None, res=None, sites=None, atoms=('CA',), lines=None,
                           residue_agg=False, overwrite=False)

        table = self.table_merger(mmcif_table, dssp_table, sifts_table)
        # Chain level
        self.assertIn('label_asym_id', table)
        self.assertIn('CHAIN_FULL', table)
        self.assertIn('PDB_entityId', table)
        # Res level
        self.assertIn('label_seq_id_full', table)
        self.assertIn('RES', table)
        self.assertIn('PDB_dbResNum', table)
        # values
        self.assertEqual('CA', table.loc[0, 'label_atom_id'])
        self.assertEqual('A', table.loc[0, 'label_asym_id'])
        self.assertEqual('118', table.loc[0, 'RES'])
        self.assertEqual('VAL', table.loc[0, 'PDB_dbResName'])
        self.assertEqual('V', table.loc[0, 'UniProt_dbResName'])

    def test_table_generator_bio(self):
        mmcif_table, dssp_table, sifts_table, variants_table, annot_table = \
            self.generator(uniprot_id=None, pdb_id=self.pdbid, bio_unit=True,
                           sifts=True, dssp=True, variants=False, annotations=False,
                           chains=None, res=None, sites=None, atoms=('CA',), lines=None,
                           residue_agg=False, overwrite=False)

        table = self.table_merger(mmcif_table, dssp_table, sifts_table)
        # Chain level
        self.assertIn('label_asym_id', table)
        self.assertIn('CHAIN_FULL', table)
        self.assertIn('PDB_entityId', table)
        # Res level
        self.assertIn('label_seq_id_full', table)
        self.assertIn('RES', table)
        self.assertIn('PDB_dbResNum', table)
        # values
        self.assertEqual('CA', table.loc[329, 'label_atom_id'])
        self.assertEqual('AA', table.loc[329, 'label_asym_id'])
        self.assertEqual('118', table.loc[329, 'RES'])
        self.assertEqual('VAL', table.loc[329, 'PDB_dbResName'])
        self.assertEqual('V', table.loc[329, 'UniProt_dbResName'])

    def test_table_merger_bio(self):
        s = self.tables()
        s.generate(pdb_id=self.pdbid, atoms=('CA',), bio_unit=True,
                   dssp=True, sifts=True)
        table = s.merge()
        # Chain level
        self.assertIn('label_asym_id', table)
        self.assertIn('CHAIN_FULL', table)
        self.assertIn('PDB_entityId', table)
        # Res level
        self.assertIn('label_seq_id_full', table)
        self.assertIn('RES', table)
        self.assertIn('PDB_dbResNum', table)
        # values
        self.assertEqual('CA', table.loc[329, 'label_atom_id'])
        self.assertEqual('AA', table.loc[329, 'label_asym_id'])
        self.assertEqual('118', table.loc[329, 'RES'])
        self.assertEqual('VAL', table.loc[329, 'PDB_dbResName'])
        self.assertEqual('V', table.loc[329, 'UniProt_dbResName'])

    def test_uni_ens_vars_merger(self):
        table = self.uni_ens_vars(self.uni_vars, self.ens_vars)
        # UniProt
        self.assertNotIn('translation', list(self.uni_vars))
        self.assertNotIn('allele', list(self.uni_vars))
        self.assertIn('taxid', list(self.uni_vars))
        self.assertIn('description', list(self.uni_vars))
        # Ensembl
        self.assertIn('translation', list(self.ens_vars))
        self.assertIn('allele', list(self.ens_vars))
        self.assertNotIn('taxid', list(self.ens_vars))
        self.assertNotIn('description', list(self.ens_vars))
        # Merged Table
        self.assertIn('translation', list(table))
        self.assertIn('allele', list(table))
        self.assertIn('taxid', list(table))
        self.assertIn('description', list(table))
        self.assertIn('begin', list(table))
        self.assertIn('end', list(table))


if __name__ == '__main__':
    logging.basicConfig(stream=sys.stderr)
    logging.getLogger("proteofav").setLevel(logging.DEBUG)
    suite = unittest.TestLoader().loadTestsFromTestCase(TestMerger)
    unittest.TextTestRunner(verbosity=2).run(suite)
