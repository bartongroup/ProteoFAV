# -*- coding: utf-8 -*-

import os
import sys
import logging
import unittest

try:
    from mock import patch
except ImportError:
    # python 3.5
    from unittest.mock import patch

from proteofav.config import defaults
from proteofav.structures import parse_mmcif_atoms, mmCIF
from proteofav.dssp import parse_dssp_residues, DSSP, filter_dssp
from proteofav.sifts import parse_sifts_residues, SIFTS
from proteofav.variants import Variants
from proteofav.validation import Validation
from proteofav.annotation import Annotation

from proteofav.mergers import (mmcif_dssp_table_merger,
                               mmcif_sifts_table_merger,
                               mmcif_validation_table_merger,
                               sifts_annotation_table_merger,
                               sifts_variants_table_merger,
                               uniprot_vars_ensembl_vars_merger,
                               table_merger, table_generator,
                               Tables)

root = os.path.dirname(__file__)
defaults.db_cif = os.path.join(root, "testdata", "mmcif")
defaults.db_sifts = os.path.join(root, "testdata", "sifts")
defaults.db_dssp = os.path.join(root, "testdata", "dssp")
defaults.db_validation = os.path.join(root, "testdata", "validation")
defaults.db_annotation = os.path.join(root, "testdata", "annotation")


@patch("proteofav.structures.defaults", defaults)
@patch("proteofav.dssp.defaults", defaults)
@patch("proteofav.sifts.defaults", defaults)
@patch("proteofav.validation.defaults", defaults)
@patch("proteofav.annotation.defaults", defaults)
class TestMerger(unittest.TestCase):
    """Test table merging methods."""

    def setUp(self):
        """Initialize the framework for testing."""

        self.cif_to_table = parse_mmcif_atoms
        self.sifts_to_table = parse_sifts_residues
        self.dssp_to_table = parse_dssp_residues

        self.pdbid = TestMerger.pdbid
        self.uniprotid = TestMerger.uniprotid
        self.inputbiodssp = TestMerger.inputbiodssp

        self.mmcif = TestMerger.mmcif
        self.mmcif_bio = TestMerger.mmcif_bio
        self.dssp = TestMerger.dssp
        self.dssp_bio = TestMerger.dssp_bio
        self.sifts = TestMerger.sifts
        self.validation = TestMerger.validation
        self.annotation = TestMerger.annotation
        self.uni_vars = TestMerger.uni_vars
        self.ens_vars = TestMerger.ens_vars
        self.variants = TestMerger.variants

        self.mmcif_sifts = mmcif_sifts_table_merger
        self.mmcif_dssp = mmcif_dssp_table_merger
        self.mmcif_validation = mmcif_validation_table_merger
        self.sifts_annotation = sifts_annotation_table_merger
        self.sifts_variants = sifts_variants_table_merger
        self.uni_ens_vars = uniprot_vars_ensembl_vars_merger

        self.table_merger = table_merger
        self.table_generator = table_generator
        self.Tables = Tables

    def tearDown(self):
        """Remove testing framework by cleaning the namespace."""

        self.cif_to_table = None
        self.sifts_to_table = None
        self.dssp_to_table = None

        self.pdbid = None
        self.uniprotid = None
        self.inputbiodssp = None

        self.mmcif = None
        self.mmcif_bio = None
        self.dssp = None
        self.dssp_bio = None
        self.sifts = None
        self.validation = None
        self.annotation = None
        self.uni_vars = None
        self.ens_vars = None
        self.variants = None

        self.mmcif_sifts = None
        self.mmcif_dssp = None
        self.mmcif_validation = None
        self.sifts_annotation = None
        self.sifts_variants = None
        self.uni_ens_vars = None

        self.table_merger = None
        self.table_generator = None
        self.Tables = None

    @classmethod
    def setUpClass(cls):
        # to be run only once
        super(TestMerger, cls).setUpClass()

        cls.pdbid = '2pah'
        cls.uniprotid = 'P00439'

        cls.inputbiodssp = os.path.join(root, "testdata", defaults.db_dssp,
                                        "{}_bio.dssp".format(cls.pdbid))

        cls.mmcif = mmCIF.select(identifier=cls.pdbid,
                                 add_res_full=True, atoms=('CA',))

        cls.mmcif_bio = mmCIF.select(identifier=cls.pdbid, bio_unit=True,
                                     bio_unit_preferred=True,
                                     add_res_full=True, atoms=('CA',))

        cls.dssp = DSSP.select(identifier=cls.pdbid,
                               add_rsa_class=True, add_ss_reduced=True)

        dssp_bio = DSSP.read(filename=cls.inputbiodssp)
        cls.dssp_bio = filter_dssp(dssp_bio, add_rsa_class=True, add_ss_reduced=True)

        cls.sifts = SIFTS.select(identifier=cls.pdbid, add_regions=True, add_dbs=False)

        cls.validation = Validation.select(identifier=cls.pdbid, add_res_full=True)

        cls.annotation = Annotation.select(identifier=cls.uniprotid, annotation_agg=True)

        cls.uni_vars, cls.ens_vars = Variants.select(cls.uniprotid,
                                                     id_source='uniprot',
                                                     synonymous=True,
                                                     uniprot_vars=True,
                                                     ensembl_germline_vars=True,
                                                     ensembl_somatic_vars=True)
        cls.variants = uniprot_vars_ensembl_vars_merger(cls.uni_vars, cls.ens_vars)

    @classmethod
    def tearDownClass(cls):
        cls.pdbid = None
        cls.uniprotid = None
        cls.inputbiodssp = None

        cls.mmcif = None
        cls.mmcif_bio = None
        cls.dssp = None
        cls.dssp_bio = None
        cls.sifts = None
        cls.validation = None
        cls.annotation = None
        cls.uni_vars = None
        cls.ens_vars = None
        cls.variants = None

    def test_empty(self):
        """Test no argument cases."""
        with self.assertRaises(ValueError):
            self.Tables.generate(pdb_id=None)
            self.Tables.generate(uniprot_id=None)

    def test_camKIV_ca_atom(self):
        """Test table merger for a simple protein example."""
        data = self.Tables.generate(pdb_id="2w4o", chains="A", dssp=True, atoms='CA',
                                    merge_tables=True)
        self.assertIsNotNone(data)
        self.assertFalse(data.empty)

        self.assertEqual(data.label_atom_id.dropna().unique()[0], 'CA', 'Other atoms than CA')

        self.assertEqual(data.PDB_dbChainId.unique()[0], 'A', 'Other chain')
        self.assertEqual(data.CHAIN.dropna().unique()[0], 'A', 'Other chain')
        self.assertEqual(data.label_asym_id.dropna().unique()[0], 'A', 'Other chain')

        self.assertEqual(data.shape[0], 278, 'wrong number of rows')
        self.assertEqual(data.AA[~data.AA.isnull()].shape[0], 278, 'wrong number of residues')
        self.assertEqual(data.UniProt_dbResName[~data.UniProt_dbResName.isnull()].shape[0],
                         278, 'wrong number of residues')

    def test_merge_4ibw_A_with_alt_loc(self):
        """
        Test case in a structure with alt locations."""
        data = self.Tables.generate(pdb_id="4ibw", chains="A",
                                    merge_tables=True)
        self.assertFalse(data.empty)

    def test_merge_3mn5_with_insertion_code(self):
        """
        Test case with insertion code
        """
        self.cif_path = os.path.join(os.path.dirname(__file__), "testdata",
                                     "mmcif", "3mn5.cif")
        self.sifts_path = os.path.join(os.path.dirname(__file__), "testdata",
                                       "sifts", "3mn5.xml")
        self.dssp_path = os.path.join(os.path.dirname(__file__), "testdata",
                                      "dssp", "3mn5.dssp")

        self.cif = self.cif_to_table(self.cif_path)
        self.sifts = self.sifts_to_table(self.sifts_path)
        self.dssp = self.dssp_to_table(self.dssp_path)

        self.assertFalse(self.cif.empty)
        self.assertFalse(self.sifts.empty)
        self.assertFalse(self.dssp.empty)

        data = self.Tables.generate(pdb_id="3mn5", chains="A",
                                    merge_tables=True)
        self.assertFalse(data.empty)

    def test_merge_3fqd_A_no_pdbe_label_seq_id(self):
        self.data = self.Tables.generate(pdb_id='3fqd', chains='A',
                                         merge_tables=True)
        self.assertFalse(self.data.empty)

    def test_merge_3ehk_D_lowercased_dssp(self):
        self.data = self.Tables.generate(pdb_id='3ehk', chains='D',
                                         merge_tables=True)
        self.assertFalse(self.data.empty)

    # @unittest.expectedFailure
    def test_merge_4v9d_BD_excessive_chains(self):
        """
        DSSP files does not have BD chain, since its chain naming only support one character.
        Although its possible to map and reference the BD chain into the mmCIF table,
        it is currently unsupported by merge_tables.
        """
        data = self.Tables.generate(pdb_id='4v9d', chains='BD',
                                    merge_tables=True)
        self.assertFalse(data.empty)

    def test_merge_4abo_A_DSSP_missing_first_residue(self):
        data = self.Tables.generate(pdb_id='4abo', chains='A',
                                    merge_tables=True)
        self.assertFalse(data.empty)

    def test_merge_4why_K_DSSP_index_as_object(self):
        data = self.Tables.generate(pdb_id='4why', chains='K',
                                    merge_tables=True)
        self.assertFalse(data.empty)

    def test_merge_2pm7_D_missing_residue_DSSP(self):
        data = self.Tables.generate(pdb_id='2pm7', chains='D',
                                    merge_tables=True)
        self.assertFalse(data.empty)

    def test_merge_4myi_A_fail(self):
        data = self.Tables.generate(pdb_id='2pm7', chains='D',
                                    merge_tables=True)
        self.assertFalse(data.empty)

    def test_camKIV_wrong_chain(self):
        with self.assertRaises(ValueError):
            self.Tables.generate(pdb_id='2w4o', chains='D',
                                 merge_tables=True)

    def test_camKIV_wrong_atom(self):
        with self.assertRaises(ValueError):
            self.Tables.generate(pdb_id='2w4o', chains='A', atoms=['CC'],
                                 merge_tables=True)

    def test_camKIV_atom_list(self):
        # TODO test_camIV_list_mode(self):
        data = self.Tables.generate(pdb_id='2w4o', chains='A', atoms=['CA', 'CB'],
                                    merge_tables=True)
        self.assertFalse(data.empty)

    def test_camKIV_atom_centroid(self):
        # TODO test_camIV_centroid_mode(self):
        data = self.Tables.generate(pdb_id='2w4o', chains='A', atoms='centroid',
                                    merge_tables=True)
        self.assertFalse(data.empty)

    def test_3edv_string_index(self):
        # TODO def test_sift_3edv(self): Example dbResNum is a string, therefore was not merging.
        pass

    def test_camKIV_from_uniprot_id(self):
        data = self.Tables.generate(uniprot_id='Q16566', merge_tables=True)
        self.assertFalse(data.empty)

    def test_sequence_check_raise(self):
        # The first and the second residues in the cif file were swapped to GLY
        # so they can't be checked with dssp and sifts sequences
        badcif_path = os.path.join(os.path.dirname(__file__), "testdata",
                                   "mmcif", "2w4o_with_error.cif")
        baddata = self.cif_to_table(badcif_path)

        with patch("proteofav.structures.parse_mmcif_atoms", return_value=baddata):
            # with self.assertRaises(ValueError):
            data = self.Tables.generate(pdb_id='2w4o', merge_tables=True)
            self.assertFalse(data.empty)

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

    def test_mmcif_validation_merger(self):
        table = self.mmcif_validation(self.mmcif, self.validation)
        self.assertIn('auth_asym_id', table)
        self.assertIn('validation_rsr', table)
        self.assertEqual(table.loc[0, 'auth_asym_id'], 'A')
        self.assertEqual(table.loc[0, 'validation_rsr'], 0.242)
        self.assertEqual(table.loc[0, 'auth_seq_id_full'], '118')
        self.assertEqual(table.loc[0, 'validation_resnum_full'], '118')

    def test_sifts_annotation_merger(self):
        table = self.sifts_annotation(self.sifts, self.annotation)
        self.assertIn('UniProt_dbResNum', table)
        self.assertIn('site', table)
        self.assertEqual(table.loc[0, 'UniProt_dbResNum'], '118')
        self.assertEqual(table.loc[3, 'UniProt_dbResNum'], '121')
        self.assertEqual(table.loc[3, 'site'], '121')
        self.assertIn('Natural variant:', table.loc[3, 'annotation'])

    def test_sifts_variants_merger(self):
        table = self.sifts_variants(self.sifts, self.variants)
        self.assertIn('UniProt_dbAccessionId', table)
        self.assertIn('UniProt_dbResNum', table)
        self.assertIn('xrefs_id', table)
        self.assertIn('accession', table)
        self.assertIn('begin', table)
        self.assertEqual(table.loc[0, 'UniProt_dbResNum'], '118')
        self.assertEqual(table.loc[0, 'UniProt_dbAccessionId'], 'P00439')
        self.assertEqual(table.loc[0, 'accession'], 'P00439')
        self.assertEqual(table.loc[0, 'xrefs_id'], 'rs776442422')
        self.assertEqual(table.loc[0, 'siftScore'], 0.14)

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

    def test_table_merger(self):
        table = self.Tables.merge(self.mmcif, self.dssp, self.sifts,
                                  self.validation, self.annotation, self.variants)
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
        table = self.table_merger(self.mmcif, self.dssp, self.sifts,
                                  self.validation, self.annotation, self.variants)
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
        table = self.table_merger(self.mmcif_bio, self.dssp_bio, self.sifts,
                                  self.validation, self.annotation, self.variants)
        # Chain level
        self.assertIn('label_asym_id', table)
        self.assertIn('CHAIN_FULL', table)
        self.assertIn('PDB_entityId', table)
        # Res level
        self.assertIn('label_seq_id_full', table)
        self.assertIn('RES', table)
        self.assertIn('PDB_dbResNum', table)
        # values
        self.assertEqual('CA', table.loc[407, 'label_atom_id'])
        self.assertEqual('AA', table.loc[407, 'label_asym_id'])
        self.assertEqual('118', table.loc[407, 'RES'])
        self.assertEqual('VAL', table.loc[407, 'PDB_dbResName'])
        self.assertEqual('V', table.loc[407, 'UniProt_dbResName'])

    def test_table_generator(self):
        mmcif_table, dssp_table, sifts_table, valid_table, annot_table, vars_table = \
            self.table_generator(uniprot_id=None, pdb_id=self.pdbid, bio_unit=True,
                                 sifts=True, dssp=True, validation=False, annotations=False, variants=False,
                                 chains=None, res=None, sites=None, atoms=('CA',), lines=None,
                                 residue_agg=False, overwrite=False)

        table = self.table_merger(mmcif_table, dssp_table, sifts_table,
                                  valid_table, annot_table, vars_table)
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
        mmcif_table, dssp_table, sifts_table, valid_table, annot_table, vars_table = \
            self.table_generator(uniprot_id=None, pdb_id=self.pdbid, bio_unit=True,
                                 sifts=True, dssp=True, variants=False, annotations=False,
                                 chains=None, res=None, sites=None, atoms=('CA',), lines=None,
                                 residue_agg=False, overwrite=False)

        table = self.table_merger(mmcif_table, dssp_table, sifts_table,
                                  valid_table, annot_table, vars_table)
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
        # self.assertEqual('118', table.loc[329, 'RES'])
        self.assertEqual('VAL', table.loc[329, 'PDB_dbResName'])
        self.assertEqual('V', table.loc[329, 'UniProt_dbResName'])

    def test_table_merger_bio(self):
        self.Tables.generate(pdb_id=self.pdbid, atoms=('CA',), bio_unit=True,
                             dssp=True, sifts=True)
        table = self.Tables.merge()
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
        # self.assertEqual('118', table.loc[329, 'RES'])
        self.assertEqual('VAL', table.loc[329, 'PDB_dbResName'])
        self.assertEqual('V', table.loc[329, 'UniProt_dbResName'])

    def test_table_generator_full(self):
        mmcif_table, dssp_table, sifts_table, valid_table, annot_table, vars_table = \
            self.table_generator(uniprot_id=None, pdb_id=self.pdbid, bio_unit=True,
                                 sifts=True, dssp=True, validation=True, annotations=True, variants=True,
                                 chains=None, res=None, sites=None, atoms=('CA',), lines=None,
                                 residue_agg=False, overwrite=False)

        table = self.table_merger(mmcif_table, dssp_table, sifts_table,
                                  valid_table, annot_table, vars_table)

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
        self.assertEqual('AA', table.loc[407, 'label_asym_id'])
        self.assertEqual('118', table.loc[0, 'RES'])
        self.assertEqual('118', table.loc[407, 'RES'])
        self.assertEqual('VAL', table.loc[0, 'PDB_dbResName'])
        self.assertEqual('V', table.loc[0, 'UniProt_dbResName'])
        # validation
        self.assertEqual('118', table.loc[0, 'validation_resnum_full'])
        self.assertEqual('118', table.loc[407, 'validation_resnum_full'])
        self.assertEqual('A', table.loc[0, 'validation_chain'])
        self.assertEqual('A', table.loc[407, 'validation_chain'])
        # annotations
        self.assertIn('Natural variant:', table.loc[4, 'annotation'])
        # variants
        self.assertEqual('large_scale_study', table.loc[0, 'sourceType'])
        self.assertEqual('ENSP00000448059', table.loc[0, 'translation'])
        self.assertEqual('rs776442422', table.loc[0, 'xrefs_id'])
        self.assertEqual('ExAC', table.loc[0, 'xrefs_name'])
        self.assertEqual(0.14, table.loc[0, 'siftScore'])


if __name__ == '__main__':
    logging.basicConfig(stream=sys.stderr)
    logging.getLogger("proteofav.config").setLevel(logging.CRITICAL)
    suite = unittest.TestLoader().loadTestsFromTestCase(TestMerger)
    unittest.TextTestRunner(verbosity=2).run(suite)
