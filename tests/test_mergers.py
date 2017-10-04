#!/local/bin/python
# -*- coding: utf-8 -*-
import logging
import unittest
from os import path

try:
    from mock import patch
except ImportError:
    # python 3.5
    from unittest.mock import patch

from proteofav.config import Defaults
from proteofav.main import merge_tables
from proteofav.structures import parse_mmcif_atoms
from proteofav.sifts import _sifts_residues_regions
from proteofav.dssp import parse_dssp_residues

logging.getLogger('proteofav').setLevel(logging.CRITICAL)  # turn off logging
defaults = Defaults(path.join(path.dirname(__file__), "config.txt"))
defaults.db_cif = path.join(path.dirname(__file__), "testdata", "mmcif/")
defaults.db_sifts = path.join(path.dirname(__file__), "testdata", "sifts/")
defaults.db_dssp = path.join(path.dirname(__file__), "testdata", "dssp/")


@patch("proteofav.structures.defaults", defaults)
@patch("proteofav.dssp.defaults", defaults)
@patch("proteofav.sifts.defaults", defaults)
class TestTableMerger(unittest.TestCase):
    """Test table merging methods."""

    def setUp(self):
        """Initialize the framework for testing."""

        self.cif_to_table = parse_mmcif_atoms
        self.sifts_to_table = _sifts_residues_regions
        self.dssp_to_table = parse_dssp_residues

        self.merge_table = merge_tables

    def tearDown(self):
        """Remove testing framework by cleaning the namespace."""
        self.defaults = None
        self.merge_table = None

        self.cif_to_table = None
        self.sifts_to_table = None
        self.dssp_to_table = None

        self.merge_table = None

    def test_empty(self):
        """Test no argument cases."""
        with self.assertRaises(TypeError):
            self.merge_table(pdb_id=None)
            self.merge_table(uniprot_id=None)

    def test_camKIV_ca_atom(self):
        """Test table merger for a simple protein example."""
        data = self.merge_table(pdb_id="2w4o", chain="A")
        self.assertIsNotNone(data)
        self.assertFalse(data.empty)

        self.assertEqual(data.label_atom_id.dropna().unique()[0], 'CA', 'Other atoms than CA')

        self.assertEqual(data.PDB_dbChainId.unique()[0], 'A', 'Other chain')
        self.assertEqual(data.CHAIN.dropna().unique()[0], 'A', 'Other chain')
        self.assertEqual(data.label_asym_id.dropna().unique()[0], 'A', 'Other chain')

        self.assertEqual(data.shape[0], 349, 'wrong number of rows')
        self.assertEqual(data.AA[~data.AA.isnull()].shape[0], 278, 'wrong number of residues')
        self.assertEqual(
            data.UniProt_dbResName[~data.UniProt_dbResName.isnull()].shape[0],
            326,
            'wrong number of residues')

    def test_merge_4ibw_A_with_alt_loc(self):
        """
        Test case in a structure with alt locations."""
        data = self.merge_table(pdb_id="4ibw", chain="A")
        self.assertFalse(data.empty)

    def test_merge_3mn5_with_insertion_code(self):
        """
        Test case with insertion code
        """
        self.cif_path = path.join(path.dirname(__file__), "testdata", "mmcif/3mn5.cif")
        self.sifts_path = path.join(path.dirname(__file__), "testdata", "sifts/3mn5.xml")
        self.dssp_path = path.join(path.dirname(__file__), "testdata", "dssp/3mn5.dssp")

        self.cif = self.cif_to_table(self.cif_path)
        self.sifts = self.sifts_to_table(self.sifts_path)
        self.dssp = self.dssp_to_table(self.dssp_path)

        self.assertFalse(self.cif.empty)
        self.assertFalse(self.sifts.empty)
        self.assertFalse(self.dssp.empty)

        data = self.merge_table(pdb_id="3mn5", chain="A")
        self.assertFalse(data.empty)

    def test_merge_3fqd_A_no_pdbe_label_seq_id(self):
        self.data = self.merge_table(pdb_id='3fqd', chain='A',
                                     sequence_check='ignore')
        self.assertFalse(self.data.empty)

    def test_merge_3ehk_D_lowercased_dssp(self):
        self.data = self.merge_table(pdb_id='3ehk', chain='D')
        self.assertFalse(self.data.empty)

    @unittest.expectedFailure
    def test_merge_4v9d_BD_excessive_chains(self):
        """
        DSSP files does not have BD chain, since its chain naming only support one character.
        Although its possible to map and reference the BD chain into the mmCIF table,
        it is currently unsupported by merge_tables.
        """
        data = self.merge_table(pdb_id='4v9d', chain='BD')
        self.assertFalse(data.empty)

    def test_merge_4abo_A_DSSP_missing_first_residue(self):
        data = self.merge_table(pdb_id='4abo', chain='A')
        self.assertFalse(data.empty)

    def test_merge_4why_K_DSSP_index_as_object(self):
        data = self.merge_table(pdb_id='4why', chain='K')
        self.assertFalse(data.empty)

    def test_merge_2pm7_D_missing_residue_DSSP(self):
        data = self.merge_table(pdb_id='2pm7', chain='D')
        self.assertFalse(data.empty)

    def test_merge_4myi_A_fail(self):
        data = self.merge_table(pdb_id='2pm7', chain='D')
        self.assertFalse(data.empty)

    def test_camKIV_wrong_chain(self):
        with self.assertRaises(ValueError):
            self.merge_table(pdb_id='2w4o', chain='D')

    def test_camKIV_wrong_atom(self):
        with self.assertRaises(ValueError):
            self.merge_table(pdb_id='2w4o', chain='A', atoms=['CC'])

    def test_camKIV_atom_list(self):
        # TODO test_camIV_list_mode(self):
        data = self.merge_table(pdb_id='2w4o', chain='A', atoms=['CA', 'CB'])
        self.assertFalse(data.empty)

    def test_camKIV_atom_centroid(self):
        # TODO test_camIV_centroid_mode(self):
        data = self.merge_table(pdb_id='2w4o', chain='A', atoms='centroid')
        self.assertFalse(data.empty)

    def test_3edv_string_index(self):
        # TODO def test_sift_3edv(self): Example dbResNum is a string, therefore was not merging.
        pass

    def test_camKIV_from_uniprot_id(self):
        defaults.api_pdbe = 'http://www.ebi.ac.uk/pdbe/api/'
        with patch('proteofav.structures.defaults', defaults):
            data = self.merge_table(uniprot_id='Q16566')
            self.assertFalse(data.empty)

    def test_sequence_check_raise(self):
        # The first and the second residues in the cif file were swaped to GLY
        # so they can't be checked with dssp and sift sequences
        badcif_path = path.join(path.dirname(__file__), "testdata", "mmcif/2w4o_with_error.cif")
        baddata = self.cif_to_table(badcif_path)

        with patch("proteofav.structures.parse_mmcif_atoms", return_value=baddata):
            with self.assertRaises(ValueError):
                self.merge_table(pdb_id='2w4o')
                # data = self.merge_table(pdb_id='2w4o', sequence_check='warn') # todo try capture warn
                # self.assertFalse(data.empty)


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestTableMerger)
    unittest.TextTestRunner(verbosity=1).run(suite)
