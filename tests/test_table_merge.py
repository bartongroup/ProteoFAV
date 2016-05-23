#!/local/bin/python
# -*- coding: utf-8 -*-


import logging
import unittest
from os import path

from mock import patch

from proteofav.config import Defaults
from proteofav.main import merge_tables
from proteofav.structures import _dssp, _sifts_residues, _mmcif_atom

log = logging.getLogger(__name__)

defaults = Defaults("config.txt")


@patch("proteofav.structures.defaults", defaults)
class TestTableMerger(unittest.TestCase):
    """Test table merging methodsthe DSSP parser methods."""

    def setUp(self):
        """Initialize the framework for testing."""

        self.cif_to_table = _mmcif_atom
        self.sifts_to_table = _sifts_residues
        self.dssp_to_table = _dssp

        self.merge_table = merge_tables

    def tearDown(self):
        """Remove testing framework by cleaning the namespace."""
        self.defaults = None
        self.merge_table = None

        self.cif_to_table = None
        self.sifts_to_table = None
        self.dssp_to_table = None

        self.merge_table = None

    def test_camIV_ca_mode(self):
        """
        Test table merger for a protein example.
        :return:
        """
        data = self.merge_table(pdb_id="2w4o", chain="A")
        self.assertIsNotNone(data)

        self.cif_path = path.join(path.dirname(__file__), "CIF/2w4o.cif")
        self.sifts_path = path.join(path.dirname(__file__), "SIFTS/2w4o.xml")
        self.dssp_path = path.join(path.dirname(__file__), "DSSP/2w4o.dssp")
        self.cif = self.cif_to_table(self.cif_path)
        self.sifts = self.sifts_to_table(self.sifts_path)
        self.dssp = self.dssp_to_table(self.dssp_path)

        self.assertFalse(self.cif.empty)
        self.assertFalse(self.sifts.empty)
        self.assertFalse(self.dssp.empty)

        # number of rows (residues) is determined by the number of residues in
        # cif file
        chain = self.sifts
        n_residues = chain.shape[0]

        self.assertEqual(data.shape[0], n_residues)

        # number of columns is given buy sum of cols in cif, dssp and sifts
        groupby = {'label_comp_id': 'unique', 'label_atom_id': 'unique',
                   'label_asym_id': 'unique', 'Cartn_x': 'unique',
                   'Cartn_y': 'unique', 'Cartn_z': 'unique',
                   'occupancy': 'unique', 'B_iso_or_equiv': 'unique',
                   'id': 'unique'}
        self.cif = self.cif.groupby('auth_seq_id').agg(groupby)

    def test_merge_4ibw_A_with_alt_loc(self):
        """
        Test case in a structure with alt locations.
        """
        data = self.merge_table(pdb_id="4ibw", chain="A")
        self.assertFalse(data.empty)

    def test_merge_3mn5_with_insertion_code(self):
        """
        Test case with insertion code
        """
        self.cif_path = path.join(path.dirname(__file__), "CIF/3mn5.cif")
        self.sifts_path = path.join(path.dirname(__file__), "SIFTS/3mn5.xml")
        self.dssp_path = path.join(path.dirname(__file__), "DSSP/3mn5.dssp")

        self.cif = self.cif_to_table(self.cif_path)
        self.sifts = self.sifts_to_table(self.sifts_path)
        self.dssp = self.dssp_to_table(self.dssp_path)

        self.assertFalse(self.cif.empty)
        self.assertFalse(self.sifts.empty)
        self.assertFalse(self.dssp.empty)

        data = self.merge_table(pdb_id="3mn5", chain="A")
        self.assertFalse(data.empty)

    def test_merge_3fqd_A_no_pdbe_label_seq_id(self):
        self.data = self.merge_table(pdb_id='3fqd', chain='A')
        self.assertFalse(self.data.empty)

    def test_merge_3ehk_D_lowercased_dssp(self):
        self.data = self.merge_table(pdb_id='3ehk', chain='D')
        self.assertFalse(self.data.empty)

    @unittest.expectedFailure
    def test_merge_4v9d_BD_excessive_chains(self):
        """
        DSSP files does not have BD chain. What can be done is mapping the Cif sequence to the
        order of chain in DSSP but merge_tables won't support that by default.
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

        # TODO assert fails graciously without data
        # TODO sequence_check='raise' assert raise
        # TODO sequence_check='raise' assert warn
        # TODO test_camIV_list_mode(self):
        # TODO test_camIV_centroid_mode(self):
        # TODO def test_dssp_3ovv(self):
        # TODO def test_sift_3edv(self): Example dbResNum is a string, therefore was not merging.


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestTableMerger)
    unittest.TextTestRunner(verbosity=2).run(suite)
