#!/local/bin/python
# -*- coding: utf-8 -*-

__author__ = 'tbrittoborges'
__version__ = "1.0"

import unittest

import to_table
from utils import merge_tables
from config import Defaults

class TestTableMeger(unittest.TestCase):
    """Test table merging methodsthe DSSP parser methods."""

    def setUp(self):
        """Initialize the framework for testing."""
        self.defaults = Defaults("config.txt")
        self.defaults.db_mmcif = "CIF"
        self.defaults.db_dssp = "DSSP"
        self.defaults.db_sifts = "SIFTS"

        self.cif_to_table = to_table._mmcif_atom_to_table
        self.sifts_to_table = to_table._sifts_residues_to_table
        self.dssp_to_table = to_table._dssp_to_table

        self.merge_table = merge_tables


    def tearDown(self):
        """Remove testing framework by cleaning the namespace."""
        self.defaults = None

        self.cif_to_table = None
        self.sifts_to_table = None
        self.dssp_to_table = None

        self.merge_table = None

    def test_camIV_ca_mode(self):
        """
        Test table merger for a protein example.
        :return:
        """
        data = self.merge_table("Q16566", defaults=self.defaults)
        self.assertIsNotNone(data)
        cif = self.cif_to_table("CIF/2w4o.cif")
        sifts = self.sifts_to_table("SIFTS/2w4o.xml")
        dssp = self.dssp_to_table("DSSP/2w4o.dssp")

        # number of rows (residues) is determiend by the number of residues in
        # cif file
        chain = cif.query("label_asym_id == 'A' & group_PDB == 'ATOM' & label_atom_id == 'CA'")
        n_residues = chain.shape[0]

        self.assertEqual(data.shape[0], n_residues)

        # number of columns is given buy sum of cols in cif, dssp and sifts
        n_cols = cif.shape[1] + sifts.shape[1] + dssp.shape[1]
        self.assertEqual(data.shape[1], n_cols, "Incorrect number of cols in camIV table in CA mode.")


    def test_camIV_list_mode(self):
        pass

    def test_camIV_centroid_mode(self):
        pass

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestTableMeger)
    unittest.TextTestRunner(verbosity=2).run(suite)