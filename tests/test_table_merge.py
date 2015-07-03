#!/local/bin/python
# -*- coding: utf-8 -*-

__author__ = 'tbrittoborges'
__version__ = "1.0"

from os import path
import unittest

from to_table import _dssp_to_table, _sifts_residues_to_table, _mmcif_atom_to_table

from main import merge_tables
from config import Defaults


class TestTableMeger(unittest.TestCase):
    """Test table merging methodsthe DSSP parser methods."""

    def setUp(self):
        """Initialize the framework for testing."""
        self.defaults = Defaults("config.txt")
        self.defaults.db_mmcif = "CIF"
        self.defaults.db_dssp = "DSSP"
        self.defaults.db_sifts = "SIFTS"

        self.cif_to_table = _mmcif_atom_to_table
        self.sifts_to_table = _sifts_residues_to_table
        self.dssp_to_table = _dssp_to_table

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

        data = self.merge_table("Q16566", "2w4o", "A", defaults=self.defaults)
        self.assertIsNotNone(data)

        self.cif_path = path.join(path.dirname(__file__), "CIF/2w4o.cif")
        self.sifts_path = path.join(path.dirname(__file__), "SIFTS/2w4o.xml")
        self.dssp_path = path.join(path.dirname(__file__), "DSSP/2w4o.dssp")
        self.cif = self.cif_to_table(self.cif_path)
        self.sifts = self.sifts_to_table(self.sifts_path)
        self.dssp = self.dssp_to_table(self.dssp_path)

        # number of rows (residues) is determined by the number of residues in
        # cif file
        chain = self.cif.query(
            "label_asym_id == 'A' & group_PDB == 'ATOM' & label_atom_id == 'CA'")
        n_residues = chain.shape[0]

        self.assertEqual(data.shape[0], n_residues)

        # number of columns is given buy sum of cols in cif, dssp and sifts
        n_cols = self.sifts.shape[1] + self.dssp.shape[1] + 6  # dssp.shape[1] = 8 - 2 inds

        self.assertEqual(data.shape[1], n_cols,
                         "Incorrect number of cols in camIV table in CA mode:"
                         "{} instead {}.".format(data.shape[1], n_cols))

        self.assertTrue(any(self.cif.auth_seq_id))
        self.assertTrue(any(self.sifts.UniProt_dbResNum))
        self.assertTrue(any(self.dssp.ss))

    def test_camIV_list_mode(self):
        pass

    def test_camIV_centroid_mode(self):
        pass

    def test_dssp_3ovv(self):
        pass


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestTableMeger)
    unittest.TextTestRunner(verbosity=2).run(suite)
