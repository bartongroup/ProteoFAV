#!/local/bin/python
# -*- coding: utf-8 -*-

"""
Created on 10/06/2015

"""

__version__ = "1.0"

import unittest
from main import to_table


class TestMMCIFParser(unittest.TestCase):
    """Test the mmCIF parser methods."""

    def setUp(self):
        """Initialize the framework for testing."""
        self.example_mmcif = "CIF/2pah.cif"
        self.mmcif_parser = to_table._mmcif_atom_to_table

    def tearDown(self):
        """Remove testing framework."""

        self.example_mmcif = None
        self.mmcif_parser = None

    def test_to_table_mmcif(self):
        """
        Tests the parsing real mmCIF files.

        Some checks are made to whether the parsed keys and values
        are the ones we are expecting.
        """

        data = self.mmcif_parser(self.example_mmcif)

        # number of values per column (or rows)
        self.assertEqual(len(data), 5317)

        # number of keys (or columns)
        self.assertEqual(len(data.columns.values), 26)

        # check whether there are particular keys
        self.assertIn('label_asym_id', data.columns.values)

        # check the values for particular entries
        self.assertTrue(data.loc[1, 'label_asym_id'] == 'A')
        self.assertEqual(data.loc[1, 'pdbx_PDB_model_num'], 1)
        self.assertEqual(data['group_PDB'][0], 'ATOM')


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestMMCIFParser)
    unittest.TextTestRunner(verbosity=2).run(suite)
