#!/local/bin/python

"""
Created on 09/06/2015

"""

__version__ = "1.0"

import unittest
from structs import to_table


class TestDSSPParser(unittest.TestCase):
    """Test the DSSP parser methods."""

    def setUp(self):
        """Initialize the framework for testing."""
        self.example_dssp = "DSSP/1iej.dssp"
        self.residues_parser = to_table._dssp_to_table

    def tearDown(self):
        """Remove testing framework."""

        self.example_dssp = None
        self.residues_parser = None

    def test_to_table_dssp_residues(self):
        """
        Tests the parsing real DSSP files.
        This test focuses on the method that parses the residue entries.

        Some checks are made to whether the parsed keys and values
        are the ones we are expecting.
        """

        data = self.residues_parser(self.example_dssp)

        # number of values per column (or rows)
        self.assertEqual(len(data), 329)

        # number of keys (or columns)
        self.assertEqual(len(data.columns.values), 7)

        # check whether there are particular keys
        self.assertIn('chain_id', data.columns.values)

        # check the values for particular entries
        self.assertTrue(data.loc[1, 'chain_id'] == 'A')
        self.assertEqual(data.loc[8, 'ss'], 'E')
        self.assertEqual(data.loc[1, 'acc'], 179)


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestDSSPParser)
    unittest.TextTestRunner(verbosity=2).run(suite)
