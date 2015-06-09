#!/local/bin/python

"""
Created on 09/06/2015

"""

__version__ = "1.0"

import unittest
from structs import to_table


class TestSIFTSParser(unittest.TestCase):
    """Test the SIFTS parser methods."""

    def setUp(self):
        """Initialize the framework for testing."""
        self.example_xml = "SIFTS/2pah.xml"
        self.residues_parser = to_table._sifts_to_table_residues
        self.regions_parser = to_table._sifts_to_table_regions

    def tearDown(self):
        """Remove testing framework."""

        self.example_xml = None
        self.residue_parser = None
        self.region_parser = None

    def test_to_table_sifts_residues(self):
        """
        Tests the parsing real SIFTS xml files.
        This test focuses on the method that parses the residue entries.

        Some checks are made to whether the parsed keys and values
        are the ones we are expecting.
        """

        data = self.residues_parser(self.example_xml)

        # number of values per column (or rows)
        self.assertEqual(len(data), 335)

        # number of keys (or columns)
        self.assertEqual(len(data.columns.values), 38)

        # check whether there are particular keys
        self.assertTrue('CATH_dbAccessionId' in data.columns.values)

        # check the values for particular entries
        self.assertTrue(data['CATH_dbAccessionId'][0] == '1.10.800.10')
        self.assertEqual(data['CATH_dbChainId'][0], 'A')
        self.assertEqual(data['CATH_dbResName'][0], 'VAL')

    def test_to_table_sifts_regions(self):
        """
        Tests the parsing real SIFTS xml files.
        This test focuses on the method that parses the region entries.

        Some checks are made to whether the parsed keys and values
        are the ones we are expecting.
        """

        data = self.regions_parser(self.example_xml)

        # number of values per column (or rows)
        self.assertEqual(len(data), 13)

        # number of keys (or columns)
        self.assertEqual(len(data.columns.values), 24)

        # check whether there are particular keys
        self.assertTrue('PDB_dbAccessionId' in data.columns.values)

        # check the values of particular entries
        self.assertTrue(data['PDB_dbAccessionId'][0] == '2pah')
        self.assertEqual(data['Start'][0], '1')
        self.assertEqual(data['End'][0], '335')

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestSIFTSParser)
    unittest.TextTestRunner(verbosity=2).run(suite)
