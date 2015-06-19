#!/local/bin/python
# -*- coding: utf-8 -*-

"""
Created on 09/06/2015

"""

__version__ = "1.0"

import unittest
from os import path
from sifts import to_table


class TestSIFTSParser(unittest.TestCase):
    """Test the SIFTS parser methods."""

    def setUp(self):
        """Initialize the framework for testing."""
        self.example_xml = path.join(path.dirname(__file__), "SIFTS/2pah.xml")
        self.residues_parser = to_table._sifts_residues_to_table
        self.regions_parser = to_table._sifts_regions_to_table

        self.pdb_id = '2pah'
        self.uniprot_id = 'P00439'
        self.pdb_uniprot = to_table._pdb_uniprot_sifts_mapping_to_table
        self.uniprot_pdb = to_table._uniprot_pdb_sifts_mapping_to_table

    def tearDown(self):
        """Remove testing framework."""

        self.example_xml = None
        self.residue_parser = None
        self.region_parser = None

        self.pdb_id = None
        self.uniprot_id = None
        self.pdb_uniprot = None
        self.uniprot_uniprot = None

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
        self.assertIn('CATH_dbAccessionId', data.columns.values)

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
        self.assertIn('PDB_dbAccessionId', data.columns.values)

        # check the values of particular entries
        self.assertTrue(data['PDB_dbAccessionId'][0] == '2pah')
        self.assertEqual(data['Start'][0], '1')
        self.assertEqual(data['End'][0], '335')

    def test_to_table_pdb_uniprot_sifts_mapping(self):
        """
        Testing the PDBe API for mapping between PDB and UniProt
        accession identifiers.
        """

        data = self.pdb_uniprot(self.pdb_id)

        # number of values per column (or rows)
        self.assertEqual(len(data), 1)

        # number of keys (or columns)
        self.assertEqual(len(data.columns.values), 1)

        # check whether there are particular keys
        self.assertIn('uniprot_id', data.columns.values)

        # check the values of particular entries
        self.assertTrue(data['uniprot_id'][0] == 'P00439')

    def test_to_table_uniprot_pdb_sifts_mapping(self):
        """
        Testing the PDBe API for mapping between UniProt and PDB
        accession identifiers.
        """

        data = self.uniprot_pdb(self.uniprot_id)

        # number of values per column (or rows)
        self.assertEqual(len(data), 17)

        # number of keys (or columns)
        self.assertEqual(len(data.columns.values), 10)

        # check whether there are particular keys
        self.assertIn('pdb_id', data.columns.values)

        # check the values of particular entries
        self.assertTrue(data['pdb_id'][0] == '2pah')
        self.assertEqual(data['chain_id'][0], 'A')
        self.assertEqual(data['experimental_method'][0], 'X-ray diffraction')
        self.assertTrue(type(data['coverage'][0]), float)
        self.assertTrue(type(data['resolution'][0]), float)
        self.assertTrue(type(data['tax_id'][0]), int)

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestSIFTSParser)
    unittest.TextTestRunner(verbosity=2).run(suite)
