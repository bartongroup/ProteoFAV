#!/local/bin/python
# -*- coding: utf-8 -*-


import unittest
from os import path

from proteofav.structures import sifts_best
from proteofav.parsers import parse_sifts_residues_from_file
from proteofav.utils import (_pdb_uniprot_sifts_mapping,
                             _uniprot_pdb_sifts_mapping)


class TestSIFTSParser(unittest.TestCase):
    """Test the SIFTS parser methods."""

    def setUp(self):
        """Initialize the framework for testing."""
        self.example_xml = path.join(path.dirname(__file__), "SIFTS/2pah.xml")
        self.residues_parser = parse_sifts_residues_from_file

        self.pdb_id = '2pah'
        self.uniprot_id = 'P00439'
        self.pdb_uniprot = _pdb_uniprot_sifts_mapping
        self.uniprot_pdb = _uniprot_pdb_sifts_mapping
        self.pdb_best = sifts_best

    def tearDown(self):
        """Remove testing framework."""

        self.example_xml = None
        self.residue_parser = None

        self.pdb_id = None
        self.uniprot_id = None
        self.pdb_uniprot = None
        self.uniprot_uniprot = None
        self.pdb_best = None

    def test_to_table_sifts_residues_and_regions(self):
        """
        Tests the parsing real SIFTS xml files.
        This test focuses on the method that parses the residue entries
        and the region entries.

        Some checks are made to whether the parsed keys and values
        are the ones we are expecting.
        """

        data = self.residues_parser(self.example_xml,
                                    sources=('CATH', 'SCOP', 'Pfam'))

        # number of values per column (or rows)
        self.assertEqual(len(data), 670)

        # number of keys (or columns)
        self.assertEqual(len(data.columns.values), 17)

        # check whether there are particular keys
        self.assertIn('CATH_dbAccessionId', data.columns.values)

        # check the values for particular entries
        self.assertTrue(data['CATH_dbAccessionId'][0] == '1.10.800.10')
        self.assertEqual(data['PDB_dbChainId'][0], 'A')
        self.assertEqual(data['PDB_dbResName'][0], 'VAL')

        # check whether there are particular keys
        self.assertIn('UniProt_regionId', data.columns.values)
        self.assertIn('CATH_regionId', data.columns.values)

        # check the values of particular entries
        self.assertTrue(data['CATH_regionId'][0] == '1')

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

        # check whether there are particular keys
        self.assertIn('pdb_id', data.columns.values)

        # check the values of particular entries
        self.assertEqual(data['pdb_id'].unique()[0], '2pah')
        self.assertIn('A', data['chain_id'].unique())
        self.assertIn('X-ray diffraction', data['experimental_method'].unique())
        self.assertTrue(type(data['coverage'][0]), float)
        self.assertTrue(type(data['resolution'][0]), float)
        self.assertTrue(type(data['tax_id'][0]), int)

    def test_to_table_uniprot_pdb_sifts_best(self):
        """
        Testing the PDBe API for mapping between UniProt and PDB
        accession identifiers.
        """

        data = self.pdb_best(self.uniprot_id, first=True)

        # number of values
        self.assertEqual(len(data), 10)

        # number of keys (or columns)
        self.assertIn('coverage', data)
        self.assertIn('pdb_id', data)
        self.assertIn('chain_id', data)


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestSIFTSParser)
    unittest.TextTestRunner(verbosity=2).run(suite)
