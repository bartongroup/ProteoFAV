#!/local/bin/python
# -*- coding: utf-8 -*-

"""
Created on 09/06/2015

"""
from os import path
import unittest
import numpy

from to_table import _dssp_to_table

__version__ = "1.0"


class TestDSSPParser(unittest.TestCase):
    """Test the DSSP parser methods."""

    def setUp(self):
        """Initialize the framework for testing."""
        self.example_dssp = path.join(path.dirname(__file__), "DSSP/1iej.dssp")
        self.dssp_ins_code = path.join(path.dirname(__file__), "DSSP/3mg7.dssp")
        self.residues_parser = _dssp_to_table

    def tearDown(self):
        """Remove testing framework."""

        self.example_dssp = None
        self.residues_parser = None

    def test_to_table_dssp_residues(self):
        """
        Tests the parsing real DSSP files.‰

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

    # def test_to_table_dssp_3j3q(self):
    #     pass
    #
    # def test_to_table_dssp_4v60(self):
    #     pass
    #
    # def test_to_table_dssp_3kic(self):
    #     pass
    #     # test if the number of residues per chain is correct
    #     # check if all chains are there

    def test_to_table_dssp_3mg7(self):
        """
        Getting some letters in the accession column.
        """
        self.data = self.residues_parser(self.dssp_ins_code)
        self.assertFalse(self.data.empty, 'Empty DataFrame for example with '
                                          'insertion code.')

        self.assertEqual(self.data.ix[3303].icode, '102A')
        self.assertEqual(self.data.ix[3303].aa, 'I')
        self.assertEqual(self.data.ix[3303].ss, 'H')
        self.assertEqual(self.data.ix[3303].acc, 55)
        self.assertEqual(self.data.ix[3303].phi, -92.9)
        self.assertEqual(self.data.ix[3303].psi, -50.1)
        self.assertNotEqual(self.data.ix[3303].icode, 102)

        self.assertEqual(self.data.ix[6402].icode, '187J')
        self.assertEqual(self.data.ix[6402].aa, 'L')
        self.assertTrue(numpy.isnan(self.data.ix[6402].ss))
        self.assertEqual(self.data.ix[6402].acc, 92)
        self.assertEqual(self.data.ix[6402].phi, -57.8)
        self.assertEqual(self.data.ix[6402].psi, 360)
        self.assertNotEqual(self.data.ix[6402].psi, 800)

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestDSSPParser)
    unittest.TextTestRunner(verbosity=2).run(suite)
