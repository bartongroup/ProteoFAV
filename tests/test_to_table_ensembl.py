#!/local/bin/python
# -*- coding: utf-8 -*-

"""
Created on 11/06/2015

"""

__version__ = "1.0"

import unittest
from variants import to_table
from utils import utils


class TestENSEMBLParser(unittest.TestCase):
    """Test Ensembl fetcher/parser methods."""

    def setUp(self):
        """Initialize the framework for testing."""
        self.ensembl_id = 'ENSP00000326864'
        self.ensembl_id_error1 = ''
        self.ensembl_id_error2 = '1234 ads'
        self.ensembl_id_error3 = 1234567890102340
        self.ensembl_id_error4 = ()
        self.ensembl_id_error5 = []
        self.variant_id = 'rs557625940'
        self.variant_id_error1 = ''
        self.variant_id_error2 = 123456
        self.variant_id_error3 = []
        self.ensembl_trascript = to_table._transcript_variants_ensembl_to_table
        self.ensembl_somatic = to_table._somatic_variants_ensembl_to_table
        self.ensembl_variant = to_table._ensembl_variant_to_table
        self.isvalid = utils.isvalid_ensembl

    def tearDown(self):
        """Remove testing framework."""

        self.ensembl_id = None
        self.ensembl_id_error1 = None
        self.ensembl_id_error2 = None
        self.ensembl_id_error3 = None
        self.ensembl_id_error4 = None
        self.ensembl_id_error5 = None
        self.variant_id = None
        self.variant_id_error1 = None
        self.variant_id_error2 = None
        self.variant_id_error3 = None
        self.ensembl_trascript = None
        self.ensembl_somatic = None
        self.ensembl_variant = None

    def test_ensembl_ids(self):
        """
        Testing input of invalid Ensembl/Variant identifiers.
        """

        self.assertTrue(self.isvalid(self.ensembl_id))
        self.assertFalse(self.isvalid(self.ensembl_id_error1))
        self.assertFalse(self.isvalid(self.ensembl_id_error2))
        self.assertFalse(self.isvalid(self.ensembl_id_error3))
        self.assertFalse(self.isvalid(self.ensembl_id_error4))
        self.assertFalse(self.isvalid(self.ensembl_id_error5))

        self.assertTrue(self.isvalid(self.variant_id, variant=True))
        self.assertFalse(self.isvalid(self.variant_id_error1, variant=True))
        self.assertFalse(self.isvalid(self.variant_id_error2, variant=True))
        self.assertFalse(self.isvalid(self.variant_id_error3, variant=True))

    def test_to_table_transcript_variants_ensembl(self):
        """
        Tests the fetching and parsing real Ensembl Protein ids.

        Some checks are made to whether the parsed keys and values
        are the ones we are expecting.
        """

        data = self.ensembl_trascript(self.ensembl_id, verbose=False)

        # number of values per column (or rows)
        self.assertEqual(len(data), 127)

        # number of keys (or columns)
        self.assertEqual(len(data.columns.values), 14)

        # check whether there are particular keys
        self.assertIn('Parent', data.columns.values)

        # check the values for particular entries
        self.assertTrue(data['Parent'][0] == 'ENST00000321944')
        self.assertTrue(data['feature_type'][0] == "transcript_variation")
        self.assertEqual(data['codons'][0], 'Gcg/Acg')
        self.assertEqual(data['start'][0], 328)

    def test_to_table_somatic_variants_ensembl(self):
        """
        Tests the fetching and parsing real Ensembl Protein ids.

        Some checks are made to whether the parsed keys and values
        are the ones we are expecting.
        """

        data = self.ensembl_somatic(self.ensembl_id, verbose=False)

        # number of values per column (or rows)
        self.assertEqual(len(data), 59)

        # number of keys (or columns)
        self.assertEqual(len(data.columns.values), 14)

        # check whether there are particular keys
        self.assertIn('Parent', data.columns.values)

        # check the values for particular entries
        self.assertTrue(data['Parent'][0] == 'ENST00000321944')
        self.assertTrue(data['feature_type'][0] == "somatic_transcript_variation")
        self.assertEqual(data['codons'][0], 'gaC/gaT')
        self.assertEqual(data['start'][0], 466)

    def test_to_table_ensembl_variant(self):
        """
        Tests the fetching and parsing real Ensembl Variant ids.

        Some checks are made to whether the parsed keys and values
        are the ones we are expecting.
        """

        data = self.ensembl_variant(self.variant_id, verbose=False)

        # number of values per column (or rows)
        self.assertEqual(len(data), 1)

        # number of keys (or columns)
        self.assertEqual(len(data.columns.values), 17)

        # check whether there are particular keys
        self.assertIn('location', data.columns.values)

        # check the values for particular entries
        self.assertEqual(data['location'][0], '19:39175331-39175331')
        self.assertTrue(data['name'][0] == "rs557625940")
        self.assertEqual(data['seq_region_name'][0], '19')
        self.assertEqual(data['evidence'][0], ['Frequency', '1000Genomes'])

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestENSEMBLParser)
    unittest.TextTestRunner(verbosity=2).run(suite)
