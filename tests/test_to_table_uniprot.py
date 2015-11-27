#!/local/bin/python
# -*- coding: utf-8 -*-


import unittest

import utils
from variants.to_table import _uniprot_ensembl_mapping, _uniprot_info


class TestUNIPROTParser(unittest.TestCase):
    """Test UniProt fetcher/parser methods."""

    def setUp(self):
        """Initialize the framework for testing."""
        self.uniprot_id = 'O96013'
        self.uniprot_id_error1 = ''
        self.uniprot_id_error2 = '1234 ads'
        self.uniprot_id_error3 = 1234
        self.uniprot_id_error4 = ()
        self.uniprot_id_error5 = []
        self.uniprot_info = _uniprot_info
        self.uniprot_ensembl = _uniprot_ensembl_mapping
        self.isvalid = utils.is_valid

    def tearDown(self):
        """Remove testing framework."""

        self.uniprot_id = None
        self.uniprot_id_error1 = None
        self.uniprot_id_error2 = None
        self.uniprot_id_error3 = None
        self.uniprot_id_error4 = None
        self.uniprot_id_error5 = None
        self.uniprot_info = None
        self.uniprot_ensembl = None

    def test_uniprot_ids(self):
        """
        Testing input of invalid UniProt identifiers.
        """
        self.assertTrue(self.isvalid(self.uniprot_id, 'uniprot'))
        self.assertFalse(self.isvalid(self.uniprot_id_error1, 'uniprot'))
        self.assertFalse(self.isvalid(self.uniprot_id_error2, 'uniprot'))
        self.assertFalse(self.isvalid(self.uniprot_id_error3, 'uniprot'))
        self.assertFalse(self.isvalid(self.uniprot_id_error4, 'uniprot'))
        self.assertFalse(self.isvalid(self.uniprot_id_error5, 'uniprot'))

    def test_to_table_uniprot_info(self):
        """
        Tests the fetching and parsing real UniProt ids.
        This test focuses on the method that parses the residue entries.

        Some checks are made to whether the parsed keys and values
        are the ones we are expecting.
        """

        data = self.uniprot_info(self.uniprot_id)

        # number of values per column (or rows)
        self.assertEqual(len(data), 1)

        # number of keys (or columns)
        self.assertEqual(len(data.columns.values), 7)

        # check whether there are particular keys
        self.assertIn('Sequence', data.columns.values)

        # check the values for particular entries
        self.assertTrue(data['Length'][0] == 591)
        self.assertTrue(data['Status'][0] == "reviewed")
        self.assertEqual(data['Entry name'][0], 'PAK4_HUMAN')

    def test_to_table_uniprot_ensembl_mapping(self):
        """
        Tests the fetching and parsing real UniProt ids.
        This test focuses on the method that parses the residue entries.

        Some checks are made to whether the parsed keys and values
        are the ones we are expecting.
        """

        data = self.uniprot_ensembl(self.uniprot_id)

        # number of values per column (or rows)
        self.assertEqual(len(data), 1)

        # number of keys (or columns)
        self.assertEqual(len(data.columns.values), 3)

        # check whether there are particular keys
        self.assertIn('TRANSCRIPT', data.columns.values)

        # check the values for particular entries
        self.assertEqual(data['GENE'][0], 'ENSG00000130669')


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestUNIPROTParser)
    unittest.TextTestRunner(verbosity=2).run(suite)
