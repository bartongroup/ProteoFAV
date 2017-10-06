# -*- coding: utf-8 -*-

import os
import unittest
import pandas as pd

from proteofav.variants import _fetch_ensembl_variants


class TestENSEMBLParser(unittest.TestCase):
    """Test Ensembl fetcher/parser methods."""

    def setUp(self):
        """Initialize the framework for testing."""
        self.uniprot_id = 'O96013'
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
        self.ensembl_trascript = _fetch_ensembl_variants
        self.fetch_variant = _fetch_ensembl_variants

    def tearDown(self):
        """Remove testing framework."""

        self.uniprot_id = None
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
        self.fetch_variant = None
        self.uniprot_variants = None

    def test_querying_ensembl_transcript_variants(self):
        data = self.fetch_variant(self.ensembl_id, feature='transcript_variation')

        # doesn't check anything about the output
        self.assertIsInstance(data, pd.DataFrame)

    def test_querying_ensembl_somatic_variants(self):
        data = self.fetch_variant(
            self.ensembl_id, feature='somatic_transcript_variation')

        # doesn't check anything about the output
        self.assertIsInstance(data, pd.DataFrame)

    def test_to_table_transcript_variants_ensembl(self):
        """
        Tests the fetching and parsing real Ensembl Protein ids.

        Some checks are made to whether the parsed keys and values
        are the ones we are expecting.

        Using mocking because variant data varies over time.
        For simplicity here loading the resulting pandas table from a
        dumped csv file.
        """

        # data = self.ensembl_trascript(self.ensembl_id, species='human')
        # mocking the resulting pandas dataframe
        variants = os.path.join(os.path.dirname(__file__), "testdata",
                                "variation", "transcript_variation_{}.csv".format(self.ensembl_id))
        data = pd.read_csv(variants)

        # number of values per column (or rows)
        self.assertEqual(len(data), 206)

        # number of keys (or columns)
        self.assertEqual(len(data.columns.values), 16)

        # check whether there are particular keys
        self.assertIn('Parent', data.columns.values)

        # check the values for particular entries
        self.assertTrue(data['Parent'][0] == 'ENST00000321944')
        self.assertTrue(data['feature_type'][0] == "transcript_variation")
        self.assertEqual(data['codons'][0], 'Cat/Tat')
        self.assertEqual(data['start'][0], 135)

    def test_to_table_somatic_variants_ensembl(self):
        """
        Tests the fetching and parsing real Ensembl Protein ids.

        Some checks are made to whether the parsed keys and values
        are the ones we are expecting.

        Using mocking because variant data varies over time.
        For simplicity here loading the resulting pandas table from a
        dumped csv file.
        """

        # data = self.ensembl_somatic(self.ensembl_id, species='human')
        # mocking the resulting pandas dataframe
        variants = os.path.join(os.path.dirname(__file__), "testdata", "variation",
                                "somatic_variation_{}.csv".format(self.ensembl_id))
        data = pd.read_csv(variants)

        # number of values per column (or rows)
        self.assertEqual(len(data), 40)

        # number of keys (or columns)
        self.assertEqual(len(data.columns.values), 16)

        # check whether there are particular keys
        self.assertIn('Parent', data.columns.values)

        # check the values for particular entries
        self.assertTrue(data['Parent'][0] == 'ENST00000321944')
        self.assertTrue(data['feature_type'][0] == "somatic_transcript_variation")
        self.assertEqual(data['codons'][0], 'Gag/Cag')
        self.assertEqual(data['start'][0], 119)

    def test_to_table_ensembl_variant(self):
        """
        Tests the fetching and parsing real Ensembl Variant ids.

        Some checks are made to whether the parsed keys and values
        are the ones we are expecting.

        Using mocking because variant data varies over time.
        For simplicity here loading the resulting pandas table from a
        dumped csv file.
        """

        # data = self.ensembl_variant(self.variant_id, species='human')
        # mocking the resulting pandas dataframe
        variants = os.path.join(os.path.dirname(__file__), "testdata", "variation",
                                "ensembl_variation_{}.csv".format(self.variant_id))
        data = pd.read_csv(variants)

        # number of values per column (or rows)
        self.assertEqual(len(data), 1)

        # number of keys (or columns)
        self.assertEqual(len(data.columns.values), 19)

        # check whether there are particular keys
        self.assertIn('location', data.columns.values)

        # check the values for particular entries
        self.assertEqual(data['location'][0], '19:39175331-39175331')
        self.assertTrue(data['name'][0] == "rs557625940")
        self.assertEqual(data['seq_region_name'][0], 19)
        self.assertEqual(data['most_severe_consequence'][0], 'missense_variant')

    def test_to_table_uniprot_ensembl_variants(self):
        """
        Tests the wrapper method that goes from a uniprot to a list
        of ensembl protein transcripts and to variants.

        Using mocking because variant data varies over time.
        For simplicity here loading the resulting pandas table from a
        dumped csv file.
        """

        # querying the various ensembl endpoints
        # data = self.uniprot_variants(self.uniprot_id)
        # mocking the resulting pandas dataframe
        variants = os.path.join(os.path.dirname(__file__), "testdata", "variation",
                                "uniprot_variants_{}.csv".format(self.uniprot_id))
        data = pd.read_csv(variants)

        # number of values per column (or rows)
        self.assertEqual(len(data), 903)

        # number of keys (or columns)
        self.assertEqual(len(data.columns.values), 5)

        # check whether there are particular keys
        self.assertIn('translation', data.columns.values)
        self.assertIn('id', data.columns.values)
        self.assertIn('start', data.columns.values)
        self.assertIn('residues', data.columns.values)

        # check the values for particular entries
        self.assertEqual(data['translation'][0], 'ENSP00000351049')
        self.assertTrue(data['id'][0] == "rs769098772")
        self.assertTrue(data['start'][0] == 295)
        self.assertTrue(data['residues'][0] == "E/A")


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestENSEMBLParser)
    unittest.TextTestRunner(verbosity=2).run(suite)
