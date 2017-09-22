#!/local/bin/python
# -*- coding: utf-8 -*-


import unittest
from os import path

import logging
import numpy

from proteofav.library import scop_3to1
from proteofav.structures import _import_dssp_chains_ids
from proteofav.parsers import parse_dssp_from_file

logging.getLogger('proteofav').setLevel(logging.CRITICAL)  # turn off logging


class TestDSSPParser(unittest.TestCase):
    """Test the DSSP parser methods."""

    def setUp(self):
        """Initialize the framework for testing."""
        self.example_dssp = path.join(path.dirname(__file__), "DSSP/1iej.dssp")
        self.dssp_ins_code = path.join(path.dirname(__file__), "DSSP/3mg7.dssp")
        self.residues_parser = parse_dssp_from_file
        self.fix_dssp_ignoring_chains_ids = _import_dssp_chains_ids

    def tearDown(self):
        """Remove testing framework."""

        self.example_dssp = None
        self.dssp_ins_code = None
        self.residues_parser = None
        self.fix_dssp_ignoring_chains_ids = None

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

    @unittest.expectedFailure  # FIXME
    def test_fix_dssp_ignoring_chains_ids_has_as_many_chains(self):
        fix = self.fix_dssp_ignoring_chains_ids
        chains_4v9d = (
            'AB', 'AC', 'AD', 'AE', 'AF', 'AG', 'AH', 'AI', 'AJ', 'AK', 'AL', 'AM', 'AN', 'AO',
            'AP', 'AQ', 'AR', 'AS', 'AT', 'AU', 'AY', 'BB', 'BC', 'BD', 'BE', 'BF', 'BG', 'BH',
            'BI', 'BJ', 'BK', 'BL', 'BM', 'BN', 'BO', 'BP', 'BQ', 'BR', 'BS', 'BT', 'BU', 'CC',
            'CD', 'CE', 'CF', 'CG', 'CH', 'CI', 'CJ', 'CK', 'CL', 'CM', 'CN', 'CO', 'CP', 'CQ',
            'CR', 'CS', 'CT', 'CU', 'CV', 'CW', 'CX', 'CY', 'CZ', 'C0', 'C1', 'C2', 'C3', 'C4',
            'DC', 'DD', 'DE', 'DF', 'DG', 'DH', 'DI', 'DJ', 'DK', 'DL', 'DM', 'DN', 'DO', 'DP',
            'DQ', 'DR', 'DS', 'DT', 'DU', 'DV', 'DW', 'DX', 'DY', 'DZ', 'D0', 'D1', 'D2', 'D3',
            'D4')
        table = fix('4v9d')
        dssp_has_seq = table.aa.isin(scop_3to1.values())
        try:
            self.assertItemsEqual(table.loc[dssp_has_seq, 'chain_id'].unique(), chains_4v9d)
        except NameError:
            # python 3.5
            self.assertCountEqual(table.loc[dssp_has_seq, 'chain_id'].unique(), chains_4v9d)

    def test_to_table_dssp_3mg7(self):
        """
        Getting some letters in the accession column.
        """
        self.data = self.residues_parser(self.dssp_ins_code)
        self.assertFalse(self.data.empty, 'Empty DataFrame for example with '
                                          'insertion code.')

        self.assertEqual(self.data.loc[3303, 'icode'], '102A')
        self.assertEqual(self.data.loc[3303, 'aa'], 'I')
        self.assertEqual(self.data.loc[3303, 'ss'], 'H')
        self.assertEqual(self.data.loc[3303, 'acc'], 55)
        self.assertEqual(self.data.loc[3303, 'phi'], -92.9)
        self.assertEqual(self.data.loc[3303, 'psi'], -50.1)
        self.assertNotEqual(self.data.loc[3303, 'icode'], 102)

        self.assertEqual(self.data.loc[6402, 'icode'], '187J')
        self.assertEqual(self.data.loc[6402, 'aa'], 'L')
        self.assertTrue(numpy.isnan(self.data.loc[6402, 'ss']))
        self.assertEqual(self.data.loc[6402, 'acc'], 92)
        self.assertEqual(self.data.loc[6402, 'phi'], -57.8)
        self.assertEqual(self.data.loc[6402, 'psi'], 360)
        self.assertNotEqual(self.data.loc[6402, 'psi'], 800)

    def test_empty(self):
        with self.assertRaises(ValueError):
            self.residues_parser(path.join(path.dirname(__file__), "DSSP/empty.dssp"))


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestDSSPParser)
    unittest.TextTestRunner(verbosity=2).run(suite)
