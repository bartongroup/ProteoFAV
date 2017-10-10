# -*- coding: utf-8 -*-

import os
import sys
import logging
import unittest

try:
    from mock import patch
except ImportError:
    # python 3.5
    from unittest.mock import patch

from proteofav.library import scop_3to1
from proteofav.config import defaults
from proteofav.dssp import (parse_dssp_residues, _import_dssp_chains_ids,
                            select_dssp, filter_dssp, get_rsa, get_rsa_class,
                            _add_dssp_rsa, _add_dssp_rsa_class,
                            _add_dssp_ss_reduced, _add_dssp_full_chain,
                            download_dssp, DSSP)

root = os.path.abspath(os.path.dirname(__file__))
defaults.db_dssp = os.path.join(root, "testdata", "dssp")


@patch("proteofav.dssp.defaults", defaults)
class TestDSSPParser(unittest.TestCase):
    """Test the DSSP parser methods."""

    def setUp(self):
        """Initialize the framework for testing."""
        self.example_dssp = os.path.join(os.path.dirname(__file__), "testdata",
                                         "dssp", "1iej.dssp")
        self.example_dssp2 = os.path.join(os.path.dirname(__file__), "testdata",
                                          "dssp", "2pah.dssp")
        self.example_dssp_bio = os.path.join(os.path.dirname(__file__), "testdata",
                                             "dssp", "2pah_bio.dssp")
        self.dssp_ins_code = os.path.join(os.path.dirname(__file__), "testdata",
                                          "dssp", "3mg7.dssp")
        self.output_dssp = os.path.join(os.path.dirname(__file__), "testdata",
                                        "2pah.dssp")
        self.residues_parser = parse_dssp_residues
        self.fix_dssp_ignoring_chains_ids = _import_dssp_chains_ids
        self.pdbid = '2pah'
        self.filter_dssp = filter_dssp
        self.add_full_chain = _add_dssp_full_chain
        self.add_rsa = _add_dssp_rsa
        self.add_rsa_class = _add_dssp_rsa_class
        self.add_ss_reduced = _add_dssp_ss_reduced
        self.get_rsa = get_rsa
        self.get_rsa_class = get_rsa_class
        self.download_dssp = download_dssp
        self.DSSP = DSSP

    def tearDown(self):
        """Remove testing framework."""

        self.example_dssp = None
        self.example_dssp2 = None
        self.example_dssp_bio = None
        self.dssp_ins_code = None
        self.output_dssp = None
        self.residues_parser = None
        self.fix_dssp_ignoring_chains_ids = None
        self.filter_dssp = None
        self.add_full_chain = None
        self.add_rsa = None
        self.add_rsa_class = None
        self.add_ss_reduced = None
        self.get_rsa = None
        self.get_rsa_class = None
        self.download_dssp = None
        self.DSSP = None

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
        self.assertEqual(len(data.columns.values), 12)
        # check whether there are particular keys
        self.assertIn('CHAIN', data.columns.values)
        # check the values for particular entries
        self.assertTrue(data.loc[0, 'CHAIN'] == 'A')
        self.assertEqual(data.loc[7, 'SS'], 'E')
        self.assertEqual(data.loc[0, 'ACC'], 179)

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
            self.assertItemsEqual(table.loc[dssp_has_seq, 'CHAIN'].unique(), chains_4v9d)
        except NameError:
            # python 3.5
            self.assertCountEqual(table.loc[dssp_has_seq, 'CHAIN'].unique(), chains_4v9d)

    def test_to_table_dssp_3mg7(self):
        """
        Getting some letters in the accession column.
        """
        self.data = self.residues_parser(self.dssp_ins_code)
        self.assertFalse(self.data.empty, 'Empty DataFrame for example with '
                                          'insertion code.')

        self.assertEqual(self.data.loc[3302, 'RES_FULL'], '102A')
        self.assertEqual(self.data.loc[3302, 'AA'], 'I')
        self.assertEqual(self.data.loc[3302, 'SS'], 'H')
        self.assertEqual(self.data.loc[3302, 'ACC'], 55)
        self.assertEqual(self.data.loc[3302, 'PHI'], -92.9)
        self.assertEqual(self.data.loc[3302, 'PSI'], -50.1)
        self.assertNotEqual(self.data.loc[3302, 'RES'], 102)
        self.assertEqual(self.data.loc[3302, 'RES'], '102')
        self.assertEqual(self.data.loc[3302, 'INSCODE'], 'A')

        self.assertEqual(self.data.loc[6401, 'RES_FULL'], '187J')
        self.assertEqual(self.data.loc[6401, 'AA'], 'L')
        self.assertEqual(self.data.loc[6401, 'SS'], '')
        self.assertEqual(self.data.loc[6401, 'ACC'], 92)
        self.assertEqual(self.data.loc[6401, 'PHI'], -57.8)
        self.assertEqual(self.data.loc[6401, 'PSI'], 360)
        self.assertNotEqual(self.data.loc[6401, 'PSI'], 800)

    def test_empty(self):
        with self.assertRaises(ValueError):
            self.residues_parser(os.path.join(os.path.dirname(__file__),
                                              "testdata", "dssp", "empty.dssp"))

    def test_parser_dssp(self):
        data = self.residues_parser(self.example_dssp2)
        self.assertEqual(data.loc[0, 'CHAIN'], 'A')
        self.assertEqual(data.loc[0, 'RES'], '118')
        self.assertEqual(data.loc[331, 'CHAIN'], 'B')
        self.assertEqual(data.loc[331, 'RES'], '118')

    def test_parser_dssp_bio(self):
        data = self.residues_parser(self.example_dssp_bio)
        self.assertEqual(data.loc[0, 'CHAIN'], 'A')
        self.assertEqual(data.loc[0, 'RES'], '118')
        self.assertEqual(data.loc[662, 'CHAIN'], 'B')
        self.assertEqual(data.loc[662, 'RES'], '118')

    def test_filter_dssp_excluded_cols(self):
        data = self.residues_parser(self.example_dssp2, excluded_cols=())
        keys = [k for k in data]
        self.assertIn("LINE", keys)
        self.assertIn("STRUCTURE", keys)
        self.assertIn("BP1", keys)
        self.assertIn("BP2", keys)
        self.assertIn("BP2_CHAIN", keys)
        self.assertIn("X-CA", keys)
        self.assertIn("Y-CA", keys)
        self.assertIn("Z-CA", keys)
        exc_cols = ("LINE", "STRUCTURE", "BP1", "BP2", "BP2_CHAIN", "X-CA", "Y-CA", "Z-CA")
        data = self.filter_dssp(data, excluded_cols=exc_cols)
        keys = [k for k in data]
        self.assertNotIn("LINE", keys)
        self.assertNotIn("STRUCTURE", keys)
        self.assertNotIn("BP1", keys)
        self.assertNotIn("BP2", keys)
        self.assertNotIn("BP2_CHAIN", keys)
        self.assertNotIn("X-CA", keys)
        self.assertNotIn("Y-CA", keys)
        self.assertNotIn("Z-CA", keys)

    def test_filter_dssp_add_chain_full(self):
        data = self.residues_parser(self.example_dssp_bio)
        data = self.add_full_chain(data)
        self.assertIn("AA", data.CHAIN_FULL.unique())
        self.assertIn("BA", data.CHAIN_FULL.unique())
        data = self.residues_parser(self.example_dssp_bio)
        data = self.filter_dssp(data, add_full_chain=True)
        self.assertIn("AA", data.CHAIN_FULL.unique())
        self.assertIn("BA", data.CHAIN_FULL.unique())

    def test_filter_dssp_add_rsa(self):
        data = self.residues_parser(self.example_dssp_bio)
        data = self.add_rsa(data)
        self.assertEqual(52.863, data.loc[2, 'RSA'])
        data = self.residues_parser(self.example_dssp_bio)
        data = self.filter_dssp(data, add_rsa=True)
        self.assertEqual(52.863, data.loc[2, 'RSA'])

    def test_filter_dssp_add_ss_reduced(self):
        data = self.residues_parser(self.example_dssp_bio)
        data = self.add_ss_reduced(data)
        self.assertEqual('H', data.loc[29, 'SS_CLASS'])
        data = self.residues_parser(self.example_dssp_bio)
        data = self.filter_dssp(data, add_ss_reduced=True)
        self.assertEqual('H', data.loc[29, 'SS_CLASS'])

    def test_filter_dssp_add_rsa_class(self):
        data = self.residues_parser(self.example_dssp_bio)
        data = self.add_rsa(data)
        data = self.add_rsa_class(data)
        self.assertEqual('Surface', data.loc[2, 'RSA_CLASS'])
        data = self.residues_parser(self.example_dssp_bio)
        data = self.filter_dssp(data, add_rsa=True, add_rsa_class=True)
        self.assertEqual('Surface', data.loc[2, 'RSA_CLASS'])

    def test_filter_dssp_chain(self):
        data = self.residues_parser(self.example_dssp2)
        data = self.filter_dssp(data, chains=('A',))
        self.assertNotIn("B", data.CHAIN.unique())

    def test_filter_dssp_chain_full(self):
        data = self.residues_parser(self.example_dssp_bio)
        data = self.filter_dssp(data, chains_full=('BA',))
        self.assertIn("B", data.CHAIN.unique())
        self.assertNotIn("B", data.CHAIN_FULL.unique())

    def test_filter_dssp_res(self):
        data = self.residues_parser(self.example_dssp2)
        data = self.filter_dssp(data, res=('118',))
        self.assertNotIn('119', data.RES.unique())

    def test_get_rsa(self):
        rsa = self.get_rsa(10.0, "A", method="Sander")
        self.assertEqual(9.434, rsa)
        rsa = self.get_rsa(20.0, "A", method="Miller")
        self.assertEqual(17.699, rsa)
        rsa = self.get_rsa(30.0, "A", method="Wilke")
        self.assertEqual(23.256, rsa)

    def test_get_rsa_class(self):
        rsa_class = self.get_rsa_class(25.5)
        self.assertEqual('Surface', rsa_class)
        rsa_class = self.get_rsa_class(7.5)
        self.assertEqual('Part. Exposed', rsa_class)
        rsa_class = self.get_rsa_class(1.5)
        self.assertEqual('Core', rsa_class)

    def test_download_dssp(self):
        self.download_dssp(self.pdbid, filename=self.output_dssp,
                           overwrite=True)
        if os.path.exists(self.output_dssp):
            os.remove(self.output_dssp)

    def test_main_DSSP(self):
        # read
        data = self.DSSP.read(self.example_dssp2)
        self.assertEqual(data.loc[0, 'CHAIN'], 'A')
        self.assertEqual(data.loc[0, 'RES'], '118')
        self.assertEqual(data.loc[331, 'CHAIN'], 'B')
        self.assertEqual(data.loc[331, 'RES'], '118')
        # download
        self.DSSP.download(self.pdbid, filename=self.output_dssp,
                           overwrite=True)
        if os.path.exists(self.output_dssp):
            os.remove(self.output_dssp)
        # select
        self.DSSP.select(self.pdbid)


if __name__ == '__main__':
    logging.basicConfig(stream=sys.stderr)
    logging.getLogger("proteofav.config").setLevel(logging.CRITICAL)
    suite = unittest.TestLoader().loadTestsFromTestCase(TestDSSPParser)
    unittest.TextTestRunner(verbosity=2).run(suite)
