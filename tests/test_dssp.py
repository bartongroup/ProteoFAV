# -*- coding: utf-8 -*-


import os
import sys
import logging
import unittest

from proteofav.dssp import (DSSP,
                            _add_dssp_full_chain,
                            _add_dssp_rsa,
                            _add_dssp_rsa_class,
                            _add_dssp_ss_reduced,
                            get_rsa,
                            get_rsa_class)

from proteofav.config import defaults as config


class TestDSSP(unittest.TestCase):
    """Test the DSSP parser methods."""

    def setUp(self):
        """Initialize the framework for testing."""

        root = os.path.abspath(os.path.dirname(__file__))

        self.pdbid = '2pah'
        self.inputdssp = os.path.join(root, "testdata", config.db_dssp,
                                      "{}.dssp".format(self.pdbid))
        self.inputbiodssp = os.path.join(root, "testdata", config.db_dssp,
                                         "{}_bio.dssp".format(self.pdbid))
        self.outputdssp = os.path.join(root, "testdata", config.db_tmp,
                                       "{}.dssp".format(self.pdbid))
        self.excluded = ()
        self.dssp = DSSP
        self.add_full_chain = _add_dssp_full_chain
        self.add_rsa = _add_dssp_rsa
        self.add_rsa_class = _add_dssp_rsa_class
        self.add_ss_reduced = _add_dssp_ss_reduced
        self.get_rsa = get_rsa
        self.get_rsa_class = get_rsa_class

        logging.disable(logging.DEBUG)

    def tearDown(self):
        """Remove testing framework."""

        self.pdbid = None
        self.inputdssp = None
        self.inputbiodssp = None
        self.outputdssp = None
        self.excluded = None
        self.dssp = None
        self.add_full_chain = None
        self.add_rsa = None
        self.add_rsa_class = None
        self.add_ss_reduced = None
        self.get_rsa = None
        self.get_rsa_class = None

        logging.disable(logging.NOTSET)

    def test_reader_cif_data(self):
        data = self.dssp().read(self.inputdssp)
        self.assertEqual(data.loc[0, 'CHAIN'], 'A')
        self.assertEqual(data.loc[0, 'RES'], '118')
        self.assertEqual(data.loc[331, 'CHAIN'], 'B')
        self.assertEqual(data.loc[331, 'RES'], '118')

    def test_reader_biocif_data(self):
        data = self.dssp().read(self.inputbiodssp)
        self.assertEqual(data.loc[0, 'CHAIN'], 'A')
        self.assertEqual(data.loc[0, 'RES'], '118')
        self.assertEqual(data.loc[662, 'CHAIN'], 'B')
        self.assertEqual(data.loc[662, 'RES'], '118')

    def test_reader_default_excluded(self):
        data = self.dssp().read(self.inputdssp)
        keys = list(data)
        self.assertNotIn("LINE", keys)
        self.assertNotIn("STRUCTURE", keys)
        self.assertNotIn("BP1", keys)
        self.assertNotIn("BP2", keys)
        self.assertNotIn("BP2_CHAIN", keys)
        self.assertNotIn("X-CA", keys)
        self.assertNotIn("Y-CA", keys)
        self.assertNotIn("Z-CA", keys)

    def test_reader_new_excluded(self):
        data = self.dssp().read(self.inputdssp,
                                excluded_cols=self.excluded)
        keys = list(data)
        self.assertIn("LINE", keys)
        self.assertIn("STRUCTURE", keys)
        self.assertIn("BP1", keys)
        self.assertIn("BP2", keys)
        self.assertIn("BP2_CHAIN", keys)
        self.assertIn("X-CA", keys)
        self.assertIn("Y-CA", keys)
        self.assertIn("Z-CA", keys)

    def test_reader_add_chain_full(self):
        data = self.dssp().read(self.inputbiodssp, add_full_chain=True)
        self.assertIn("AA", data.CHAIN_FULL.unique())
        self.assertIn("BA", data.CHAIN_FULL.unique())

    def test_reader_add_rsa(self):
        data = self.dssp().read(self.inputbiodssp, add_rsa=True)
        self.assertEqual(52.863, data.loc[2, 'RSA'])

    def test_reader_add_ss_reduced(self):
        data = self.dssp().read(self.inputbiodssp, add_ss_reduced=True)
        self.assertEqual('H', data.loc[29, 'SS_CLASS'])

    def test_reader_add_rsa_class(self):
        data = self.dssp().read(self.inputbiodssp, add_rsa=True, add_rsa_class=True)
        self.assertEqual('Surface', data.loc[2, 'RSA_CLASS'])

    def test_download_dssp(self):
        self.dssp().download(self.pdbid, filename=self.outputdssp,
                             overwrite=True, decompress=False)
        if os.path.exists(self.outputdssp):
            os.remove(self.outputdssp)

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


if __name__ == '__main__':
    logging.basicConfig(stream=sys.stderr)
    logging.getLogger("proteofav").setLevel(logging.DEBUG)
    suite = unittest.TestLoader().loadTestsFromTestCase(TestDSSP)
    unittest.TextTestRunner(verbosity=2).run(suite)
