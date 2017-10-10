# -*- coding: utf-8 -*-

import os
import sys
import logging
import unittest

from proteofav.validation import (parse_validation_residues, select_validation,
                                  _fix_label_alt_id, _fix_pdb_ins_code,
                                  _add_validation_res_full,
                                  filter_validation, download_validation,
                                  Validation)


class TestValidationParser(unittest.TestCase):
    """Initial testing for parsing PDB validation xml files."""

    def setUp(self):
        """Initialize the framework for testing."""
        self.pdbid = "2pah"
        self.example_validation = os.path.join(os.path.dirname(__file__), "testdata",
                                               "validation", "2pah_validation.xml")
        self.output_validation = os.path.join(os.path.dirname(__file__), "testdata",
                                              "2pah_validation.xml")
        self.parser = parse_validation_residues
        self.fix_label_alt_id = _fix_label_alt_id
        self.fix_pdb_ins_code = _fix_pdb_ins_code
        self.filter_validation = filter_validation
        self.download_validation = download_validation
        self.Validation = Validation
        self.add_validation_res_full = _add_validation_res_full

    def tearDown(self):
        """Remove testing framework."""
        self.pdbid = None
        self.example_validation = None
        self.output_validation = None
        self.parser = None
        self.fix_label_alt_id = None
        self.fix_pdb_ins_code = None
        self.filter_validation = None
        self.download_validation = None
        self.Validation = None
        self.add_validation_res_full = None

    def test_to_2pah(self):
        data = self.parser(self.example_validation)
        self.assertFalse(data.empty, 'Parser returned an empty DataFrame')
        self.assertEqual(data.shape, (653, 26), 'Parser returned a DataFrame '
                                                'with unexpected shape')
        self.assertEqual(data.validation_chain.unique().tolist(), ['A', 'B'])
        self.assertEqual(data.loc[626, 'validation_resname'], 'THR')

    def test_fix_label_alt_id(self):
        data = self.parser(self.example_validation, fix_label_alt_id=False)
        self.assertEqual('118', data.loc[0, 'validation_resnum'])
        self.assertEqual(' ', data.loc[0, 'validation_altcode'])
        data = self.fix_label_alt_id(data)
        self.assertEqual('118', data.loc[0, 'validation_resnum'])
        self.assertEqual('.', data.loc[0, 'validation_altcode'])

    def test_fix_type_symbol(self):
        data = self.parser(self.example_validation, fix_ins_code=False)
        self.assertEqual('118', data.loc[0, 'validation_resnum'])
        self.assertEqual(' ', data.loc[0, 'validation_icode'])
        data = self.fix_pdb_ins_code(data)
        self.assertEqual('118', data.loc[0, 'validation_resnum'])
        self.assertEqual('?', data.loc[0, 'validation_icode'])

    def test_filter_validation_chain(self):
        data = self.parser(self.example_validation)
        data = self.filter_validation(data, chains=('A',))
        self.assertNotIn("B", data.validation_chain.unique())

    def test_filter_validation_res(self):
        data = self.parser(self.example_validation)
        data = self.filter_validation(data, res=('118',))
        self.assertNotIn('119', data.validation_resnum.unique())

    def test_filter_validation_add_res_full(self):
        data = self.parser(self.example_validation)
        self.assertIn('validation_resnum', data)
        self.assertIn('validation_icode', data)
        self.assertNotIn('validation_resnum_full', data)
        data = self.add_validation_res_full(data)
        self.assertIn('validation_resnum_full', data)
        self.assertEqual(data.loc[0, 'validation_resnum_full'], '118')
        data = self.parser(self.example_validation)
        data = self.filter_validation(data, add_res_full=True)
        self.assertEqual(data.loc[0, 'validation_resnum_full'], '118')

    def test_download_validation(self):
        self.download_validation(self.pdbid, filename=self.output_validation,
                                 overwrite=True)
        if os.path.exists(self.output_validation):
            os.remove(self.output_validation)

    def test_main_Validation(self):
        # read
        data = self.Validation.read(self.example_validation)
        self.assertEqual(data.shape, (653, 26), 'Parser returned a DataFrame '
                                                'with unexpected shape')
        self.assertEqual(data.validation_chain.unique().tolist(), ['A', 'B'])
        self.assertEqual(data.loc[626, 'validation_resname'], 'THR')
        # download
        self.Validation.download(self.pdbid, filename=self.output_validation,
                                 overwrite=True)
        if os.path.exists(self.output_validation):
            os.remove(self.output_validation)
        # select
        self.Validation.select(self.pdbid)


if __name__ == '__main__':
    logging.basicConfig(stream=sys.stderr)
    logging.getLogger("proteofav.config").setLevel(logging.CRITICAL)
    suite = unittest.TestLoader().loadTestsFromTestCase(TestValidationParser)
    unittest.TextTestRunner(verbosity=2).run(suite)
