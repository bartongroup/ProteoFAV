# -*- coding: utf-8 -*-

import os
import sys
import logging
import unittest

from proteofav.validation import (parse_validation_residues, select_validation,
                                  _fix_label_alt_id, _fix_pdb_ins_code)


class TestValidationParser(unittest.TestCase):
    """Initial testing for parsing PDB validation xml files."""

    def setUp(self):
        """Initialize the framework for testing."""
        self.example_validation = os.path.join(os.path.dirname(__file__), "testdata",
                                               "validation", "2pah_validation.xml")
        self.parser = parse_validation_residues
        self.fix_label_alt_id = _fix_label_alt_id
        self.fix_pdb_ins_code = _fix_pdb_ins_code

    def tearDown(self):
        """Remove testing framework."""
        self.example_validation = None
        self.parser = None
        self.fix_label_alt_id = None
        self.fix_pdb_ins_code = None

    def test_to_2pah(self):
        data = self.parser(self.example_validation)
        self.assertFalse(data.empty, 'Parser returned an empty DataFrame')
        self.assertEqual(data.shape, (653, 26), 'Parser returned a DataFrame '
                                                'with unexpected shape')
        self.assertEqual(data.chain.unique().tolist(), ['A', 'B'])
        self.assertEqual(data.loc[626, 'resname'], 'THR')

    def test_fix_label_alt_id(self):
        data = self.parser(self.example_validation, fix_label_alt_id=False)
        self.assertEqual('118', data.loc[0, 'resnum'])
        self.assertEqual(' ', data.loc[0, 'altcode'])
        data = self.fix_label_alt_id(data)
        self.assertEqual('118', data.loc[0, 'resnum'])
        self.assertEqual('.', data.loc[0, 'altcode'])

    def test_fix_type_symbol(self):
        data = self.parser(self.example_validation, fix_ins_code=False)
        self.assertEqual('118', data.loc[0, 'resnum'])
        self.assertEqual(' ', data.loc[0, 'icode'])
        data = self.fix_pdb_ins_code(data)
        self.assertEqual('118', data.loc[0, 'resnum'])
        self.assertEqual('?', data.loc[0, 'icode'])


if __name__ == '__main__':
    logging.basicConfig(stream=sys.stderr)
    logging.getLogger("proteofav.config").setLevel(logging.CRITICAL)
    suite = unittest.TestLoader().loadTestsFromTestCase(TestValidationParser)
    unittest.TextTestRunner(verbosity=2).run(suite)
