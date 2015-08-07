from os import path
from to_table import _pdb_validation_to_table

__author__ = 'tbrittoborges'

import unittest


class TestValidationParser(unittest.TestCase):
    """Initial testing for parsing PDB validation xml files."""

    def setUp(self):
        """Initialize the framework for testing."""
        self.path_2pah = path.join(path.dirname(__file__),
                                   "/validation_xml/2pah_validation.xml")
        self.parser = _pdb_validation_to_table


    def tearDown(self):
        """Remove testing framework."""
        self.path_2pah = None
        self.parser = None

    def test_to_2pah(self):
        """
        """
        data = self.parser(self.path_2pah)
        self.assertFalse(data.empty, 'Parser returned an empty DataFrame')
        self.assertEqual(data.shape, (653, 26), 'Parser returned a DataFrame '
                                                'with unespected shape')
        self.assertEqual(data.chain.unique().tolist(), ['A', 'B'])
        self.assertEqual(data.loc[626, 'resname'], 'THR')


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestValidationParser)
    unittest.TextTestRunner(verbosity=2).run(suite)
