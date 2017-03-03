# coding=utf-8
import unittest

import pandas as pd

from proteofav.visualise import make_chimera_attribute_file


class TestVisualiser(unittest.TestCase):
    def setUp(self):
        """Initialize the framework for testing."""
        # self.example_dssp = path.join(path.dirname(__file__), "DSSP/1iej.dssp")
        self.chimera_file = make_chimera_attribute_file

    def tearDown(self):
        """Remove testing framework."""
        self.chimera_file = None

    def test_make_chimera_attribute_file(self):

        example = pd.Series([32.4, 46.2, 52.5, 4.4, 22.1], index = range(1, 6), name='test')
        result = """# Generated with ProteoFAV
    attribute: test
    match mode: 1 - to - 1
    recipient: residues
    :1	32.4
    :2	46.2
    :3	52.5
    :4	4.4
    :5	22.1"""

        lines = self.chimera_file(example)
        self.assertEqual(lines, result)

    def test_raises_chimera_attribute_file(self):

        with self.assertRaises(TypeError):
            self.chimera_file(None)

        with self.assertRaises(NotImplementedError):
            lines = self.chimera_file(pd.Series(['0']), recipient='atoms')



if __name__ == '__main__':
    unittest.main()
