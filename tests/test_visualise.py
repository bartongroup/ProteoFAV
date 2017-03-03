# coding=utf-8
import unittest

import pandas as pd

from proteofav.visualise import make_chimera_attribute_file, make_chimera_command_file


class TestVisualiser(unittest.TestCase):
    def setUp(self):
        """Initialize the framework for testing."""
        # self.example_dssp = path.join(path.dirname(__file__), "DSSP/1iej.dssp")
        self.chimera_attribute = make_chimera_attribute_file
        self.chimera_command = make_chimera_command_file

    def tearDown(self):
        """Remove testing framework."""
        self.chimera_attribute = None

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

        lines = self.chimera_attribute(example)
        self.assertEqual(lines, result)

    def test_raises_chimera_attribute_file(self):

        with self.assertRaises(TypeError):
            self.chimera_attribute(None)

        with self.assertRaises(NotImplementedError):
            lines = self.chimera_attribute(pd.Series(['0']), recipient='atoms')

    def test_make_chimera_command_file(self):
        result = "open here\ncolor green,r helix\ncolor yellow,r strand\ncolor gray,r coil\n"
        lines = self.chimera_command('here')
        self.assertEqual(lines, result)

        result = "open here\n"
        lines = self.chimera_command('here', color_secondary_structure=False)
        self.assertEqual(lines, result)

        result = "open here\ntest\ncolor green,r helix\ncolor yellow,r strand\ncolor gray,r coil\n"
        lines = self.chimera_command('here', content='test\n')
        self.assertEqual(lines, result)

    def test_visualise_chimera(self):
        pass  # TODO



if __name__ == '__main__':
    unittest.main()
