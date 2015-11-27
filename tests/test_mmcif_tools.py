#!/local/bin/python
# -*- coding: utf-8 -*-


import unittest
from os import path

from structures.mmcif_tools import _bio_unit_parse_operation_expression
from structures.mmcif_tools import _bio_unit_prepare_operation
from structures.mmcif_tools import _bio_unit_to_mmcif
from structures.mmcif_tools import _bio_unit_to_table
from structures.mmcif_tools import _mmcif_info_to_dict
from structures.to_table import _mmcif_atom


class TestMMCIFTools(unittest.TestCase):
    """Test the mmCIF methods to generate biological assemblies
    and calculate contacts (by the means of distance and others)."""

    def setUp(self):
        """Initialize the framework for testing."""
        self.example_mmcif = path.join(path.dirname(__file__), "CIF/2pah.cif")
        self.mmcif_atom_parser = _mmcif_atom
        self.mmcif_info_parser = _mmcif_info_to_dict
        self.bio_unit_builder = _bio_unit_to_table
        self.bio_unit_mmcif_builder = _bio_unit_to_mmcif
        self.parse_operaton = _bio_unit_parse_operation_expression
        self.prepare_operaton = _bio_unit_prepare_operation
        pass

    def tearDown(self):
        """Remove testing framework."""
        self.example_mmcif = None
        self.mmcif_atom_parser = None
        self.mmcif_info_parser = None
        self.bio_unit_builder = None
        self.bio_unit_mmcif_builder = None
        self.parse_operaton = None
        self.prepare_operaton = None
        pass

    def test_mmcif_bio_units(self):
        # TODO
        pass

    def test_parse_operation_expression(self):
        # TODO
        pass

    def test_prepare_operation(self):
        # TODO
        pass

    def test_bio_unit_to_mmcif(self):
        # TODO
        pass

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestMMCIFTools)
    unittest.TextTestRunner(verbosity=2).run(suite)
