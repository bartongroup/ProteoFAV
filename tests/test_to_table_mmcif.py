#!/local/bin/python
# -*- coding: utf-8 -*-

"""
Created on 10/06/2015

"""

__version__ = "1.0"

from os import path
import unittest

from to_table import _mmcif_atom_to_table
from mmcif_tools import _mmcif_info_to_dict

import utils


class TestMMCIFParser(unittest.TestCase):
    """Test the mmCIF parser methods."""

    def setUp(self):
        """Initialize the framework for testing."""
        self.example_mmcif = path.join(path.dirname(__file__), "CIF/2pah.cif")
        self.mmcif_atom_parser = _mmcif_atom_to_table
        self.mmcif_info_parser = _mmcif_info_to_dict

        self.pdb_id = '2pah'
        self.pdb_id_error1 = ''
        self.pdb_id_error2 = '1234 ads'
        self.pdb_id_error3 = 1234
        self.pdb_id_error4 = ()
        self.pdb_id_error5 = []
        self.isvalid = utils.isvalid_pdb

    def tearDown(self):
        """Remove testing framework."""

        self.example_mmcif = None
        self.mmcif_atom_parser = None
        self.mmcif_info_parser = None

        self.pdb_id = None
        self.pdb_id_error1 = None
        self.pdb_id_error2 = None
        self.pdb_id_error3 = None
        self.pdb_id_error4 = None
        self.pdb_id_error5 = None
        self.isvalid = None

    def test_pdb_ids(self):
        """
        Testing input of invalid UniProt identifiers.
        """
        #
        self.assertTrue(self.isvalid(self.pdb_id))
        self.assertFalse(self.isvalid(self.pdb_id_error1))
        self.assertFalse(self.isvalid(self.pdb_id_error2))
        self.assertFalse(self.isvalid(self.pdb_id_error3))
        self.assertFalse(self.isvalid(self.pdb_id_error4))
        self.assertFalse(self.isvalid(self.pdb_id_error5))

    def test_atom_to_table_mmcif(self):
        """
        Tests the parsing real mmCIF files.

        Some checks are made to whether the parsed keys and values
        are the ones we are expecting.

        Parsing of ATOM and HETATOM lines.
        """

        data = self.mmcif_atom_parser(self.example_mmcif)

        # number of values per column (or rows)
        self.assertEqual(len(data), 5317)

        # number of keys (or columns)
        self.assertEqual(len(data.columns.values), 26)

        # check whether there are particular keys
        self.assertIn('label_asym_id', data.columns.values)

        # check the values for particular entries
        self.assertTrue(data.loc[1, 'label_asym_id'] == 'A')
        self.assertEqual(data.loc[1, 'pdbx_PDB_model_num'], 1)
        self.assertEqual(data['group_PDB'][0], 'ATOM')

    def test_info_to_table_mmcif(self):
        """
        Tests the parsing real mmCIF files.

        Some checks are made to whether the parsed keys and values
        are the ones we are expecting.

        Parsing all lines other than ATOM/HETATOM lines.
        """

        data = self.mmcif_info_parser(self.example_mmcif)

        # number of main dict keys
        self.assertEqual(len(data), 66)

        # check whether there are particular main keys
        self.assertIn('pdbx_struct_assembly', data)

        # number of leading keys for a testing main key
        self.assertEqual(len(data['exptl_crystal']), 5)

        # check whether there are particular leading keys
        self.assertIn('oligomeric_details', data['pdbx_struct_assembly'])

        # check the values for particular entries
        self.assertEqual(data['pdbx_struct_assembly']['oligomeric_details'],
                         'tetrameric')
        self.assertEqual(data['pdbx_struct_assembly_gen']['asym_id_list'],
                         'A,C,B,D')
        self.assertEqual(data['pdbx_struct_oper_list']['symmetry_operation'],
                         ['x,y,z', '-y,-x,-z+2/3'])

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestMMCIFParser)
    unittest.TextTestRunner(verbosity=2).run(suite)
