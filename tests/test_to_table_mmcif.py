#!/local/bin/python
# -*- coding: utf-8 -*-


import logging
import sys
import unittest
from os import path, remove

from proteofav.mmcif_tools import _bio_unit_to_table
from proteofav.structures import _mmcif_atom
from proteofav.mmcif_tools import _mmcif_info_to_dict
from proteofav.utils import is_valid


class TestMMCIFParser(unittest.TestCase):
    """Test the mmCIF parser methods."""

    def setUp(self):
        """Initialize the framework for testing."""
        self.example_mmcif = path.join(path.dirname(__file__), "CIF/2pah.cif")
        self.mmcif_atom_parser = _mmcif_atom
        self.mmcif_info_parser = _mmcif_info_to_dict
        self.bio_unit_builder = _bio_unit_to_table
        self.example_tsv_out = path.join(path.dirname(__file__), "CIF/2pah-bio.tsv")

        self.pdb_id = '2pah'
        self.pdb_id_error1 = ''
        self.pdb_id_error2 = '1234 ads'
        self.pdb_id_error3 = 1234
        self.pdb_id_error4 = ()
        self.pdb_id_error5 = []
        self.isvalid = is_valid

    def tearDown(self):
        """Remove testing framework."""

        self.example_mmcif = None
        self.mmcif_atom_parser = None
        self.mmcif_info_parser = None
        self.bio_unit_builder = None
        self.example_tsv_out = None

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
        self.assertTrue(self.isvalid(self.pdb_id, database='pdbe'))
        self.assertFalse(self.isvalid(self.pdb_id_error1, database='pdbe'))
        self.assertFalse(self.isvalid(self.pdb_id_error2, database='pdbe'))
        self.assertFalse(self.isvalid(self.pdb_id_error3, database='pdbe'))
        self.assertFalse(self.isvalid(self.pdb_id_error4, database='pdbe'))
        self.assertFalse(self.isvalid(self.pdb_id_error5, database='pdbe'))

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

    def test_bio_unit_to_table_mmcif(self):
        """
        Tests the parsing real mmCIF files.

        Some checks are made to whether the parsed keys and values
        are the ones we are expecting.

        Parses the atom lines and the info lines and looks for
        biological assemblies. Finds the most likely and generates
        a new pandas Dataframe with the new atoms.
        """

        log = logging.getLogger("Biological.Assemblies")

        # uses _mmcif_atom, _mmcif_info_to_table, and
        # other mmcif_tools that needed to be tested

        # method == 1
        log.info("Method 1")
        data = self.bio_unit_builder(self.example_mmcif,
                                     most_likely=True,
                                     method=1)

        # number of values per column (or rows);
        # double number of rows: dimer to tetramer
        self.assertEqual(len(data), 5317 * 2)

        # number of keys (or columns); same number of columns
        self.assertEqual(len(data.columns.values), 26)

        # check whether there are particular keys
        self.assertIn('label_asym_id', data.columns.values)

        # check the values for particular entries
        self.assertTrue(data.loc[1, 'label_asym_id'] == 'A.1')
        self.assertEqual(data.loc[1, 'pdbx_PDB_model_num'], 1)
        self.assertEqual(data['group_PDB'][0], 'ATOM')

        # method == 2
        log.info("Method 2")
        data = self.bio_unit_builder(self.example_mmcif,
                                     most_likely=True,
                                     method=2)

        # number of values per column (or rows);
        # double number of rows: dimer to tetramer
        self.assertEqual(len(data), 5317 * 2)

        # number of keys (or columns); same number of columns
        self.assertEqual(len(data.columns.values), 26)

        # check whether there are particular keys
        self.assertIn('label_asym_id', data.columns.values)

        # check the values for particular entries
        self.assertTrue(data.loc[1, 'label_asym_id'] == 'A')
        self.assertEqual(data.loc[1, 'pdbx_PDB_model_num'], 1)
        self.assertEqual(data['group_PDB'][0], 'ATOM')

        # method == 3
        log.info("Method 3")
        data = self.bio_unit_builder(self.example_mmcif,
                                     most_likely=True,
                                     method=3)

        # number of values per column (or rows);
        # double number of rows: dimer to tetramer
        self.assertEqual(len(data), 5317 * 2)

        # number of keys (or columns); same number of columns
        self.assertEqual(len(data.columns.values), 27)

        # check whether there are particular keys
        self.assertIn('bio_unit_counter', data.columns.values)

        # check the values for particular entries
        self.assertTrue(data.loc[1, 'label_asym_id'] == 'A')
        self.assertEqual(data.loc[1, 'bio_unit_counter'], 1)
        self.assertEqual(data['group_PDB'][0], 'ATOM')

        # create a tsv file
        data.to_csv(self.example_tsv_out, sep='\t')
        self.assertTrue(path.isfile(self.example_tsv_out))
        remove(self.example_tsv_out)


if __name__ == '__main__':
    logging.basicConfig(stream=sys.stderr)
    logging.getLogger("Biological.Assemblies").setLevel(logging.DEBUG)
    suite = unittest.TestLoader().loadTestsFromTestCase(TestMMCIFParser)
    unittest.TextTestRunner(verbosity=2).run(suite)
