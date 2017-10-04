#!/local/bin/python
# -*- coding: utf-8 -*-


import logging
import sys
import unittest

try:
    from mock import patch
except ImportError:
    # python 3.5
    from unittest.mock import patch
from os import path

from proteofav.config import Defaults
from proteofav.structures import (parse_mmcif_atoms, _mmcif_fields, select_cif,
                                  parse_pdb_atoms, _fix_type_symbol,
                                  _fix_pdb_ins_code, _fix_label_alt_id)
from proteofav.utils import get_preferred_assembly_id

log = logging.getLogger(__name__)

defaults = Defaults("config.txt")


@patch("proteofav.structures.defaults", defaults)
class TestMMCIFParser(unittest.TestCase):
    """Test the mmCIF parser methods."""

    def setUp(self):
        """Initialize the framework for testing."""
        self.example_mmcif = path.join(path.dirname(__file__), "testdata",
                                       "mmcif/2pah.cif")
        self.example_mmcif_bio = path.join(path.dirname(__file__), "testdata",
                                           "mmcif/2pah_bio.cif")
        self.example_pdb = path.join(path.dirname(__file__), "testdata",
                                     "pdb/2pah.pdb")
        self.example_pdb2 = path.join(path.dirname(__file__), "testdata",
                                      "pdb/1ejg.pdb")
        self.mmcif_atom_parser = parse_mmcif_atoms
        self.mmcif_info_parser = _mmcif_fields
        self.example_tsv_out = path.join(path.dirname(__file__), "testdata",
                                         "mmcif/2pah-bio.tsv")
        self.select_cif = select_cif
        self.best_assembly = get_preferred_assembly_id
        self.pdbid = '2pah'
        self.pdb_atom_parser = parse_pdb_atoms
        self.fix_label_alt_id = _fix_label_alt_id
        self.fix_pdb_ins_code = _fix_pdb_ins_code
        self.fix_type_symbol = _fix_type_symbol

    def tearDown(self):
        """Remove testing framework."""

        self.example_mmcif = None
        self.example_mmcif_bio = None
        self.example_pdb = None
        self.example_pdb2 = None
        self.mmcif_atom_parser = None
        self.mmcif_info_parser = None
        self.example_tsv_out = None
        self.select_cif = None
        self.best_assembly = None
        self.pdbid = None
        self.pdb_atom_parser = None
        self.fix_label_alt_id = None
        self.fix_pdb_ins_code = None
        self.fix_type_symbol = None

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
        self.assertEqual(len(data.columns.values), 20)

        # check whether there are particular keys
        self.assertIn('label_asym_id', data.columns.values)

        # check the values for particular entries
        self.assertTrue(data.loc[1, 'label_asym_id'] == 'A')
        self.assertEqual(data.loc[1, 'pdbx_PDB_model_num'], '1')
        self.assertEqual(data['group_PDB'][0], 'ATOM')

    def test_info_to_table_mmcif(self):
        """
        Tests the parsing real mmCIF files.

        Some checks are made to whether the parsed keys and values
        are the ones we are expecting.

        Parsing all lines other than ATOM/HETATOM lines.
        """

        # describes possible macromolecular assemblies
        assembly = self.mmcif_info_parser(self.example_mmcif,
                                          field_name='_pdbx_struct_assembly.')
        # check some data values
        self.assertEqual(assembly.loc[0, 'details'], 'author_and_software_defined_assembly')
        self.assertEqual(assembly.loc[0, 'oligomeric_details'], 'tetrameric')
        self.assertEqual(assembly.loc[0, 'oligomeric_count'], 4)

        # details the generation of each macromolecular assembly
        assembly_gen = self.mmcif_info_parser(self.example_mmcif,
                                              field_name='_pdbx_struct_assembly_gen.')
        # check some data values
        self.assertEqual(assembly_gen.loc[0, 'asym_id_list'], 'A,C,B,D')

        # details translation and rotation operations required to generate/transform
        # assembly coordinates
        oper_list = self.mmcif_info_parser(self.example_mmcif,
                                           field_name='_pdbx_struct_oper_list.',
                                           require_index=True)
        # check some data values
        self.assertEqual(oper_list.loc[0, 'type'], 'identity operation')

    def test_bio_unit_to_table_mmcif(self):
        """
        Tests the parsing real mmCIF files.

        Some checks are made to whether the parsed keys and values
        are the ones we are expecting.

        Parses the atom lines and the info lines and looks for
        biological assemblies. Finds the most likely and generates
        a new pandas Dataframe with the new atoms.
        """

        # new approach that looks up the biological assemblies and
        # gets the preferred one from the PDBe server instead of
        # generating it locally

        # test loading the atom lines with biological assemblies
        best_assembly = self.best_assembly(self.pdbid)
        data = self.select_cif(self.pdbid, biounit=True,
                               assembly_id=best_assembly)

        # check some values
        # Biological assemblies from the PDBe come with updated chains
        # (auth_aym_id and label_aym_id) and two extra columns with the
        # original chain identifiers
        self.assertIn('orig_label_asym_id', data)
        self.assertIn('orig_auth_asym_id', data)

    def test_parser_cif(self):
        data = self.mmcif_atom_parser(self.example_mmcif)
        self.assertEqual(data.loc[0, 'label_asym_id'], 'A')
        self.assertEqual(data.loc[2686, 'label_asym_id'], 'B')
        self.assertEqual(data.loc[5315, 'label_asym_id'], 'C')
        self.assertEqual(data.loc[5316, 'label_asym_id'], 'D')
        self.assertEqual(data.loc[1, 'label_atom_id'], 'CA')
        self.assertEqual(data.loc[2687, 'label_atom_id'], 'CA')
        self.assertEqual(data.loc[5315, 'label_atom_id'], 'FE')
        self.assertEqual(data.loc[5316, 'label_atom_id'], 'FE')

    def test_parser_mmcif_bio(self):
        data = self.mmcif_atom_parser(self.example_mmcif_bio)
        self.assertEqual(data.loc[0, 'label_asym_id'], 'A')
        self.assertEqual(data.loc[5372, 'label_asym_id'], 'B')
        self.assertEqual(data.loc[10630, 'label_asym_id'], 'C')
        self.assertEqual(data.loc[10632, 'label_asym_id'], 'D')
        self.assertEqual(data.loc[2686, 'label_asym_id'], 'AA')
        self.assertEqual(data.loc[8001, 'label_asym_id'], 'BA')
        self.assertEqual(data.loc[10631, 'label_asym_id'], 'CA')
        self.assertEqual(data.loc[10633, 'label_asym_id'], 'DA')

    def test_parser_pdb(self):
        data = self.pdb_atom_parser(self.example_pdb)
        self.assertEqual(data.loc[0, 'label_asym_id'], 'A')
        self.assertEqual(data.loc[2686, 'label_asym_id'], 'B')
        self.assertEqual(data.loc[5315, 'label_asym_id'], 'A')
        self.assertEqual(data.loc[5316, 'label_asym_id'], 'B')
        self.assertEqual(data.loc[1, 'label_atom_id'], 'CA')

    def test_fix_pdb_ins_code(self):
        data = self.pdb_atom_parser(self.example_pdb2, fix_ins_code=False)
        self.assertEqual(1, data.loc[0, 'id'])
        self.assertEqual('', data.loc[0, 'pdbx_PDB_ins_code'])
        data = self.fix_pdb_ins_code(data)
        self.assertEqual(1, data.loc[0, 'id'])
        self.assertEqual('?', data.loc[0, 'pdbx_PDB_ins_code'])

    def test_fix_label_alt_id(self):
        data = self.pdb_atom_parser(self.example_pdb2, fix_label_alt_id=False)
        self.assertEqual(1, data.loc[0, 'id'])
        self.assertEqual('N', data.loc[0, 'label_atom_id'])
        self.assertEqual('A', data.loc[0, 'label_alt_id'])
        data = self.fix_label_alt_id(data)
        self.assertEqual(1, data.loc[0, 'id'])
        self.assertEqual('N', data.loc[0, 'label_atom_id'])
        self.assertEqual('A', data.loc[0, 'label_alt_id'])

    def test_fix_type_symbol(self):
        data = self.pdb_atom_parser(self.example_pdb2, fix_type_symbol=False)
        self.assertEqual(1, data.loc[0, 'id'])
        self.assertEqual('', data.loc[0, 'type_symbol'])
        data = self.fix_type_symbol(data)
        self.assertEqual(1, data.loc[0, 'id'])
        self.assertEqual('N', data.loc[0, 'type_symbol'])


if __name__ == '__main__':
    logging.basicConfig(stream=sys.stderr)
    logging.getLogger("mmCIF related methods").setLevel(logging.DEBUG)
    suite = unittest.TestLoader().loadTestsFromTestCase(TestMMCIFParser)
    unittest.TextTestRunner(verbosity=2).run(suite)
