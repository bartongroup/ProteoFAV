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
from proteofav.structures import _parse_mmcif_atoms_from_file, _mmcif_fields, select_cif
from proteofav.utils import get_preferred_assembly_id

log = logging.getLogger(__name__)

defaults = Defaults("config.txt")


@patch("proteofav.structures.defaults", defaults)
class TestMMCIFParser(unittest.TestCase):
    """Test the mmCIF parser methods."""

    def setUp(self):
        """Initialize the framework for testing."""
        self.example_mmcif = path.join(path.dirname(__file__), "CIF/2pah.cif")
        self.mmcif_atom_parser = _parse_mmcif_atoms_from_file
        self.mmcif_info_parser = _mmcif_fields
        self.example_tsv_out = path.join(path.dirname(__file__), "CIF/2pah-bio.tsv")
        self.select_cif = select_cif
        self.best_assembly = get_preferred_assembly_id
        self.pdbid = '2pah'

    def tearDown(self):
        """Remove testing framework."""

        self.example_mmcif = None
        self.mmcif_atom_parser = None
        self.mmcif_info_parser = None
        self.example_tsv_out = None
        self.select_cif = None
        self.best_assembly = None
        self.pdbid = None

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
        self.assertEqual(len(data.columns.values), 21)

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


if __name__ == '__main__':
    logging.basicConfig(stream=sys.stderr)
    logging.getLogger("mmCIF related methods").setLevel(logging.DEBUG)
    suite = unittest.TestLoader().loadTestsFromTestCase(TestMMCIFParser)
    unittest.TextTestRunner(verbosity=2).run(suite)
