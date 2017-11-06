# -*- coding: utf-8 -*-

import os
import sys
import logging
import unittest
import pandas as pd

try:
    from mock import patch
except ImportError:
    # python 3.5
    from unittest.mock import patch

from proteofav.config import defaults
from proteofav.structures import (parse_mmcif_atoms, _mmcif_fields, select_structures,
                                  parse_pdb_atoms, _fix_type_symbol,
                                  _fix_pdb_ins_code, _fix_label_alt_id,
                                  write_mmcif_from_table, write_pdb_from_table,
                                  _get_atom_line, residues_aggregation,
                                  filter_structures,
                                  _add_mmcif_res_full, _add_mmcif_atom_altloc,
                                  _remove_multiple_altlocs, _remove_partial_residues,
                                  read_structures, write_structures, download_structures,
                                  PDB, mmCIF, get_preferred_assembly_id,
                                  fetch_summary_properties_pdbe, get_sequence)


@patch("proteofav.structures.defaults", defaults)
class TestMMCIFParser(unittest.TestCase):
    """Test the mmCIF parser methods."""

    def setUp(self):
        """Initialize the framework for testing."""
        self.example_mmcif = os.path.join(os.path.dirname(__file__), "testdata",
                                          "mmcif", "2pah.cif")
        self.example_mmcif_bio = os.path.join(os.path.dirname(__file__), "testdata",
                                              "mmcif", "2pah_bio.cif")
        self.example_pdb = os.path.join(os.path.dirname(__file__), "testdata",
                                        "pdb", "2pah.pdb")
        self.example_pdb2 = os.path.join(os.path.dirname(__file__), "testdata",
                                         "pdb", "1ejg.pdb")
        self.output_mmcif = os.path.join(os.path.dirname(__file__), "testdata",
                                         "2pah.cif")
        self.output_pdb = os.path.join(os.path.dirname(__file__), "testdata",
                                       "2pah.pdb")
        self.mmcif_atom_parser = parse_mmcif_atoms
        self.mmcif_info_parser = _mmcif_fields
        self.example_tsv_out = os.path.join(os.path.dirname(__file__), "testdata",
                                            "mmcif", "2pah-bio.tsv")
        self.select_structures = select_structures
        self.fetch_summary_properties_pdbe = fetch_summary_properties_pdbe
        self.get_preferred_assembly_id = get_preferred_assembly_id
        self.pdbid = '2pah'
        self.pdb_atom_parser = parse_pdb_atoms
        self.fix_label_alt_id = _fix_label_alt_id
        self.fix_pdb_ins_code = _fix_pdb_ins_code
        self.fix_type_symbol = _fix_type_symbol
        self.write_mmcif_from_table = write_mmcif_from_table
        self.write_pdb_from_table = write_pdb_from_table
        self.get_atom_line = _get_atom_line
        self.residues_aggregation = residues_aggregation
        self.filter_structures = filter_structures
        self.add_mmcif_res_full = _add_mmcif_res_full
        self.add_mmcif_atom_altloc = _add_mmcif_atom_altloc
        self.remove_multiple_altlocs = _remove_multiple_altlocs
        self.remove_partial_residues = _remove_partial_residues
        self.read_structures = read_structures
        self.write_structures = write_structures
        self.download_structures = download_structures
        self.PDB = PDB
        self.mmCIF = mmCIF
        self.get_sequence = get_sequence

    def tearDown(self):
        """Remove testing framework."""

        self.example_mmcif = None
        self.example_mmcif_bio = None
        self.example_pdb = None
        self.example_pdb2 = None
        self.output_mmcif = None
        self.output_pdb = None
        self.mmcif_atom_parser = None
        self.mmcif_info_parser = None
        self.example_tsv_out = None
        self.select_structures = None
        self.fetch_summary_properties_pdbe = None
        self.get_preferred_assembly_id = None
        self.pdbid = None
        self.pdb_atom_parser = None
        self.fix_label_alt_id = None
        self.fix_pdb_ins_code = None
        self.fix_type_symbol = None
        self.write_mmcif_from_table = None
        self.write_pdb_from_table = None
        self.get_atom_line = None
        self.residues_aggregation = None
        self.filter_structures = None
        self.add_mmcif_res_full = None
        self.add_mmcif_atom_altloc = None
        self.remove_multiple_altlocs = None
        self.remove_partial_residues = None
        self.read_structures = None
        self.write_structures = None
        self.download_structures = None
        self.PDB = None
        self.mmCIF = None
        self.get_sequence = None

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
        best_assembly = self.get_preferred_assembly_id(self.pdbid)
        data = self.select_structures(self.pdbid, bio_unit=True,
                                      bio_unit_id=best_assembly)

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

    def test_writer_cif(self):
        data = self.mmcif_atom_parser(self.example_mmcif)
        self.write_mmcif_from_table(table=data, filename=self.output_mmcif,
                                    overwrite=True)
        self.assertTrue(os.path.isfile(self.output_mmcif))
        data = self.mmcif_atom_parser(self.output_mmcif)
        self.assertIn('label_asym_id', list(data))
        os.remove(self.output_mmcif)

    def test_writer_biocif(self):
        data = self.mmcif_atom_parser(self.example_mmcif_bio)
        self.write_mmcif_from_table(table=data, filename=self.output_mmcif,
                                    overwrite=True)
        self.assertTrue(os.path.isfile(self.output_mmcif))
        data = self.mmcif_atom_parser(self.output_mmcif)
        self.assertIn('label_asym_id', list(data))
        os.remove(self.output_mmcif)

    def test_writer_mmcif2pdb(self):
        data = self.mmcif_atom_parser(self.example_mmcif)
        self.write_pdb_from_table(table=data, filename=self.output_pdb,
                                  overwrite=True)
        self.assertTrue(os.path.isfile(self.output_pdb))
        data = self.pdb_atom_parser(self.output_pdb)
        self.assertIn('label_asym_id', list(data))
        os.remove(self.output_pdb)

    def test_writer_pdb2mmcif(self):
        data = self.pdb_atom_parser(self.example_pdb)
        self.write_mmcif_from_table(table=data, filename=self.output_mmcif,
                                    overwrite=True)
        self.assertTrue(os.path.isfile(self.output_mmcif))
        data = self.mmcif_atom_parser(self.output_mmcif)
        self.assertIn('label_asym_id', list(data))
        os.remove(self.output_mmcif)

    def test_get_atom_line(self):
        data = self.mmcif_atom_parser(self.example_mmcif)
        line = self.get_atom_line(data, index=0, atom_number=1)
        r = 'ATOM      1  N   VAL A 118      -7.069  21.943  18.770  1.00 56.51           N  \n'
        self.assertEqual(str(line), r)
        line = self.get_atom_line(data, index=1000, atom_number=1000)
        r = 'ATOM   1000  CD  ARG A 241       1.614   7.798  40.365  1.00 23.29           C  \n'
        self.assertEqual(str(line), r)

    def test_residues_aggregation_first(self):
        data = self.mmcif_atom_parser(self.example_mmcif)
        self.assertEqual(1, data.loc[0, 'id'])
        self.assertEqual('N', data.loc[0, 'type_symbol'])
        self.assertEqual('N', data.loc[0, 'label_atom_id'])
        self.assertEqual('VAL', data.loc[0, 'label_comp_id'])
        self.assertEqual('118', data.loc[0, 'auth_seq_id'])
        self.assertEqual(-7.069, data.loc[0, 'Cartn_x'])
        data = self.residues_aggregation(data, agg_method='first', category='auth')
        self.assertEqual(1, data.loc[0, 'id'])
        self.assertEqual('N', data.loc[0, 'type_symbol'])
        self.assertEqual('N', data.loc[0, 'label_atom_id'])
        self.assertEqual('VAL', data.loc[0, 'label_comp_id'])
        self.assertEqual('118', data.loc[0, 'auth_seq_id'])
        self.assertEqual(-7.069, data.loc[0, 'Cartn_x'])

    def test_residues_aggregation_centroid(self):
        data = self.mmcif_atom_parser(self.example_mmcif)
        self.assertEqual(1, data.loc[0, 'id'])
        self.assertEqual('N', data.loc[0, 'type_symbol'])
        self.assertEqual('N', data.loc[0, 'label_atom_id'])
        self.assertEqual('VAL', data.loc[0, 'label_comp_id'])
        self.assertEqual('118', data.loc[0, 'auth_seq_id'])
        self.assertEqual(-7.069, data.loc[0, 'Cartn_x'])
        data = self.residues_aggregation(data, agg_method='centroid', category='auth')
        self.assertEqual(1, data.loc[0, 'id'])
        self.assertEqual('N', data.loc[0, 'type_symbol'])
        self.assertEqual('N', data.loc[0, 'label_atom_id'])
        self.assertEqual('VAL', data.loc[0, 'label_comp_id'])
        self.assertEqual('118', data.loc[0, 'auth_seq_id'])
        self.assertAlmostEqual(-7.310, data.loc[0, 'Cartn_x'], places=2)

    def test_residues_aggregation_backbone_centroid(self):
        data = self.mmcif_atom_parser(self.example_mmcif)
        self.assertEqual(1, data.loc[0, 'id'])
        self.assertEqual('N', data.loc[0, 'type_symbol'])
        self.assertEqual('N', data.loc[0, 'label_atom_id'])
        self.assertEqual('VAL', data.loc[0, 'label_comp_id'])
        self.assertEqual('118', data.loc[0, 'auth_seq_id'])
        self.assertEqual(-7.069, data.loc[0, 'Cartn_x'])
        data = self.residues_aggregation(data, agg_method='backbone_centroid', category='auth')
        self.assertEqual(1, data.loc[0, 'id'])
        self.assertEqual('N', data.loc[0, 'type_symbol'])
        self.assertEqual('N', data.loc[0, 'label_atom_id'])
        self.assertEqual('VAL', data.loc[0, 'label_comp_id'])
        self.assertEqual('118', data.loc[0, 'auth_seq_id'])
        self.assertAlmostEqual(-6.312, data.loc[0, 'Cartn_x'], places=2)

    def test_residues_aggregation_mean(self):
        data = self.mmcif_atom_parser(self.example_mmcif)
        self.assertEqual(1, data.loc[0, 'id'])
        self.assertEqual('N', data.loc[0, 'type_symbol'])
        self.assertEqual('N', data.loc[0, 'label_atom_id'])
        self.assertEqual('VAL', data.loc[0, 'label_comp_id'])
        self.assertEqual('118', data.loc[0, 'auth_seq_id'])
        self.assertEqual(-7.069, data.loc[0, 'Cartn_x'])
        data = self.residues_aggregation(data, agg_method='backbone_centroid', category='auth')
        self.assertEqual(1, data.loc[0, 'id'])
        self.assertEqual('N', data.loc[0, 'type_symbol'])
        self.assertEqual('N', data.loc[0, 'label_atom_id'])
        self.assertEqual('VAL', data.loc[0, 'label_comp_id'])
        self.assertEqual('118', data.loc[0, 'auth_seq_id'])
        self.assertAlmostEqual(-6.312, data.loc[0, 'Cartn_x'], places=2)

    def test_filter_structures_residues_aggregation_centroid(self):
        data = self.mmcif_atom_parser(self.example_mmcif)
        self.assertEqual(1, data.loc[0, 'id'])
        self.assertEqual('N', data.loc[0, 'type_symbol'])
        self.assertEqual('N', data.loc[0, 'label_atom_id'])
        self.assertEqual('VAL', data.loc[0, 'label_comp_id'])
        self.assertEqual('118', data.loc[0, 'auth_seq_id'])
        self.assertEqual(-7.069, data.loc[0, 'Cartn_x'])
        data = self.filter_structures(data, residue_agg=True, agg_method='centroid')
        self.assertEqual(1, data.loc[0, 'id'])
        self.assertEqual('N', data.loc[0, 'type_symbol'])
        self.assertEqual('N', data.loc[0, 'label_atom_id'])
        self.assertEqual('VAL', data.loc[0, 'label_comp_id'])
        self.assertEqual('118', data.loc[0, 'auth_seq_id'])
        self.assertAlmostEqual(-7.310, data.loc[0, 'Cartn_x'], places=2)

    def test_filter_structures_excluded_cols(self):
        data = self.mmcif_atom_parser(self.example_mmcif, excluded_cols=())
        keys = [k for k in data]
        self.assertIn("Cartn_x_esd", keys)
        self.assertIn("Cartn_y_esd", keys)
        self.assertIn("Cartn_z_esd", keys)
        self.assertIn("occupancy_esd", keys)
        self.assertIn("B_iso_or_equiv_esd", keys)
        exc_cols = ('Cartn_x_esd', 'Cartn_y_esd', 'Cartn_z_esd',
                    'occupancy_esd', 'B_iso_or_equiv_esd')
        data = self.filter_structures(data, excluded_cols=exc_cols)
        keys = [k for k in data]
        self.assertNotIn("Cartn_x_esd", keys)
        self.assertNotIn("Cartn_y_esd", keys)
        self.assertNotIn("Cartn_z_esd", keys)
        self.assertNotIn("occupancy_esd", keys)
        self.assertNotIn("B_iso_or_equiv_esd", keys)

    def test_filter_structures_chain_id(self):
        data = self.mmcif_atom_parser(self.example_mmcif)
        data = self.filter_structures(data, chains=('A',))
        self.assertIn("A", data.label_asym_id.unique())
        self.assertNotIn("B", data.label_asym_id.unique())

    def test_filter_structures_atom_lines(self):
        data = self.mmcif_atom_parser(self.example_mmcif)
        data = self.filter_structures(data, lines=('ATOM',))
        self.assertIn("ATOM", data.group_PDB.unique())
        self.assertNotIn("HETATM", data.group_PDB.unique())

    def test_filter_structures_comps(self):
        data = self.mmcif_atom_parser(self.example_mmcif)
        data = self.filter_structures(data, comps=('PHE',))
        self.assertIn("PHE", data.label_comp_id.unique())
        self.assertNotIn("TYR", data.label_comp_id.unique())

    def test_filter_structures_res(self):
        data = self.mmcif_atom_parser(self.example_mmcif)
        data = self.filter_structures(data, res=('285',))
        self.assertIn("285", data.auth_seq_id.unique())
        self.assertNotIn("289", data.auth_seq_id.unique())

    def test_filter_structures_atoms(self):
        data = self.mmcif_atom_parser(self.example_mmcif)
        data = self.filter_structures(data, atoms=('CA',))
        self.assertIn("CA", data.label_atom_id.unique())
        self.assertNotIn("N", data.label_atom_id.unique())

    def test_filter_structures_remove_hydrogen(self):
        data = self.pdb_atom_parser(self.example_pdb2)
        self.assertIn("H", data.label_atom_id.unique())
        data = self.filter_structures(data, remove_hydrogens=True)
        self.assertNotIn("H", data.label_atom_id.unique())

    def test_filter_structures_add_res_full(self):
        data = self.pdb_atom_parser(self.example_pdb2)
        data = self.add_mmcif_res_full(data)
        self.assertIn('label_seq_id_full', data)
        self.assertIn('auth_seq_id_full', data)
        self.assertEqual(data.loc[0, 'auth_seq_id_full'], '1')
        data = self.pdb_atom_parser(self.example_pdb2)
        data = self.filter_structures(data, add_res_full=True)
        self.assertEqual(data.loc[0, 'auth_seq_id_full'], '1')

    def test_filter_structures_add_pdb_altloc(self):
        data = self.pdb_atom_parser(self.example_pdb2)
        data = self.add_mmcif_atom_altloc(data)
        self.assertIn('label_atom_altloc_id', [k for k in data])
        self.assertIn('auth_atom_altloc_id', [k for k in data])
        self.assertEqual('N.A', data.loc[0, 'label_atom_altloc_id'])
        self.assertEqual('N.B', data.loc[1, 'label_atom_altloc_id'])
        self.assertEqual('CA.A', data.loc[2, 'label_atom_altloc_id'])
        self.assertEqual('CA.B', data.loc[3, 'label_atom_altloc_id'])
        self.assertEqual('C', data.loc[4, 'label_atom_altloc_id'])
        self.assertEqual('O', data.loc[5, 'label_atom_altloc_id'])
        data = self.pdb_atom_parser(self.example_pdb2)
        data = self.filter_structures(data, remove_altloc=False, add_atom_altloc=True)
        self.assertIn('label_atom_altloc_id', [k for k in data])
        self.assertIn('auth_atom_altloc_id', [k for k in data])
        self.assertEqual('N.A', data.loc[0, 'label_atom_altloc_id'])
        self.assertEqual('N.B', data.loc[1, 'label_atom_altloc_id'])
        self.assertEqual('CA.A', data.loc[2, 'label_atom_altloc_id'])
        self.assertEqual('CA.B', data.loc[3, 'label_atom_altloc_id'])
        self.assertEqual('C', data.loc[4, 'label_atom_altloc_id'])
        self.assertEqual('O', data.loc[5, 'label_atom_altloc_id'])

    def test_filter_structures_remove_altloc(self):
        data = self.pdb_atom_parser(self.example_pdb2)
        data = self.remove_multiple_altlocs(data)
        self.assertEqual(1, data.loc[0, 'id'])
        self.assertEqual('N', data.loc[0, 'type_symbol'])
        self.assertEqual('N', data.loc[0, 'label_atom_id'])
        self.assertEqual('.', data.loc[0, 'label_alt_id'])
        data = self.pdb_atom_parser(self.example_pdb2)
        data = self.filter_structures(data, remove_altloc=True)
        self.assertEqual(1, data.loc[0, 'id'])
        self.assertEqual('N', data.loc[0, 'type_symbol'])
        self.assertEqual('N', data.loc[0, 'label_atom_id'])
        self.assertEqual('.', data.loc[0, 'label_alt_id'])

    def test_filter_structures_remove_partial_res(self):
        data = self.pdb_atom_parser(self.example_pdb2)
        # counting the number of atoms in res with label_seq_id = '25'
        self.assertEqual(54, list(data.label_seq_id.tolist()).count('25'))
        data = self.remove_partial_residues(data)
        self.assertEqual(19, list(data.label_seq_id.tolist()).count('25'))
        data = self.pdb_atom_parser(self.example_pdb2)
        self.assertEqual(54, list(data.label_seq_id.tolist()).count('25'))
        data = filter_structures(data, models=None, remove_partial_res=True,
                                 remove_altloc=False, remove_hydrogens=False,
                                 reset_atom_id=False, add_res_full=False)
        self.assertEqual(19, list(data.label_seq_id.tolist()).count('25'))

    def test_read_structures(self):
        data = self.read_structures(self.example_mmcif)
        self.assertEqual(data.loc[0, 'label_asym_id'], 'A')
        self.assertEqual(data.loc[2686, 'label_asym_id'], 'B')
        self.assertEqual(data.loc[5315, 'label_asym_id'], 'C')
        self.assertEqual(data.loc[5316, 'label_asym_id'], 'D')
        self.assertEqual(data.loc[1, 'label_atom_id'], 'CA')
        self.assertEqual(data.loc[2687, 'label_atom_id'], 'CA')
        self.assertEqual(data.loc[5315, 'label_atom_id'], 'FE')
        self.assertEqual(data.loc[5316, 'label_atom_id'], 'FE')

        data = self.read_structures(self.example_pdb)
        self.assertEqual(data.loc[0, 'label_asym_id'], 'A')
        self.assertEqual(data.loc[2686, 'label_asym_id'], 'B')
        self.assertEqual(data.loc[5315, 'label_asym_id'], 'A')
        self.assertEqual(data.loc[5316, 'label_asym_id'], 'B')
        self.assertEqual(data.loc[1, 'label_atom_id'], 'CA')

        example_dssp = os.path.join(os.path.dirname(__file__), "testdata",
                                    "dssp", "2pah.dssp")
        with self.assertRaises(ValueError, msg='DSSP format is not recognised'):
            self.read_structures(example_dssp)

    def test_write_structures(self):
        # read mmcif and write pdb
        data = self.mmcif_atom_parser(self.example_mmcif)
        self.write_structures(data, filename=self.output_pdb, overwrite=True)
        self.assertTrue(os.path.isfile(self.output_pdb))
        data = self.pdb_atom_parser(self.output_pdb)
        self.assertIn('label_asym_id', list(data))
        os.remove(self.output_pdb)

        # read pdb and write mmcif
        data = self.pdb_atom_parser(self.example_pdb)
        self.write_structures(data, filename=self.output_mmcif, overwrite=True)
        self.assertTrue(os.path.isfile(self.output_mmcif))
        data = self.mmcif_atom_parser(self.output_mmcif)
        self.assertIn('label_asym_id', list(data))
        os.remove(self.output_mmcif)

    def test_download_structures(self):
        self.download_structures(self.pdbid, self.output_mmcif,
                                 overwrite=True)
        self.assertTrue(os.path.isfile(self.output_mmcif))
        os.remove(self.output_mmcif)

        self.download_structures(self.pdbid, self.output_mmcif, bio_unit=True,
                                 overwrite=True)
        self.assertTrue(os.path.isfile(self.output_mmcif))
        os.remove(self.output_mmcif)

        self.download_structures(self.pdbid, self.output_pdb,
                                 overwrite=True)
        self.assertTrue(os.path.isfile(self.output_pdb))
        os.remove(self.output_pdb)

    def test_main_mmCIF(self):
        # read
        data = self.mmCIF.read(self.example_mmcif)
        self.assertEqual(data.loc[0, 'label_asym_id'], 'A')
        self.assertEqual(data.loc[2686, 'label_asym_id'], 'B')
        self.assertEqual(data.loc[5315, 'label_asym_id'], 'C')
        self.assertEqual(data.loc[5316, 'label_asym_id'], 'D')
        self.assertEqual(data.loc[1, 'label_atom_id'], 'CA')
        self.assertEqual(data.loc[2687, 'label_atom_id'], 'CA')
        self.assertEqual(data.loc[5315, 'label_atom_id'], 'FE')
        self.assertEqual(data.loc[5316, 'label_atom_id'], 'FE')
        # read mmcif and write pdb
        data = self.mmCIF.read(self.example_mmcif)
        self.mmCIF.write(data, filename=self.output_pdb, overwrite=True)
        self.assertTrue(os.path.isfile(self.output_pdb))
        data = self.pdb_atom_parser(self.output_pdb)
        self.assertIn('label_asym_id', list(data))
        os.remove(self.output_pdb)
        # download
        self.mmCIF.download(self.pdbid, self.output_mmcif,
                            overwrite=True)
        self.assertTrue(os.path.isfile(self.output_mmcif))
        os.remove(self.output_mmcif)
        # select
        self.mmCIF.select(self.pdbid)

    def test_main_PDB(self):
        # read
        data = self.PDB.read(self.example_pdb)
        self.assertEqual(data.loc[0, 'label_asym_id'], 'A')
        self.assertEqual(data.loc[2686, 'label_asym_id'], 'B')
        self.assertEqual(data.loc[5315, 'label_asym_id'], 'A')
        self.assertEqual(data.loc[5316, 'label_asym_id'], 'B')
        self.assertEqual(data.loc[1, 'label_atom_id'], 'CA')
        # read pdb and write mmcif
        data = self.PDB.read(self.example_pdb)
        self.PDB.write(data, filename=self.output_mmcif, overwrite=True)
        self.assertTrue(os.path.isfile(self.output_mmcif))
        data = self.mmcif_atom_parser(self.output_mmcif)
        self.assertIn('label_asym_id', list(data))
        os.remove(self.output_mmcif)
        # download
        self.PDB.download(self.pdbid, self.output_pdb,
                          overwrite=True)
        self.assertTrue(os.path.isfile(self.output_pdb))
        os.remove(self.output_pdb)
        # select
        self.PDB.select(self.pdbid)

    def test_summary_properties_pdbe(self):
        r = self.fetch_summary_properties_pdbe(self.pdbid)
        self.assertTrue(r.ok)

    def test_preferred_assembly_pdbe(self):
        r = self.get_preferred_assembly_id(self.pdbid)
        self.assertEqual("1", r)

    def test_get_sequence_structures(self):
        table = mmCIF.read(filename=self.example_mmcif)
        self.assertTrue(isinstance(table, pd.DataFrame))
        table = self.filter_structures(table, chains=('A',), lines=('ATOM',))
        table = self.residues_aggregation(table)
        seq = self.get_sequence(table)
        self.assertEqual(329, len(seq))
        self.assertEqual('VPWFPRTIQELDRFANQILDADHPG', seq[0:25])


if __name__ == '__main__':
    logging.basicConfig(stream=sys.stderr)
    logging.getLogger("proteofav.config").setLevel(logging.CRITICAL)
    suite = unittest.TestLoader().loadTestsFromTestCase(TestMMCIFParser)
    unittest.TextTestRunner(verbosity=2).run(suite)
