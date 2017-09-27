# -*- coding: utf-8 -*-

import os
import sys
import logging
import unittest
import pandas as pd

from proteofav.structures import (PDB, mmCIF,
                                  _fix_label_alt_id,
                                  _fix_pdb_ins_code,
                                  _fix_type_symbol,
                                  _add_mmcif_res_full,
                                  _add_mmcif_atom_altloc,
                                  _remove_multiple_altlocs,
                                  _remove_partial_residues,
                                  residues_aggregation,
                                  write_mmcif_from_table,
                                  _get_atom_line,
                                  write_pdb_from_table)

from proteofav.config import defaults as config


class TestStructures(unittest.TestCase):
    """Test the PDBx parser methods."""

    def setUp(self):
        """Initialize the framework for testing."""

        root = os.path.abspath(os.path.dirname(__file__))

        self.pdbid = '2pah'
        self.pdbid2 = '1ejg'
        self.inputcif = os.path.join(root, "testdata", config.db_mmcif,
                                     "{}.cif".format(self.pdbid))
        self.inputbiocif = os.path.join(root, "testdata", config.db_mmcif,
                                        "{}_bio.cif".format(self.pdbid))
        self.outputcif = os.path.join(root, "testdata", config.db_tmp,
                                      "{}.cif".format(self.pdbid))
        self.inputpdb = os.path.join(root, "testdata", config.db_pdb,
                                     "{}.pdb".format(self.pdbid))
        self.inputpdb2 = os.path.join(root, "testdata", config.db_pdb,
                                      "{}.pdb".format(self.pdbid2))
        self.outputpdb = os.path.join(root, "testdata", config.db_tmp,
                                      "{}.pdb".format(self.pdbid))
        self.excluded = ()
        self.pdb = PDB
        self.mmcif = mmCIF
        self.fix_label_alt_id = _fix_label_alt_id
        self.fix_pdb_ins_code = _fix_pdb_ins_code
        self.fix_type_symbol = _fix_type_symbol
        self.add_res_full = _add_mmcif_res_full
        self.add_mmcif_atom_altloc = _add_mmcif_atom_altloc
        self.remove_altloc = _remove_multiple_altlocs
        self.remove_partial_residues = _remove_partial_residues
        self.residues_aggregation = residues_aggregation
        self.write_mmcif_from_table = write_mmcif_from_table
        self.get_atom_line = _get_atom_line
        self.write_pdb_from_table = write_pdb_from_table

        logging.disable(logging.DEBUG)

    def tearDown(self):
        """Remove testing framework."""

        self.pdbid = None
        self.pdbid2 = None
        self.inputcif = None
        self.inputbiocif = None
        self.outputcif = None
        self.inputpdb = None
        self.inputpdb2 = None
        self.outputpdb = None
        self.excluded = None
        self.pdb = None
        self.mmcif = None
        self.fix_label_alt_id = None
        self.fix_pdb_ins_code = None
        self.fix_type_symbol = None
        self.add_res_full = None
        self.add_mmcif_atom_altloc = None
        self.remove_altloc = None
        self.remove_partial_residues = None
        self.residues_aggregation = None
        self.write_mmcif_from_table = None
        self.get_atom_line = None
        self.write_pdb_from_table = None

        logging.disable(logging.NOTSET)

    def test_reader_cif(self):
        data = self.mmcif().read(self.inputcif)
        self.assertEqual(data.loc[0, 'label_asym_id'], 'A')
        self.assertEqual(data.loc[2686, 'label_asym_id'], 'B')
        self.assertEqual(data.loc[5315, 'label_asym_id'], 'C')
        self.assertEqual(data.loc[5316, 'label_asym_id'], 'D')
        self.assertEqual(data.loc[1, 'label_atom_id'], 'CA')
        self.assertEqual(data.loc[2687, 'label_atom_id'], 'CA')
        self.assertEqual(data.loc[5315, 'label_atom_id'], 'FE')
        self.assertEqual(data.loc[5316, 'label_atom_id'], 'FE')

    def test_reader_biocif(self):
        data = self.mmcif().read(self.inputbiocif)
        self.assertEqual(data.loc[0, 'label_asym_id'], 'A')
        self.assertEqual(data.loc[5372, 'label_asym_id'], 'B')
        self.assertEqual(data.loc[10630, 'label_asym_id'], 'C')
        self.assertEqual(data.loc[10632, 'label_asym_id'], 'D')
        self.assertEqual(data.loc[2686, 'label_asym_id'], 'AA')
        self.assertEqual(data.loc[8001, 'label_asym_id'], 'BA')
        self.assertEqual(data.loc[10631, 'label_asym_id'], 'CA')
        self.assertEqual(data.loc[10633, 'label_asym_id'], 'DA')

    def test_reader_pdb(self):
        data = self.pdb().read(self.inputpdb)
        self.assertEqual(data.loc[0, 'label_asym_id'], 'A')
        self.assertEqual(data.loc[2686, 'label_asym_id'], 'B')
        self.assertEqual(data.loc[5315, 'label_asym_id'], 'A')
        self.assertEqual(data.loc[5316, 'label_asym_id'], 'B')
        self.assertEqual(data.loc[1, 'label_atom_id'], 'CA')
        self.assertEqual(data.loc[2687, 'label_atom_id'], 'CA')
        self.assertEqual(data.loc[5315, 'label_atom_id'], 'FE')
        self.assertEqual(data.loc[5316, 'label_atom_id'], 'FE')

    def test_reader_mmcif_default_excluded(self):
        data = self.mmcif().read(self.inputcif)
        keys = [k for k in data]
        self.assertNotIn("Cartn_x_esd", keys)
        self.assertNotIn("Cartn_y_esd", keys)
        self.assertNotIn("Cartn_z_esd", keys)
        self.assertNotIn("occupancy_esd", keys)
        self.assertNotIn("B_iso_or_equiv_esd", keys)

    def test_reader_mmcif_new_excluded(self):
        data = self.mmcif().read(self.inputcif, excluded_cols=self.excluded)
        keys = [k for k in data]
        self.assertIn("Cartn_x_esd", keys)
        self.assertIn("Cartn_y_esd", keys)
        self.assertIn("Cartn_z_esd", keys)
        self.assertIn("occupancy_esd", keys)
        self.assertIn("B_iso_or_equiv_esd", keys)

    def test_reader_mmcif_default_category(self):
        data = self.mmcif().read(self.inputcif, category="label")
        keys = [k for k in data.label_asym_id.unique()]
        self.assertEqual(keys, ['A', 'B', 'C', 'D'])

    def test_reader_mmcif_new_category(self):
        data = self.mmcif().read(self.inputcif, category="auth")
        keys = [k for k in data.auth_asym_id.unique()]
        self.assertEqual(keys, ['A', 'B'])

    def test_filter_mmcif_chain_id(self):
        data = self.mmcif().read(self.inputcif, chains=('A',))
        self.assertIn("A", data.label_asym_id.unique())
        self.assertNotIn("B", data.label_asym_id.unique())

    def test_filter_atom_lines(self):
        data = self.mmcif().read(self.inputcif, lines=('ATOM',))
        self.assertIn("ATOM", data.group_PDB.unique())
        self.assertNotIn("HETATM", data.group_PDB.unique())

    def test_writer_cif(self):
        data = self.mmcif().read(self.inputcif)
        self.mmcif().write(table=data, filename=self.outputcif, overwrite=True)
        self.assertTrue(os.path.isfile(self.outputcif))
        os.remove(self.outputcif)

    def test_writer_biocif(self):
        data = self.mmcif().read(self.inputbiocif)
        self.mmcif().write(table=data, filename=self.outputcif, overwrite=True)
        self.assertTrue(os.path.isfile(self.outputcif))
        os.remove(self.outputcif)

    def test_writer_mmcif2pdb(self):
        data = self.mmcif().read(self.inputcif)
        self.pdb().write(table=data, filename=self.outputpdb, overwrite=True)
        self.assertTrue(os.path.isfile(self.outputpdb))
        os.remove(self.outputpdb)

    def test_writer_pdb2mmcif(self):
        data = self.pdb().read(self.inputpdb)
        self.pdb().write(table=data, filename=self.outputcif, overwrite=True)
        self.assertTrue(os.path.isfile(self.outputcif))
        os.remove(self.outputcif)

    def test_get_atom_line(self):
        data = self.mmcif().read(self.inputcif)
        line = self.get_atom_line(data, index=0, atom_number=1)
        r = 'ATOM      1  N   VAL A 118      -7.069  21.943  18.770  1.00 56.51           N  \n'
        self.assertEqual(str(line), r)
        line = self.get_atom_line(data, index=1000, atom_number=1000)
        r = 'ATOM   1000  CD  ARG A 241       1.614   7.798  40.365  1.00 23.29           C  \n'
        self.assertEqual(str(line), r)

    def test_reader_add_res_full(self):
        data = self.mmcif().read(self.inputcif, add_res_full=True)
        self.assertIn('label_seq_id_full', data)
        self.assertIn('auth_seq_id_full', data)
        self.assertEqual(data.loc[0, 'auth_seq_id_full'], '118')

    def test_add_res_full(self):
        data = self.mmcif().read(self.inputcif)
        data = self.add_res_full(data)
        self.assertIn('label_seq_id_full', data)
        self.assertIn('auth_seq_id_full', data)
        self.assertEqual(data.loc[0, 'auth_seq_id_full'], '118')

    def test_residues_aggregation_first(self):
        data = self.mmcif().read(self.inputcif)
        self.assertEqual(1, data.loc[0, 'id'])
        self.assertEqual('N', data.loc[0, 'type_symbol'])
        self.assertEqual('N', data.loc[0, 'label_atom_id'])
        self.assertEqual('VAL', data.loc[0, 'label_comp_id'])
        self.assertEqual('118', data.loc[0, 'auth_seq_id'])
        self.assertEqual(-7.069, data.loc[0, 'Cartn_x'])
        data = self.residues_aggregation(data, agg_method='first', category='label')
        self.assertEqual(1, data.loc[0, 'id'])
        self.assertEqual('N', data.loc[0, 'type_symbol'])
        self.assertEqual('N', data.loc[0, 'label_atom_id'])
        self.assertEqual('VAL', data.loc[0, 'label_comp_id'])
        self.assertEqual('118', data.loc[0, 'auth_seq_id'])
        self.assertEqual(-7.069, data.loc[0, 'Cartn_x'])

    def test_residues_aggregation_centroid(self):
        data = self.mmcif().read(self.inputcif)
        self.assertEqual(1, data.loc[0, 'id'])
        self.assertEqual('N', data.loc[0, 'type_symbol'])
        self.assertEqual('N', data.loc[0, 'label_atom_id'])
        self.assertEqual('VAL', data.loc[0, 'label_comp_id'])
        self.assertEqual('118', data.loc[0, 'auth_seq_id'])
        self.assertEqual(-7.069, data.loc[0, 'Cartn_x'])
        data = self.residues_aggregation(data, agg_method='centroid', category='label')
        self.assertEqual(1, data.loc[0, 'id'])
        self.assertEqual('N', data.loc[0, 'type_symbol'])
        self.assertEqual('N', data.loc[0, 'label_atom_id'])
        self.assertEqual('VAL', data.loc[0, 'label_comp_id'])
        self.assertEqual('118', data.loc[0, 'auth_seq_id'])
        self.assertAlmostEqual(-7.310, data.loc[0, 'Cartn_x'], places=2)

    def test_reader_pdb_altloc(self):
        data = self.pdb().read(self.inputpdb2, remove_altloc=False, add_atom_altloc=True)
        self.assertIn('label_atom_altloc_id', [k for k in data])
        self.assertIn('auth_atom_altloc_id', [k for k in data])
        self.assertEqual('N.A', data.loc[0, 'label_atom_altloc_id'])
        self.assertEqual('N.B', data.loc[1, 'label_atom_altloc_id'])
        self.assertEqual('CA.A', data.loc[2, 'label_atom_altloc_id'])
        self.assertEqual('CA.B', data.loc[3, 'label_atom_altloc_id'])
        self.assertEqual('C', data.loc[4, 'label_atom_altloc_id'])
        self.assertEqual('O', data.loc[5, 'label_atom_altloc_id'])

    def test_pdb_altloc(self):
        data = self.pdb().read(self.inputpdb2, remove_altloc=False)
        data = self.add_mmcif_atom_altloc(data)
        self.assertIn('label_atom_altloc_id', [k for k in data])
        self.assertIn('auth_atom_altloc_id', [k for k in data])
        self.assertEqual('N.A', data.loc[0, 'label_atom_altloc_id'])
        self.assertEqual('N.B', data.loc[1, 'label_atom_altloc_id'])
        self.assertEqual('CA.A', data.loc[2, 'label_atom_altloc_id'])
        self.assertEqual('CA.B', data.loc[3, 'label_atom_altloc_id'])
        self.assertEqual('C', data.loc[4, 'label_atom_altloc_id'])
        self.assertEqual('O', data.loc[5, 'label_atom_altloc_id'])

    def test_reader_remove_altloc(self):
        data = self.pdb().read(self.inputpdb2, remove_altloc=True)
        self.assertEqual(1, data.loc[0, 'id'])
        self.assertEqual('N', data.loc[0, 'type_symbol'])
        self.assertEqual('N', data.loc[0, 'label_atom_id'])
        self.assertEqual('.', data.loc[0, 'label_alt_id'])

    def test_remove_altloc(self):
        data = self.pdb().read(self.inputpdb2, remove_altloc=False, reset_atom_id=False)
        ndata = self.remove_altloc(data)
        self.assertEqual(1, ndata.loc[0, 'id'])
        self.assertEqual('N', ndata.loc[0, 'type_symbol'])
        self.assertEqual('N', ndata.loc[0, 'label_atom_id'])
        self.assertEqual('.', ndata.loc[0, 'label_alt_id'])

    def test_remove_partial_res(self):
        data = self.pdb().read(self.inputpdb2, remove_altloc=True,
                               remove_partial_res=False, reset_atom_id=True)
        # counting the number of atoms in res with label_seq_id = '25'
        self.assertEqual(16, list(data.label_seq_id.tolist()).count('25'))
        # data = reader.atoms(format_type='pdb', remove_altloc=True,
        #                     remove_partial_res=True, reset_atom_id=True)
        data = self.remove_partial_residues(data)
        # counting the number of atoms in res with label_seq_id = '25'
        self.assertEqual(8, list(data.label_seq_id.tolist()).count('25'))

    def test_fix_label_alt_id(self):
        data = self.pdb().read(self.inputpdb2, pdb_fix_label_alt_id=False)
        data = self.fix_label_alt_id(data)
        self.assertEqual(1, data.loc[0, 'id'])
        self.assertEqual('N', data.loc[0, 'label_atom_id'])
        self.assertEqual('A', data.loc[0, 'label_alt_id'])

    def test_fix_pdb_ins_code(self):
        data = self.pdb().read(self.inputpdb2, pdb_fix_ins_code=False)
        data = self.fix_pdb_ins_code(data)
        self.assertEqual(1, data.loc[0, 'id'])
        self.assertEqual('?', data.loc[0, 'pdbx_PDB_ins_code'])

    def test_fix_type_symbol(self):
        data = self.pdb().read(self.inputpdb2, pdb_fix_type_symbol=False)
        data = self.fix_type_symbol(data)
        self.assertEqual(1, data.loc[0, 'id'])
        self.assertEqual('N', data.loc[0, 'type_symbol'])


if __name__ == '__main__':
    logging.basicConfig(stream=sys.stderr)
    logging.getLogger("proteofav").setLevel(logging.DEBUG)
    suite = unittest.TestLoader().loadTestsFromTestCase(TestStructures)
    unittest.TextTestRunner(verbosity=2).run(suite)
