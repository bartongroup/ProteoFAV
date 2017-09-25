# -*- coding: utf-8 -*-

import os
import sys
import logging
import unittest

from proteofav.parsers import (parse_mmcif_atoms_from_file,
                               parse_pdb_atoms_from_file,
                               parse_dssp_from_file,
                               parse_sifts_residues_from_file,
                               _parse_sifts_regions_from_file,
                               _parse_sifts_dbs_from_file)

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
        self.inputpdb = os.path.join(root, "testdata", config.db_pdb,
                                     "{}.pdb".format(self.pdbid))
        self.inputpdb2 = os.path.join(root, "testdata", config.db_pdb,
                                      "{}.pdb".format(self.pdbid2))
        self.inputdssp = os.path.join(root, "testdata", config.db_dssp,
                                      "{}.dssp".format(self.pdbid))
        self.inputbiodssp = os.path.join(root, "testdata", config.db_dssp,
                                         "{}_bio.dssp".format(self.pdbid))
        self.inputsifts = os.path.join(root, "testdata", config.db_sifts,
                                       "{}.xml".format(self.pdbid))

        self.parser_mmcif = parse_mmcif_atoms_from_file
        self.parser_pdb = parse_pdb_atoms_from_file
        self.parser_dssp = parse_dssp_from_file
        self.parser_sifts = parse_sifts_residues_from_file
        self.parser_sifts_regions = _parse_sifts_regions_from_file
        self.parser_sifts_dbs = _parse_sifts_dbs_from_file

        logging.disable(logging.DEBUG)

    def tearDown(self):
        """Remove testing framework."""

        self.pdbid = None
        self.pdbid2 = None
        self.inputcif = None
        self.inputbiocif = None
        self.inputpdb = None
        self.inputpdb2 = None
        self.inputdssp = None
        self.inputsifts = None

        self.parser_mmcif = None
        self.parser_pdb = None
        self.parser_dssp = None
        self.parser_sifts = None
        self.parser_sifts_regions = None
        self.parser_sifts_dbs = None

        logging.disable(logging.NOTSET)

    def test_parser_cif(self):
        data = self.parser_mmcif(self.inputcif)
        self.assertEqual(data.loc[0, 'label_asym_id'], 'A')
        self.assertEqual(data.loc[2686, 'label_asym_id'], 'B')
        self.assertEqual(data.loc[5315, 'label_asym_id'], 'C')
        self.assertEqual(data.loc[5316, 'label_asym_id'], 'D')
        self.assertEqual(data.loc[1, 'label_atom_id'], 'CA')
        self.assertEqual(data.loc[2687, 'label_atom_id'], 'CA')
        self.assertEqual(data.loc[5315, 'label_atom_id'], 'FE')
        self.assertEqual(data.loc[5316, 'label_atom_id'], 'FE')

    def test_parser_mmcif_bio(self):
        data = self.parser_mmcif(self.inputbiocif)
        self.assertEqual(data.loc[0, 'label_asym_id'], 'A')
        self.assertEqual(data.loc[5372, 'label_asym_id'], 'B')
        self.assertEqual(data.loc[10630, 'label_asym_id'], 'C')
        self.assertEqual(data.loc[10632, 'label_asym_id'], 'D')
        self.assertEqual(data.loc[2686, 'label_asym_id'], 'AA')
        self.assertEqual(data.loc[8001, 'label_asym_id'], 'BA')
        self.assertEqual(data.loc[10631, 'label_asym_id'], 'CA')
        self.assertEqual(data.loc[10633, 'label_asym_id'], 'DA')

    def test_parser_pdb(self):
        data = self.parser_pdb(self.inputpdb)
        self.assertEqual(data.loc[0, 'label_asym_id'], 'A')
        self.assertEqual(data.loc[2686, 'label_asym_id'], 'B')
        self.assertEqual(data.loc[5315, 'label_asym_id'], 'A')
        self.assertEqual(data.loc[5316, 'label_asym_id'], 'B')
        self.assertEqual(data.loc[1, 'label_atom_id'], 'CA')

    def test_parser_dssp(self):
        data = self.parser_dssp(self.inputdssp)
        self.assertEqual(data.loc[0, 'CHAIN'], 'A')
        self.assertEqual(data.loc[0, 'RES'], '118')
        self.assertEqual(data.loc[331, 'CHAIN'], 'B')
        self.assertEqual(data.loc[331, 'RES'], '118')

    def test_parser_dssp_bio(self):
        data = self.parser_dssp(self.inputbiodssp)
        self.assertEqual(data.loc[0, 'CHAIN'], 'A')
        self.assertEqual(data.loc[0, 'RES'], '118')
        self.assertEqual(data.loc[662, 'CHAIN'], 'B')
        self.assertEqual(data.loc[662, 'RES'], '118')

    def test_parser_sifts(self):
        data = self.parser_sifts(self.inputsifts)
        self.assertEqual(data.loc[0, 'PDB_Annotation'], 'Observed')

    def test_parser_sifts_regions(self):
        data = self.parser_sifts_regions(self.inputsifts)
        self.assertEqual(data['B']['PDB']['1']['dbAccessionId'], self.pdbid)
        self.assertEqual(data['B']['PDB']['1']['start'], 1)
        self.assertEqual(data['B']['PDB']['1']['end'], 335)
        self.assertEqual(data['B']['PDB']['1']['dbCoordSys'], 'PDBresnum')

    def test_parser_sifts_dbs(self):
        data = self.parser_sifts_dbs(self.inputsifts)
        self.assertEqual(data['UniProt']['dbSource'], 'UniProt')
        self.assertEqual(data['UniProt']['dbCoordSys'], 'UniProt')
        self.assertEqual(data['UniProt']['dbVersion'], '2017.03')


if __name__ == '__main__':
    logging.basicConfig(stream=sys.stderr)
    logging.getLogger("proteofav").setLevel(logging.DEBUG)
    suite = unittest.TestLoader().loadTestsFromTestCase(TestStructures)
    unittest.TextTestRunner(verbosity=2).run(suite)
