# -*- coding: utf-8 -*-

import os
import sys
import logging
import unittest

try:
    from mock import patch
except ImportError:
    # python 3.5
    from unittest.mock import patch

from proteofav.sifts import (parse_sifts_residues, sifts_best,
                             _parse_sifts_regions_from_file,
                             _parse_sifts_dbs_from_file, select_sifts,
                             filter_sifts, download_sifts, SIFTS)

from proteofav.config import defaults

root = os.path.abspath(os.path.dirname(__file__))
defaults.db_sifts = os.path.join(root, "testdata", "sifts")


@patch("proteofav.sifts.defaults", defaults)
class TestSIFTSParser(unittest.TestCase):
    """Test the SIFTS parser methods."""

    def setUp(self):
        """Initialize the framework for testing."""
        self.example_xml = os.path.join(os.path.dirname(__file__),
                                        "testdata", "sifts", "2pah.xml")
        self.residues_parser = parse_sifts_residues
        self.output_sifts = os.path.join(os.path.dirname(__file__),
                                         "testdata", "2pah.xml")
        self.pdbid = '2pah'
        self.uniprot_id = 'P00439'
        self.pdb_best = sifts_best
        self.parser_sifts_regions = _parse_sifts_regions_from_file
        self.parser_sifts_dbs = _parse_sifts_dbs_from_file
        self.filter_sifts = filter_sifts
        self.download_sifts = download_sifts
        self.SIFTS = SIFTS

    def tearDown(self):
        """Remove testing framework."""

        self.example_xml = None
        self.residue_parser = None
        self.output_sifts = None
        self.pdbid = None
        self.uniprot_id = None
        self.pdb_best = None
        self.parser_sifts_regions = None
        self.parser_sifts_dbs = None
        self.filter_sifts = None
        self.download_sifts = None
        self.SIFTS = None

    def test_to_table_sifts_residues_and_regions(self):
        """
        Tests the parsing real SIFTS xml files.
        This test focuses on the method that parses the residue entries
        and the region entries.

        Some checks are made to whether the parsed keys and values
        are the ones we are expecting.
        """

        data = self.residues_parser(self.example_xml)

        # number of values per column (or rows)
        self.assertEqual(len(data), 670)

        # number of keys (or columns)
        self.assertEqual(len(data.columns.values), 34)

        # state of the residue
        self.assertEqual(data.loc[0, 'PDB_Annotation'], 'Observed')

        # check whether there are particular keys
        self.assertIn('CATH_dbAccessionId', data.columns.values)

        # check the values for particular entries
        self.assertTrue(data['CATH_dbAccessionId'][0] == '1.10.800.10')
        self.assertEqual(data['PDB_dbChainId'][0], 'A')
        self.assertEqual(data['PDB_dbResName'][0], 'VAL')

        # check whether there are particular keys
        self.assertIn('UniProt_regionId', data.columns.values)
        self.assertIn('CATH_regionId', data.columns.values)

        # check the values of particular entries
        self.assertTrue(data['CATH_regionId'][0] == 1)

    def test_to_table_uniprot_pdb_sifts_best(self):
        """
        Testing the PDBe API for mapping between UniProt and PDB
        accession identifiers.
        """

        data = self.pdb_best(self.uniprot_id, first=True)

        # number of values
        self.assertEqual(len(data), 10)

        # number of keys (or columns)
        self.assertIn('coverage', data)
        self.assertIn('pdb_id', data)
        self.assertIn('chain_id', data)

    def test_parser_sifts_regions(self):
        data = self.parser_sifts_regions(self.example_xml)
        self.assertEqual(data['B']['PDB']['1']['dbAccessionId'], self.pdbid)
        self.assertEqual(data['B']['PDB']['1']['start'], 1)
        self.assertEqual(data['B']['PDB']['1']['end'], 335)
        self.assertEqual(data['B']['PDB']['1']['dbCoordSys'], 'PDBresnum')

    def test_parser_sifts_dbs(self):
        data = self.parser_sifts_dbs(self.example_xml)
        self.assertEqual(data['UniProt']['dbSource'], 'UniProt')
        self.assertEqual(data['UniProt']['dbCoordSys'], 'UniProt')
        self.assertEqual(data['UniProt']['dbVersion'], '2017.03')

    def test_sifts_filter_default_excluded(self):
        data = self.residues_parser(self.example_xml, excluded_cols=())
        keys = [k for k in data]
        self.assertIn("InterPro_dbAccessionId", keys)
        self.assertIn("NCBI_dbAccessionId", keys)

        exc_cols = ("InterPro", "GO", "EC", "NCBI")
        data = self.residues_parser(self.example_xml, excluded_cols=exc_cols)
        keys = [k for k in data]
        self.assertNotIn("InterPro_dbAccessionId", keys)
        self.assertNotIn("NCBI_dbAccessionId", keys)

    def test_sifts_filter_uniprot_id(self):
        data = self.residues_parser(self.example_xml)
        data = self.filter_sifts(data, uniprot=('P00439',))
        self.assertEqual("P00439", data.loc[0, 'UniProt_dbAccessionId'])

    def test_sifts_filter_chain(self):
        data = self.residues_parser(self.example_xml)
        data = self.filter_sifts(data, chains=('A',))
        self.assertIn("A", data.PDB_entityId.unique())
        self.assertNotIn("B", data.PDB_entityId.unique())

    def test_sifts_filter_chain_auth(self):
        data = self.residues_parser(self.example_xml)
        data = self.filter_sifts(data, chain_auth=('A',))
        self.assertIn("A", data.PDB_dbChainId.unique())
        self.assertNotIn("B", data.PDB_dbChainId.unique())

    def test_sifts_filter_res(self):
        data = self.residues_parser(self.example_xml)
        data = self.filter_sifts(data, res=('285',))
        self.assertIn("285", data.PDB_dbResNum.unique())
        self.assertNotIn("286", data.PDB_dbResNum.unique())

    def test_sifts_filter_site(self):
        data = self.residues_parser(self.example_xml)
        data = self.filter_sifts(data, site=('118',))
        self.assertIn("118", data.UniProt_dbResNum.unique())
        self.assertNotIn("119", data.UniProt_dbResNum.unique())

    def test_download_sifts(self):
        self.download_sifts(self.pdbid, filename=self.output_sifts,
                            overwrite=True)
        if os.path.exists(self.output_sifts):
            os.remove(self.output_sifts)

    def test_main_SIFTS(self):
        # read
        data = self.SIFTS.read(self.example_xml)
        self.assertTrue(data['CATH_dbAccessionId'][0] == '1.10.800.10')
        self.assertEqual(data['PDB_dbChainId'][0], 'A')
        self.assertEqual(data['PDB_dbResName'][0], 'VAL')
        self.assertIn('UniProt_regionId', data.columns.values)
        self.assertIn('CATH_regionId', data.columns.values)
        # download
        self.SIFTS.download(self.pdbid, filename=self.output_sifts,
                            overwrite=True)
        if os.path.exists(self.output_sifts):
            os.remove(self.output_sifts)
        # select
        self.SIFTS.select(self.pdbid)


if __name__ == '__main__':
    logging.basicConfig(stream=sys.stderr)
    logging.getLogger("proteofav.config").setLevel(logging.CRITICAL)
    suite = unittest.TestLoader().loadTestsFromTestCase(TestSIFTSParser)
    unittest.TextTestRunner(verbosity=2).run(suite)
