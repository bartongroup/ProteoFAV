# -*- coding: utf-8 -*-


import os
import sys
import logging
import unittest

from proteofav.sifts import SIFTS

from proteofav.config import defaults as config


class TestSIFTS(unittest.TestCase):
    """Test the DSSP parser methods."""

    def setUp(self):
        """Initialize the framework for testing."""

        root = os.path.abspath(os.path.dirname(__file__))

        self.pdbid = '2pah'
        self.inputsifts = os.path.join(root, "testdata", config.db_sifts,
                                       "{}.xml".format(self.pdbid))
        self.outputsifts = os.path.join(root, "testdata", config.db_tmp,
                                        "{}.xml".format(self.pdbid))
        self.excluded = ()
        self.sifts = SIFTS

        logging.disable(logging.DEBUG)

    def tearDown(self):
        """Remove testing framework."""

        self.pdbid = None
        self.inputsifts = None
        self.outputsifts = None
        self.excluded = None
        self.sifts = None

        logging.disable(logging.NOTSET)

    def test_reader_sifts_data(self):
        data = self.sifts().read(self.inputsifts)
        self.assertEqual(data.loc[0, 'PDB_Annotation'], 'Observed')

    def test_reader_sifts_add_regions(self):
        data = self.sifts().read(self.inputsifts, add_regions=True)
        self.assertEqual(data.loc[0, 'PDB_regionId'], 1)
        self.assertEqual(data.loc[0, 'PDB_regionStart'], 1)
        self.assertEqual(data.loc[0, 'PDB_regionEnd'], 335)

    def test_reader_sifts_add_dbs(self):
        data = self.sifts().read(self.inputsifts, add_dbs=True)
        self.assertIn('PDB_dbVersion', data)
        self.assertEqual(data.loc[0, 'PDB_dbVersion'], '10.17')

    def test_reader_default_excluded(self):
        data = self.sifts().read(self.inputsifts)
        keys = list(data)
        self.assertNotIn("InterPro_dbAccessionId", keys)
        self.assertNotIn("NCBI_dbAccessionId", keys)

    def test_reader_new_excluded(self):
        data = self.sifts().read(self.inputsifts, excluded_cols=self.excluded)
        keys = list(data)
        self.assertIn("InterPro_dbAccessionId", keys)
        self.assertIn("NCBI_dbAccessionId", keys)

    def test_filter_uniprot_id(self):
        data = self.sifts().read(self.inputsifts, uniprot=('P00439',))
        self.assertEqual("P00439", data.loc[0, 'UniProt_dbAccessionId'])

    def test_filter_chain(self):
        data = self.sifts().read(self.inputsifts, chains=('A',))
        self.assertIn("A", data.PDB_entityId.unique())
        self.assertNotIn("B", data.PDB_entityId.unique())

    def test_download_sifts(self):
        self.sifts().download(self.pdbid, filename=self.outputsifts,
                              overwrite=True, decompress=True)
        if os.path.exists(self.outputsifts):
            os.remove(self.outputsifts)


if __name__ == '__main__':
    logging.basicConfig(stream=sys.stderr)
    logging.getLogger("proteofav").setLevel(logging.DEBUG)
    suite = unittest.TestLoader().loadTestsFromTestCase(TestSIFTS)
    unittest.TextTestRunner(verbosity=2).run(suite)
