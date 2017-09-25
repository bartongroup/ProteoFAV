# -*- coding: utf-8 -*-

import os
import sys
import logging
import unittest
from unittest.mock import patch

from proteofav.downloaders import (Downloader,
                                   _download_structure_from_pdbe,
                                   _download_sifts_from_ebi,
                                   _download_dssp_from_cmbi)

from proteofav.config import defaults as config

root = os.path.abspath(os.path.dirname(__file__))


@patch("proteofav.config.defaults.db_pdb", os.path.join(root, "testdata", 'tmp/'))
@patch("proteofav.config.defaults.db_mmcif", os.path.join(root, "testdata", 'tmp/'))
@patch("proteofav.config.defaults.db_sifts", os.path.join(root, "testdata", 'tmp/'))
@patch("proteofav.config.defaults.db_dssp", os.path.join(root, "testdata", 'tmp/'))
class TestFetchers(unittest.TestCase):
    """Test the dssp parser methods."""

    def setUp(self):
        """Initialize the framework for testing."""

        self.pdbid = "2pah"
        self.Downloader = Downloader
        self.download_structure_from_pdbe = _download_structure_from_pdbe
        self.download_sifts_from_ebi = _download_sifts_from_ebi
        self.download_dssp_from_cmbi = _download_dssp_from_cmbi

        logging.disable(logging.DEBUG)

    def tearDown(self):
        """Remove testing framework."""

        self.pdbid = None
        self.Downloader = None
        self.download_structure_from_pdbe = None
        self.download_sifts_from_ebi = None
        self.download_dssp_from_cmbi = None

        logging.disable(logging.NOTSET)

    def test_downloader_generic(self):
        url_root = config.ftp_dssp
        url_endpoint = "{}.dssp".format(self.pdbid)
        url = url_root + url_endpoint
        filename = "{}.dssp".format(self.pdbid)
        outputfile = os.path.join(config.db_dssp, filename)
        os.makedirs(os.path.join(config.db_dssp), exist_ok=True)
        self.Downloader(url=url, filename=outputfile,
                        decompress=False, override=False)
        os.remove(os.path.join(config.db_dssp, filename))

    def test_downloader_generic_decompress(self):
        url_root = config.ftp_sifts
        url_endpoint = "{}.xml.gz".format(self.pdbid)
        url = url_root + url_endpoint
        filename = "{}.xml.gz".format(self.pdbid)
        outputfile = os.path.join(config.db_sifts, filename)
        os.makedirs(os.path.join(config.db_sifts), exist_ok=True)
        self.Downloader(url=url, filename=outputfile,
                        decompress=True, override=False)
        os.remove(os.path.join(config.db_sifts, filename.rstrip('.gz')))

    def test_download_structure_from_pdbe_pdb(self):
        self.download_structure_from_pdbe(self.pdbid, pdb=True)
        os.remove(os.path.join(config.db_pdb, "{}.pdb".format(self.pdbid)))

    def test_download_structure_from_pdbe_mmcif(self):
        self.download_structure_from_pdbe(self.pdbid, pdb=False)
        os.remove(os.path.join(config.db_mmcif, "{}.cif".format(self.pdbid)))

    def test_download_structure_from_pdbe_mmcif_bio(self):
        self.download_structure_from_pdbe(self.pdbid, pdb=False, bio=True)
        os.remove(os.path.join(config.db_mmcif, "{}_bio.cif".format(self.pdbid)))

    def test_download_dssp_from_cmbi(self):
        self.download_dssp_from_cmbi(self.pdbid)
        os.remove(os.path.join(config.db_dssp, "{}.dssp".format(self.pdbid)))

    def test_download_sifts_from_ebi(self):
        self.download_sifts_from_ebi(self.pdbid)
        os.remove(os.path.join(config.db_sifts, "{}.xml".format(self.pdbid)))


if __name__ == '__main__':
    logging.basicConfig(stream=sys.stderr)
    logging.getLogger("proteofav").setLevel(logging.DEBUG)
    suite = unittest.TestLoader().loadTestsFromTestCase(TestFetchers)
    unittest.TextTestRunner(verbosity=2).run(suite)
