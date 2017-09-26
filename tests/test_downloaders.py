# -*- coding: utf-8 -*-

import os
import sys
import logging
import unittest

from proteofav.downloaders import Downloader

from proteofav.config import defaults as config


class TestFetchers(unittest.TestCase):
    """Test the dssp parser methods."""

    def setUp(self):
        """Initialize the framework for testing."""

        root = os.path.abspath(os.path.dirname(__file__))

        self.pdbid = "2pah"
        self.outputcif = os.path.join(os.path.join(root, "testdata",
                                                   'tmp', "{}.cif".format(self.pdbid)))
        self.outputsifts = os.path.join(os.path.join(root, "testdata",
                                                     'tmp', "{}.xml".format(self.pdbid)))
        self.Downloader = Downloader

        logging.disable(logging.DEBUG)

    def tearDown(self):
        """Remove testing framework."""

        self.pdbid = None
        self.outputcif = None
        self.outputsifts = None
        self.Downloader = None

        logging.disable(logging.NOTSET)

    def test_download(self):
        url = config.http_pdbe + "entry-files/download/{}_updated.cif".format(self.pdbid)
        self.Downloader(url=url, filename=self.outputcif,
                        decompress=False, overwrite=True)
        os.remove(self.outputcif)

    def test_download_decompress(self):
        url = config.ftp_sifts + "{}.xml.gz".format(self.pdbid)
        self.Downloader(url=url, filename=self.outputsifts,
                        decompress=True, overwrite=True)
        os.remove(self.outputsifts)


if __name__ == '__main__':
    logging.basicConfig(stream=sys.stderr)
    logging.getLogger("proteofav").setLevel(logging.DEBUG)
    suite = unittest.TestLoader().loadTestsFromTestCase(TestFetchers)
    unittest.TextTestRunner(verbosity=2).run(suite)
