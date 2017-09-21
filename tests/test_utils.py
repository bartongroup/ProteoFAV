# -*- coding: utf-8 -*-

import os
import sys
import logging
import unittest
import pandas as pd
import requests_cache

try:
    from mock import patch, MagicMock
except ImportError:
    from unittest.mock import patch, MagicMock

from proteofav.config import Defaults
from proteofav.utils import (row_selector, InputFileHandler, OutputFileHandler)

log = logging.getLogger(__name__)

defaults = Defaults("config.txt")


@patch("proteofav.structures.defaults", defaults)
class TestUTILS(unittest.TestCase):
    """Test the utility methods."""

    def setUp(self):
        """Initialize the framework for testing."""

        requests_cache.uninstall_cache()
        self.row_selector = row_selector
        self.InputFileHandler = InputFileHandler
        self.OutputFileHandler = OutputFileHandler

        logging.disable(logging.DEBUG)

    def tearDown(self):
        """Remove testing framework."""

        self.row_selector = None
        self.InputFileHandler = None
        self.OutputFileHandler = None

        logging.disable(logging.NOTSET)

    def test_row_selector(self):
        data = pd.DataFrame([{'label': '1', 'value': 1},
                             {'label': '2', 'value': 1},
                             {'label': '3', 'value': 2},
                             {'label': '4', 'value': 3},
                             {'label': '5', 'value': 5}])
        d = self.row_selector(data, key='value', value=3, method='equals')
        self.assertEqual(len(d.index), 1)
        d = self.row_selector(data, key='value', value=3, method='diffs')
        self.assertEqual(len(d.index), 4)
        d = self.row_selector(data, key='value', value=None, method='first')
        self.assertEqual(len(d.index), 2)
        d = self.row_selector(data, key='value', value=(2, 3), method='isin')
        self.assertEqual(len(d.index), 2)

    def test_inputfilehandler(self):
        valid = os.path.join(os.path.dirname(__file__), "CIF/2pah.cif")
        self.InputFileHandler(valid)

        invalid = os.path.join(os.path.dirname(__file__), "CIF/null.cif")
        with self.assertRaises(IOError) or self.assertRaises(OSError):
            self.InputFileHandler(invalid)

    def test_outputfilehandler(self):
        valid = os.path.join(os.path.dirname(__file__), "CIF/2pah.cif")
        self.OutputFileHandler(valid, overwrite=True)

        with self.assertRaises(OSError):
            self.OutputFileHandler(valid, overwrite=False)

        invalid = os.path.join(os.path.dirname(__file__), "CIF/NEW_DIR/null.cif")
        with self.assertRaises(OSError):
            self.InputFileHandler(invalid)


if __name__ == '__main__':
    logging.basicConfig(stream=sys.stderr)
    logging.getLogger("proteofav").setLevel(logging.DEBUG)
    suite = unittest.TestLoader().loadTestsFromTestCase(TestUTILS)
    unittest.TextTestRunner(verbosity=2).run(suite)
