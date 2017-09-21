# -*- coding: utf-8 -*-

import sys
import logging
import unittest
import pandas as pd
import requests_cache

try:
    from mock import patch, MagicMock
except ImportError:
    from unittest.mock import patch, MagicMock

from proteofav.utils import row_selector


class TestUTILS(unittest.TestCase):
    """Test the utility methods."""

    def setUp(self):
        """Initialize the framework for testing."""

        requests_cache.uninstall_cache()
        self.row_selector = row_selector

        logging.disable(logging.DEBUG)

    def tearDown(self):
        """Remove testing framework."""

        self.row_selector = None

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


if __name__ == '__main__':
    logging.basicConfig(stream=sys.stderr)
    logging.getLogger("proteofav").setLevel(logging.DEBUG)
    suite = unittest.TestLoader().loadTestsFromTestCase(TestUTILS)
    unittest.TextTestRunner(verbosity=2).run(suite)
