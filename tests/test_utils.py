# -*- coding: utf-8 -*-

import os
import sys
import logging
import unittest
import numpy as np
import pandas as pd
import requests_cache

from proteofav.utils import (row_selector,
                             InputFileHandler,
                             OutputFileHandler,
                             constrain_column_types,
                             exclude_columns,
                             splitting_up_by_key,
                             merging_down_by_key,
                             flatten_nested_structure,
                             refactor_key_val_singletons)

from proteofav.config import defaults as config

log = logging.getLogger(__name__)


class TestUTILS(unittest.TestCase):
    """Test the utility methods."""

    def setUp(self):
        """Initialize the framework for testing."""

        requests_cache.uninstall_cache()
        self.row_selector = row_selector
        self.InputFileHandler = InputFileHandler
        self.OutputFileHandler = OutputFileHandler
        self.constrain_column_types = constrain_column_types
        self.exclude_columns = exclude_columns
        self.mock_df = pd.DataFrame(
            [{'label': '1', 'value': 1, 'type': 23.4},
             {'label': '2', 'value': 1, 'type': 1},
             {'label': '3', 'value': 2, 'type': np.nan},
             {'label': '4', 'value': 3, 'type': 123.1},
             {'label': '5', 'value': 5, 'type': 0.32}])

        self.merging_down_by_key = merging_down_by_key
        self.splitting_up_by_key = splitting_up_by_key
        self.flatten_nested_structure = flatten_nested_structure
        self.refactor_key_val_singletons = refactor_key_val_singletons
        self.vars_mock = pd.DataFrame([{'xrefs_id': 'id1', 'other_key': []},
                                       {'xrefs_id': 'id2', 'other_key': 123},
                                       {'xrefs_id': 'id1', 'other_key': []},
                                       {'xrefs_id': 'id3', 'other_key': 'string'},
                                       {'xrefs_id': 'id2', 'other_key': 245},
                                       {'xrefs_id': 'id3', 'other_key': np.nan}])
        self.vars_mock2 = pd.DataFrame([{'xrefs_id': ['id1', 'id2'], 'other_key': 123},
                                        {'xrefs_id': ['id1', 'id2', 'id3'], 'other_key': 456}])
        self.json_mock = {
            'M1': [{'d1': 1},
                   {'d1': 2},
                   {'d2': 3},
                   {'d3': {'dd3': 'dd3'}}],
            'M2': 'value',
            'M3': {'x1':
                       {'x2': 1, 'x3': 2}
                   },
            'M4': 'four',
            'M5': [1, 2, 3],
            'M6': {'z1': 'z1'}
        }
        logging.disable(logging.DEBUG)

    def tearDown(self):
        """Remove testing framework."""

        self.row_selector = None
        self.InputFileHandler = None
        self.OutputFileHandler = None
        self.constrain_column_types = None
        self.exclude_columns = None
        self.mock_df = None
        self.merging_down_by_key = None
        self.splitting_up_by_key = None
        self.flatten_nested_structure = None
        self.refactor_key_val_singletons = None
        self.vars_mock = None
        self.vars_mock2 = None
        self.json_mock = None

        logging.disable(logging.NOTSET)

    def test_row_selector(self):
        d = self.row_selector(self.mock_df, key='value', value=3)
        self.assertEqual(len(d.index), 1)
        d = self.row_selector(self.mock_df, key='value', value=3, reverse=True)
        self.assertEqual(len(d.index), 4)
        d = self.row_selector(self.mock_df, key='value', value='first')
        self.assertEqual(len(d.index), 2)
        d = self.row_selector(self.mock_df, key='value', value=(2, 3))
        self.assertEqual(len(d.index), 2)

    def test_inputfilehandler(self):
        valid = os.path.join(os.path.dirname(__file__), "testdata",
                             config.db_mmcif, "2pah.cif")
        self.InputFileHandler(valid)

        invalid = os.path.join(os.path.dirname(__file__), "testdata",
                               config.db_mmcif, "null.cif")
        with self.assertRaises(IOError) or self.assertRaises(OSError):
            self.InputFileHandler(invalid)

    def test_outputfilehandler(self):
        valid = os.path.join(os.path.dirname(__file__), "testdata",
                             config.db_mmcif, "2pah.cif")
        self.OutputFileHandler(valid, overwrite=True)

        with self.assertRaises(OSError):
            self.OutputFileHandler(valid, overwrite=False)

        invalid = os.path.join(os.path.dirname(__file__), "testdata",
                               config.db_tmp, "NEW_DIR", "null.cif")
        with self.assertRaises(OSError):
            self.InputFileHandler(invalid)

    def test_constrain_column_types(self):
        dtypes = {'type': 'float64',
                  'value': 'int64',
                  'label': 'object'}
        dnans = {'type': 0.0}
        dreplaces = {'label': ['1', 'i']}

        self.mock_df = self.constrain_column_types(self.mock_df, dtypes)
        self.assertEqual(self.mock_df["type"].dtype, np.float64)
        self.assertEqual(self.mock_df["value"].dtype, np.int64)
        self.assertEqual(self.mock_df["label"].dtype, np.object)

        self.mock_df = self.constrain_column_types(self.mock_df, dtypes,
                                                   nan_value_dict=dnans)
        self.assertEqual(self.mock_df["type"].dtype, np.float64)
        self.assertEqual(self.mock_df.loc[2, "type"], 0.0)

        self.mock_df = self.constrain_column_types(self.mock_df, dtypes, dnans,
                                                   replace_value_dict=dreplaces)
        self.assertEqual(self.mock_df["label"].dtype, np.object)
        self.assertEqual(self.mock_df.loc[0, "label"], "i")

    def test_exclude_columns(self):
        self.assertEqual(len(self.mock_df.columns), 3)
        self.mock_df = self.exclude_columns(self.mock_df, excluded=("type",))
        self.assertEqual(len(self.mock_df.columns), 2)
        self.assertNotIn("type", self.mock_df)

    def test_merging_down_by_key(self):
        table = self.merging_down_by_key(self.vars_mock, key='xrefs_id')
        self.assertEqual(len(self.vars_mock), 6)
        self.assertEqual(len(table), 3)
        self.assertEqual(table.loc[0, 'xrefs_id'], 'id1')
        self.assertEqual(table.loc[1, 'other_key'], (123, 245))
        self.assertEqual(table.loc[2, 'other_key'], 'string')

    def test_splitting_up_by_key(self):
        table = self.splitting_up_by_key(self.vars_mock2, key='xrefs_id')
        self.assertEqual(len(self.vars_mock2), 2)
        self.assertEqual(len(table), 5)
        self.assertEqual(table.loc[0, 'xrefs_id'], 'id1')
        self.assertEqual(table.loc[1, 'xrefs_id'], 'id2')
        self.assertEqual(table.loc[0, 'other_key'], 123)
        self.assertEqual(table.loc[1, 'other_key'], 123)
        self.assertEqual(table.loc[2, 'xrefs_id'], 'id1')
        self.assertEqual(table.loc[3, 'xrefs_id'], 'id2')
        self.assertEqual(table.loc[4, 'xrefs_id'], 'id3')
        self.assertEqual(table.loc[2, 'other_key'], 456)
        self.assertEqual(table.loc[3, 'other_key'], 456)
        self.assertEqual(table.loc[4, 'other_key'], 456)

    def test_flatten_nested_structure(self):
        data = {}
        self.flatten_nested_structure(self.json_mock, data)
        # keys
        self.assertIn('M1_d1', data)
        self.assertIn('M1_d2', data)
        self.assertIn('M1_d3_dd3', data)
        self.assertNotIn('M1', data)
        self.assertIn('M2', data)
        self.assertIn('M3_x1_x2', data)
        self.assertIn('M3_x1_x3', data)
        self.assertNotIn('M3', data)
        self.assertIn('M4', data)
        self.assertIn('M5', data)
        self.assertIn('M6_z1', data)
        # values
        self.assertEqual(data['M1_d1'], [1, 2])
        self.assertEqual(data['M1_d2'], [3])
        self.assertEqual(data['M1_d3_dd3'], ['dd3'])
        self.assertEqual(data['M2'], ['value'])
        self.assertEqual(data['M3_x1_x2'], [1])
        self.assertEqual(data['M3_x1_x3'], [2])
        self.assertEqual(data['M4'], ['four'])
        self.assertEqual(data['M5'], [[1, 2, 3]])
        self.assertEqual(data['M6_z1'], ['z1'])

    def test_refactor_key_val_singletons(self):
        data = {}
        self.flatten_nested_structure(self.json_mock, data)
        data = self.refactor_key_val_singletons(data)
        # keys
        self.assertIn('M1_d1', data)
        self.assertIn('M1_d2', data)
        self.assertIn('M1_d3_dd3', data)
        self.assertNotIn('M1', data)
        self.assertIn('M2', data)
        self.assertIn('M3_x1_x2', data)
        self.assertIn('M3_x1_x3', data)
        self.assertNotIn('M3', data)
        self.assertIn('M4', data)
        self.assertIn('M5', data)
        self.assertIn('M6_z1', data)
        # values
        self.assertEqual(data['M1_d1'], [1, 2])
        self.assertEqual(data['M1_d2'], 3)
        self.assertEqual(data['M1_d3_dd3'], 'dd3')
        self.assertEqual(data['M2'], 'value')
        self.assertEqual(data['M3_x1_x2'], 1)
        self.assertEqual(data['M3_x1_x3'], 2)
        self.assertEqual(data['M4'], 'four')
        self.assertEqual(data['M5'], [1, 2, 3])
        self.assertEqual(data['M6_z1'], 'z1')


if __name__ == '__main__':
    logging.basicConfig(stream=sys.stderr)
    logging.getLogger("proteofav").setLevel(logging.DEBUG)
    suite = unittest.TestLoader().loadTestsFromTestCase(TestUTILS)
    unittest.TextTestRunner(verbosity=2).run(suite)
