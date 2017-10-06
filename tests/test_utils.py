# -*- coding: utf-8 -*-

import os
import re
import sys
import json
import logging
import requests
import responses
import unittest
import numpy as np
import pandas as pd
import requests_cache

try:
    from mock import patch, MagicMock
except ImportError:
    from unittest.mock import patch, MagicMock

from proteofav.utils import (fetch_from_url_or_retry, row_selector,
                             InputFileHandler, OutputFileHandler,
                             constrain_column_types,
                             exclude_columns, Downloader, GenericInputs)

from proteofav.config import defaults as config


def response_mocker(kwargs, base_url, endpoint_url, status=200,
                    content_type='application/json', post=False, data=None):
    """
    Generates a mocked requests response for a given set of
    kwargs, base url and endpoint url
    """

    url = re.sub('\{\{(?P<m>[a-zA-Z_]+)\}\}', lambda m: "%s" % kwargs.get(m.group(1)),
                 base_url + endpoint_url)
    with responses.RequestsMock() as rsps:
        if post:
            rsps.add(responses.POST, url,
                     body=b'{"data": "some json formatted output"}',
                     status=status, content_type='application/json')
            response = requests.post(url, data=data)

        elif content_type == 'application/json':
            rsps.add(responses.GET, url,
                     body=b'{"data": "some json formatted output"}',
                     status=status, content_type='application/json')
            response = requests.get(url)
        elif content_type == 'text/plain':
            rsps.add(responses.GET, url,
                     body="Some text-based content\n spanning multiple lines",
                     status=status, content_type='text/plain')
            response = requests.get(url)
        else:
            rsps.add(responses.GET, url,
                     body=b"Some other binary stuff...",
                     status=status, content_type='application/octet-stream')
            response = requests.get(url)
    return response


class TestUTILS(unittest.TestCase):
    """Test the utility methods."""

    def setUp(self):
        """Initialize the framework for testing."""

        requests_cache.uninstall_cache()

        self.fetch_from_url_or_retry = fetch_from_url_or_retry
        self.row_selector = row_selector
        self.mock_df = pd.DataFrame(
            [{'label': '1', 'value': 1, 'type': 23.4},
             {'label': '2', 'value': 1, 'type': 1},
             {'label': '3', 'value': 2, 'type': np.nan},
             {'label': '4', 'value': 3, 'type': 123.1},
             {'label': '5', 'value': 5, 'type': 0.32}])

        self.InputFileHandler = InputFileHandler
        self.OutputFileHandler = OutputFileHandler
        self.constrain_column_types = constrain_column_types
        self.exclude_columns = exclude_columns

        self.pdbid = "2pah"
        root = os.path.abspath(os.path.dirname(__file__))
        self.outputcif = os.path.join(os.path.join(root, "testdata",
                                                   "{}.cif".format(self.pdbid)))
        self.outputsifts = os.path.join(os.path.join(root, "testdata",
                                                     "{}.xml".format(self.pdbid)))
        self.Downloader = Downloader
        self.GenericInputs = GenericInputs

    def tearDown(self):
        """Remove testing framework."""

        self.fetch_from_url_or_retry = None
        self.row_selector = None
        self.mock_df = None

        self.InputFileHandler = None
        self.OutputFileHandler = None
        self.constrain_column_types = None
        self.exclude_columns = None
        self.pdbid = None
        self.outputcif = None
        self.outputsifts = None
        self.Downloader = None
        self.GenericInputs = None

    def test_fetch_from_url_or_retry_get_json(self):
        # mocked requests
        identifier = "2pah"
        base_url = config.api_pdbe
        endpoint_url = "pdb/entry/summary/"
        response = response_mocker(kwargs={identifier}, base_url=base_url,
                                   endpoint_url=endpoint_url,
                                   content_type='application/json')
        self.fetch_from_url_or_retry = MagicMock(return_value=response)
        url = base_url + endpoint_url + identifier
        r = self.fetch_from_url_or_retry(url, json=True,
                                         header={'application/json'}).json()
        self.assertEqual(r, json.loads('{"data": "some json formatted output"}'))

    def test_fetch_from_url_or_retry_get_binary(self):
        # mocked requests
        identifier = "P00439"
        base_url = config.api_uniprot
        endpoint_url = "{}.fasta".format(identifier)
        response = response_mocker(kwargs={"P00439.fasta"}, base_url=base_url,
                                   endpoint_url="",
                                   content_type='application/octet-stream')
        self.fetch_from_url_or_retry = MagicMock(return_value=response)
        url = base_url + endpoint_url
        r = self.fetch_from_url_or_retry(url, json=True,
                                         header={'application/octet-stream'},
                                         retry_in=None, wait=0,
                                         n_retries=10, stream=False).content
        self.assertEqual(r, b"Some other binary stuff...")

    def test_fetch_from_url_or_retry_post_json(self):
        # mocked requests
        identifier = "1csb, 2pah"
        base_url = config.api_pdbe
        endpoint_url = "pdb/entry/summary/"
        response = response_mocker(kwargs={}, base_url=base_url,
                                   endpoint_url=endpoint_url,
                                   content_type='application/octet-stream',
                                   post=True, data=identifier)
        self.fetch_from_url_or_retry = MagicMock(return_value=response)
        url = base_url + endpoint_url + identifier
        r = self.fetch_from_url_or_retry(url, json=True, post=True, data=identifier,
                                         header={'application/octet-stream'},
                                         retry_in=None, wait=0,
                                         n_retries=10, stream=False).json()
        self.assertEqual(r, json.loads('{"data": "some json formatted output"}'))

    def test_fetch_from_url_or_retry_get_404(self):
        # mocked requests
        identifier = "P00439"
        base_url = config.api_uniprot
        endpoint_url = "{}.fasta".format(identifier)
        response = response_mocker(kwargs={"P00439.fasta"}, base_url=base_url,
                                   endpoint_url="",
                                   content_type='text/plain', status=404)
        self.fetch_from_url_or_retry = MagicMock(return_value=response)
        url = base_url + endpoint_url
        r = self.fetch_from_url_or_retry(url, json=True, header={'text/plain'},
                                         retry_in=None, wait=0,
                                         n_retries=10, stream=False)
        self.assertEqual(r.status_code, 404)
        self.assertFalse(r.ok)

    def test_fetch_from_url_or_retry_get_500(self):
        # mocked requests
        identifier = "P00439"
        base_url = config.api_uniprot
        endpoint_url = "{}.fasta".format(identifier)
        response = response_mocker(kwargs={"P00439.fasta"}, base_url=base_url,
                                   endpoint_url="",
                                   content_type='text/plain', status=500)
        self.fetch_from_url_or_retry = MagicMock(return_value=response)
        url = base_url + endpoint_url
        r = self.fetch_from_url_or_retry(url, json=True, header={'text/plain'},
                                         retry_in=(500,), wait=1,
                                         n_retries=10, stream=False)
        self.assertEqual(r.status_code, 500)
        self.assertFalse(r.ok)

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
        config.db_mmcif = "mmcif"
        valid = os.path.join(os.path.dirname(__file__), "testdata",
                             config.db_mmcif, "2pah.cif")
        self.InputFileHandler(valid)

        invalid = os.path.join(os.path.dirname(__file__), "testdata",
                               config.db_mmcif, "null.cif")
        with self.assertRaises(IOError) or self.assertRaises(OSError):
            self.InputFileHandler(invalid)

    def test_outputfilehandler(self):
        config.db_mmcif = "mmcif"
        valid = os.path.join(os.path.dirname(__file__), "testdata",
                             config.db_mmcif, "2pah.cif")
        self.OutputFileHandler(valid, overwrite=True)

        with self.assertRaises(OSError):
            self.OutputFileHandler(valid, overwrite=False)

        invalid = os.path.join(os.path.dirname(__file__), "testdata",
                               "NEW_DIR", "null.cif")
        with self.assertRaises(OSError):
            self.InputFileHandler(invalid)

    def test_download(self):
        url = config.cif_fetch + "download/{}_updated.cif".format(self.pdbid)
        self.Downloader(url=url, filename=self.outputcif,
                        decompress=False, overwrite=True)
        os.remove(self.outputcif)

    def test_download_decompress(self):
        url = config.sifts_fetch + "{}.xml.gz".format(self.pdbid)
        self.Downloader(url=url, filename=self.outputsifts,
                        decompress=True, overwrite=True)
        os.remove(self.outputsifts)

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

    def test_generic_inputs(self):
        # identifier
        t = self.GenericInputs(identifier="test1")
        v = t._get_identifier()
        self.assertEqual(v, "test1")
        v = t._get_identifier(identifier="test2")
        self.assertEqual(v, "test2")
        t = self.GenericInputs()
        v = t._get_identifier(identifier="test2")
        self.assertEqual(v, "test2")
        v = self.GenericInputs(identifier="test1")._get_identifier()
        self.assertEqual(v, "test1")
        v = self.GenericInputs()._get_identifier(identifier="test2")
        self.assertEqual(v, "test2")
        v = self.GenericInputs(identifier="test1")._get_identifier(identifier="test2")
        self.assertEqual(v, "test2")

        # filename
        v = self.GenericInputs(filename="test1")._get_filename()
        self.assertEqual(v, "test1")
        v = self.GenericInputs()._get_filename(filename="test2")
        self.assertEqual(v, "test2")

        # table
        t = self.GenericInputs(table="test1")
        v = t._get_table(table=None)
        self.assertEqual(v, "test1")
        v = t._get_table(table="test2")
        self.assertEqual(v, "test2")


if __name__ == '__main__':
    logging.basicConfig(stream=sys.stderr)
    logging.getLogger("proteofav.config").setLevel(logging.CRITICAL)
    suite = unittest.TestLoader().loadTestsFromTestCase(TestUTILS)
    unittest.TextTestRunner(verbosity=2).run(suite)
