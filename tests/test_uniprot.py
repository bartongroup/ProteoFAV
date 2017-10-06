# -*- coding: utf-8 -*-

import logging
import unittest
import numpy as np

try:
    import mock
except ImportError:
    import unittest.mock as mock
from requests.exceptions import HTTPError

from proteofav.variants import _uniprot_ensembl_mapping

from proteofav.uniprot import (fetch_uniprot_sequence, fetch_uniprot_formal_specie, _uniprot_info,
                               _fetch_uniprot_gff, map_gff_features_to_sequence,
                               _uniprot_to_ensembl_xref)
from proteofav.utils import fetch_from_url_or_retry

logging.getLogger('proteofav').setLevel(logging.CRITICAL)  # turn off logging


class TestUNIPROTParser(unittest.TestCase):
    """Test UniProt fetcher/parser methods."""

    def setUp(self):
        """Initialize the framework for testing."""
        self.uniprot_id = 'O96013'
        self.uniprot_id_error1 = ''
        self.uniprot_id_error2 = '1234 ads'
        self.uniprot_id_error3 = 1234
        self.uniprot_id_error4 = ()
        self.uniprot_id_error5 = []
        self.uniprot_info = _uniprot_info
        self.uniprot_ensembl = _uniprot_ensembl_mapping

        self.fetch_from_url_or_retry = fetch_from_url_or_retry
        self.get_uniprot_sequence = fetch_uniprot_sequence
        self.get_uniprot_organism = fetch_uniprot_formal_specie
        self.fetch_uniprot_gff = _fetch_uniprot_gff
        self.uniprot_info = _uniprot_info
        self.map_gff_features_to_sequence = map_gff_features_to_sequence
        self.uniprot_to_ensembl_xref = _uniprot_to_ensembl_xref
        self.mock_url = 'www.mockurl.com'
        self.ccc2_sequence = ("MREVILAVHGMTCSACTNTINTQLRALKGVTKCDISLVTNECQVTYDNEVTADSIKEIIE"
                              "DCGFDCEILRDSEITAISTKEGLLSVQGMTCGSCVSTVTKQVEGIEGVESVVVSLVTEEC"
                              "HVIYEPSKTTLETAREMIEDCGFDSNIIMDGNGNADMTEKTVILKVTKAFEDESPLILSS"
                              "VSERFQFLLDLGVKSIEISDDMHTLTIKYCCNELGIRDLLRHLERTGYKFTVFSNLDNTT"
                              "QLRLLSKEDEIRFWKKNSIKSTLLAIICMLLYMIVPMMWPTIVQDRIFPYKETSFVRGLF"
                              "YRDILGVILASYIQFSVGFYFYKAAWASLKHGSGTMDTLVCVSTTCAYTFSVFSLVHNMF"
                              "HPSSTGKLPRIVFDTSIMIISYISIGKYLETLAKSQTSTALSKLIQLTPSVCSIISDVER"
                              "NETKEIPIELLQVNDIVEIKPGMKIPADGIITRGESEIDESLMTGESILVPKKTGFPVIA"
                              "GSVNGPGHFYFRTTTVGEETKLANIIKVMKEAQLSKAPIQGYADYLASIFVPGILILAVL"
                              "TFFIWCFILNISANPPVAFTANTKADNFFICLQTATSVVIVACPCALGLATPTAIMVGTG"
                              "VGAQNGVLIKGGEVLEKFNSITTFVFDKTGTLTTGFMVVKKFLKDSNWVGNVDEDEVLAC"
                              "IKATESISDHPVSKAIIRYCDGLNCNKALNAVVLESEYVLGKGIVSKCQVNGNTYDICIG"
                              "NEALILEDALKKSGFINSNVDQGNTVSYVSVNGHVFGLFEINDEVKHDSYATVQYLQRNG"
                              "YETYMITGDNNSAAKRVAREVGISFENVYSDVSPTGKCDLVKKIQDKEGNNKVAVVGDGI"
                              "NDAPALALSDLGIAISTGTEIAIEAADIVILCGNDLNTNSLRGLANAIDISLKTFKRIKL"
                              "NLFWALCYNIFMIPIAMGVLIPWGITLPPMLAGLAMAFSSVSVVLSSLMLKKWTPPDIES"
                              "HGISDFKSKFSIGNFWSRLFSTRAIAGEQDIESQAGLMSNEEVL")

    def tearDown(self):
        """Remove testing framework."""

        self.uniprot_id = None
        self.uniprot_id_error1 = None
        self.uniprot_id_error2 = None
        self.uniprot_id_error3 = None
        self.uniprot_id_error4 = None
        self.uniprot_id_error5 = None
        self.uniprot_info = None
        self.uniprot_ensembl = None

        self.fetch_from_url_or_retry = None
        self.get_uniprot_sequence = None
        self.get_uniprot_organism = None
        self.fetch_uniprot_gff = None
        self.uniprot_info = None
        self.map_gff_features_to_sequence = None
        self.uniprot_to_ensembl_xref = None
        self.mock_url = None
        self.ccc2_sequence = None

    def test_to_table_uniprot_info(self):
        """
        Tests the fetching and parsing real UniProt ids.
        This test focuses on the method that parses the residue entries.

        Some checks are made to whether the parsed keys and values
        are the ones we are expecting.
        """

        data = self.uniprot_info(self.uniprot_id)

        # number of values per column (or rows)
        self.assertEqual(len(data), 1)

        # number of keys (or columns)
        self.assertEqual(len(data.columns.values), 8)

        # check whether there are particular keys
        self.assertIn('Sequence', data.columns.values)

        # check the values for particular entries
        self.assertTrue(data['Length'][0] == 591)
        self.assertTrue(data['Status'][0] == "reviewed")
        self.assertEqual(data['Entry name'][0], 'PAK4_HUMAN')

    def test_to_table_uniprot_ensembl_mapping(self):
        """
        Tests the fetching and parsing real UniProt ids.
        This test focuses on the method that parses the residue entries.

        Some checks are made to whether the parsed keys and values
        are the ones we are expecting.
        """

        data = self.uniprot_ensembl(self.uniprot_id, 'homo_sapiens')

        # number of values per column (or rows)
        self.assertEqual(len(data), 1)

        # number of keys (or columns)
        self.assertEqual(len(data.columns.values), 3)

        # check whether there are particular keys
        self.assertIn('TRANSCRIPT', data.columns.values)

        # check the values for particular entries
        self.assertEqual(data['GENE'][0], 'ENSG00000130669')

    def test_obsolete_uniprot_accession(self):
        data = self.uniprot_info('Q91887')
        self.assertEqual(len(data), 1)

        # number of keys (or columns)
        self.assertEqual(len(data.columns.values), 8)

        # check the values for particular entries
        self.assertTrue(np.isnan(data.iloc[0, -1]))  # 'Length'
        self.assertTrue(np.isnan(data.iloc[0, -1]))  # 'Status'

    def test_get(self):
        with mock.patch('proteofav.utils.requests') as mock_get:
            mock_get.get.return_value.ok = True
            mock_get.get.return_value.status = 200
            self.fetch_from_url_or_retry(self.mock_url)
            mock_get.get.assert_called_once_with(mock.ANY, headers=mock.ANY,
                                                 params=mock.ANY, stream=mock.ANY)

    @unittest.expectedFailure
    def test_raise_for_not_found(self):
        with mock.patch('proteofav.utils.requests.get') as mock_get:
            with self.assertRaises(HTTPError) as context:
                mock_get.get.return_value.ok = False
                mock_get.get.return_value.status = 404
                response = self.fetch_from_url_or_retry(self.mock_url)
                mock_get.side_effect = HTTPError(mock.Mock(status=404), 'not found')

        # self.assertTrue('not found' in context.exception)
        self.assertIsNone(response)
        mock_get.get.assert_called_once_with(mock.ANY, headers=mock.ANY, params=mock.ANY)

    def test_mock_a_ptn_sequence(self):
        """
        Test _uniprot_info table parsing with a Mock request.

         ..note:
        query=accession%3AP38995&contact=%20"&columns=id%2Csequence%2C&format=tab
        """

        table = """Entry	Sequence\nP38995	""" + self.ccc2_sequence

        with mock.patch('proteofav.utils.requests') as mock_get:
            mock_get.get.return_value.ok = True
            mock_get.get.return_value.status = 200
            mock_get.get.return_value.content = table
            self.assertEqual(self.get_uniprot_sequence('P38995'), self.ccc2_sequence)

            mock_get.get.assert_called_once_with('http://www.uniprot.org/uniprot/',
                                                 headers=mock.ANY,
                                                 params=mock.ANY,
                                                 stream=mock.ANY)

    def test_mock_a_ptn_organism_name(self):
        """
        Test whether the UniProt info table is parsed correctly by mocking the request.

        ..note:: query=accession%3AP38995&contact=%20"&columns=id%2Corganism%2C&format=tab
        """
        ccc2_organism = 'Saccharomyces cerevisiae'
        table = """Entry	Organism\nP38995	Saccharomyces cerevisiae (strain ATCC 204508 / S288c) (Baker's yeast)"""

        with mock.patch('proteofav.utils.requests') as mock_get:
            mock_get.get.return_value.ok = True
            mock_get.get.return_value.status = 200
            mock_get.get.return_value.content = table
            self.assertEqual(self.get_uniprot_organism('P38995'), ccc2_organism)

            mock_get.get.assert_called_once_with('http://www.uniprot.org/uniprot/',
                                                 headers=mock.ANY,
                                                 params=mock.ANY,
                                                 stream=mock.ANY)

    def test_uniprot_info_with_defaults(self):
        """
        Test whether Uniprot query logic still valid.

        ..note:
         Intentionally get info from Uniprot. An Uniprot update may break this.
        """

        cols = ("Entry", "Entry name", "Status", "Protein names", "Gene names", "Organism",
                "Sequence", "Length")

        first_field = "P38995"
        last_field = 1004

        table = self.uniprot_info('P38995')
        self.assertEqual(set(table.columns), set(cols))
        self.assertEqual(table.iloc[0, 0], first_field)
        self.assertEqual(table.iloc[0, -1], last_field)

    def test_fetch_uniprot_gff(self):
        """
        Test whether Uniprot GFF logic still valid.
        """
        table = self.fetch_uniprot_gff('P38995')

        self.assertTrue(set("NAME SOURCE TYPE START END SCORE STRAND FRAME GROUP".split())
                        .issubset(table.columns))
        self.assertTrue(table['empty'].isnull().all())
        self.assertTrue(table['END'].max() <= len(self.ccc2_sequence))

    def test_map_uniprot_gff(self):
        table = self.map_gff_features_to_sequence('P38995')

        self.assertEqual(list(table.columns), ['annotation'])
        self.assertTrue((~table.annotation.str.contains('Chain')).all())
        self.assertEqual(table.shape[0], len(self.ccc2_sequence))

    def test_map_uniprot_gff_un_grouped(self):
        table = map_gff_features_to_sequence('P38995', group_residues=False)
        self.assertEqual(set(table.columns), {'annotation', 'idx'})
        self.assertTrue(table.idx.max() >= len(self.ccc2_sequence))

    def test_uniprot_to_ensembl_xref(self):
        self.assertTrue(self.uniprot_to_ensembl_xref('P38995').empty)
        table = self.uniprot_to_ensembl_xref('P38995', 'saccharomyces_cerevisiae')
        self.assertEqual(set(table.columns), {'id', 'type'})


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestUNIPROTParser)
    unittest.TextTestRunner(verbosity=2).run(suite)
