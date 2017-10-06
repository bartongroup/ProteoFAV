# -*- coding: utf-8 -*-

import sys
import logging
import unittest

try:
    import mock
except ImportError:
    import unittest.mock as mock

from proteofav.annotation import (_fetch_uniprot_gff, map_gff_features_to_sequence)


class TestUNIPROTParser(unittest.TestCase):
    """Test UniProt fetcher/parser methods."""

    def setUp(self):
        """Initialize the framework for testing."""

        self.fetch_uniprot_gff = _fetch_uniprot_gff
        self.map_gff_features_to_sequence = map_gff_features_to_sequence
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

        self.fetch_uniprot_gff = None
        self.map_gff_features_to_sequence = None
        self.ccc2_sequence = None

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


if __name__ == '__main__':
    logging.basicConfig(stream=sys.stderr)
    logging.getLogger("proteofav.config").setLevel(logging.CRITICAL)
    suite = unittest.TestLoader().loadTestsFromTestCase(TestUNIPROTParser)
    unittest.TextTestRunner(verbosity=2).run(suite)
