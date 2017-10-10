# -*- coding: utf-8 -*-

import os
import sys
import logging
import unittest

try:
    import mock
except ImportError:
    import unittest.mock as mock

from proteofav.annotation import (parse_gff_features, select_annotation,
                                  annotation_aggregation, filter_annotation,
                                  download_annotation, Annotation)


class TestUNIPROTParser(unittest.TestCase):
    """Test UniProt fetcher/parser methods."""

    def setUp(self):
        """Initialize the framework for testing."""

        self.uniprotid = "P38995"
        self.example_annotation = os.path.join(os.path.dirname(__file__), "testdata",
                                               "annotation", "P38995.gff")
        self.output_annotation = os.path.join(os.path.dirname(__file__), "testdata",
                                              "P38995.gff")
        self.parse_gff_features = parse_gff_features
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
        self.annotation_aggregation = annotation_aggregation
        self.filter_annotation = filter_annotation
        self.download_annotation = download_annotation
        self.Annotation = Annotation

    def tearDown(self):
        """Remove testing framework."""

        self.uniprotid = None
        self.example_annotation = None
        self.output_annotation = None
        self.map_gff_features_to_sequence = None
        self.ccc2_sequence = None
        self.annotation_aggregation = None
        self.filter_annotation = None
        self.download_annotation = None
        self.Annotation = None

    def test_parse_uniprot_gff(self):
        """
        Test whether Uniprot GFF logic still valid.
        """
        table = self.parse_gff_features(self.example_annotation)
        self.assertTrue(set("NAME SOURCE TYPE START END SCORE STRAND FRAME GROUP".split())
                        .issubset(table.columns))
        self.assertNotIn('empty', table)
        self.assertTrue(table['END'].max() <= len(self.ccc2_sequence))

    def test_filter_annotation_uniprot_gff(self):
        table = self.parse_gff_features(self.example_annotation)
        table = self.annotation_aggregation(table)
        self.assertEqual(set(table.columns), {'annotation', 'site', 'accession'})
        self.assertTrue((~table.annotation.str.contains('Chain')).all())
        self.assertEqual(table.shape[0], len(self.ccc2_sequence))

        table = self.parse_gff_features(self.example_annotation)
        table = self.filter_annotation(table, annotation_agg=True)
        self.assertEqual(set(table.columns), {'annotation', 'site', 'accession'})
        self.assertTrue((~table.annotation.str.contains('Chain')).all())
        self.assertEqual(table.shape[0], len(self.ccc2_sequence))

    def test_parse_uniprot_gff_un_grouped(self):
        table = self.parse_gff_features(self.example_annotation)
        table = self.annotation_aggregation(table, group_residues=False)
        self.assertEqual(set(table.columns), {'annotation', 'idx', 'site', 'accession'})
        self.assertTrue(table.idx.max() >= len(self.ccc2_sequence))

        table = self.parse_gff_features(self.example_annotation)
        table = self.filter_annotation(table, annotation_agg=True, group_residues=False)
        self.assertEqual(set(table.columns), {'annotation', 'idx', 'site', 'accession'})
        self.assertTrue(table.idx.max() >= len(self.ccc2_sequence))

    def test_download_uniprot_gff(self):
        self.download_annotation(self.uniprotid, filename=self.output_annotation,
                                 overwrite=True)
        if os.path.exists(self.output_annotation):
            os.remove(self.output_annotation)

    def test_main_Annotation(self):
        # read
        table = self.Annotation.read(self.example_annotation)
        self.assertTrue(set("NAME SOURCE TYPE START END SCORE STRAND FRAME GROUP".split())
                        .issubset(table.columns))
        self.assertNotIn('empty', table)
        self.assertTrue(table['END'].max() <= len(self.ccc2_sequence))
        # download
        self.Annotation.download(self.uniprotid, filename=self.output_annotation,
                                 overwrite=True)
        if os.path.exists(self.output_annotation):
            os.remove(self.output_annotation)
        # select
        self.Annotation.select(self.uniprotid)


if __name__ == '__main__':
    logging.basicConfig(stream=sys.stderr)
    logging.getLogger("proteofav.config").setLevel(logging.CRITICAL)
    suite = unittest.TestLoader().loadTestsFromTestCase(TestUNIPROTParser)
    unittest.TextTestRunner(verbosity=2).run(suite)
