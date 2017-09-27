# -*- coding: utf-8 -*-

import sys
import json
import logging
import unittest
import requests_cache
from unittest import mock

from proteofav.fetchers import (Fetcher,
                                fetch_best_structures_pdbe,
                                fetch_uniprot_variants_ebi,
                                fetch_summary_properties_pdbe,
                                fetch_uniprot_id_from_name,
                                fetch_uniprot_species_from_id,
                                fetch_ensembl_uniprot_ensembl_mapping,
                                fetch_ensembl_ensembl_uniprot_mapping,
                                fetch_ensembl_transcript_variants,
                                fetch_ensembl_somatic_variants,
                                fetch_ensembl_variants_by_id,
                                fetch_ensembl_sequence_from_id,
                                get_preferred_assembly_id,

                                get_ensembl_protein_id_from_mapping,
                                get_uniprot_id_from_mapping,
                                get_preferred_uniprot_id_from_mapping,
                                get_preferred_ensembl_id_from_mapping,
                                get_ensembl_species_from_uniprot,

                                _fetch_sequence_from_ensembl_protein,
                                _fetch_icgc_variants,
                                _fetch_ebi_variants,
                                _fetch_ensembl_variants)

from proteofav.config import defaults as config


class TestFetchers(unittest.TestCase):
    """Test the dssp parser methods."""

    def setUp(self):
        """Initialize the framework for testing."""

        requests_cache.uninstall_cache()

        self.uniprotid = "P00439"
        self.uniprotid2 = "P00439"
        self.pdbid = "2pah"
        self.cathid = "1.50.10.100_1318"
        self.pfamid = "PF08124"
        self.ensemblid = "ENSP00000448059"
        self.varid = "rs750420403"
        self.Fetcher = Fetcher
        self.fetch_best_structures_pdbe = fetch_best_structures_pdbe
        self.fetch_summary_properties_pdbe = fetch_summary_properties_pdbe
        self.fetch_preferred_assembly_id = get_preferred_assembly_id
        self.fetch_uniprot_variants_ebi = fetch_uniprot_variants_ebi
        self.fetch_uniprot_id_from_name = fetch_uniprot_id_from_name
        self.fetch_uniprot_species_from_id = fetch_uniprot_species_from_id
        self.fetch_ensembl_uniprot_ensembl_mapping = fetch_ensembl_uniprot_ensembl_mapping
        self.fetch_ensembl_ensembl_uniprot_mapping = fetch_ensembl_ensembl_uniprot_mapping
        self.fetch_ensembl_transcript_variants = fetch_ensembl_transcript_variants
        self.fetch_ensembl_somatic_variants = fetch_ensembl_somatic_variants
        self.fetch_ensembl_variants_by_id = fetch_ensembl_variants_by_id
        self.fetch_ensembl_sequence_from_id = fetch_ensembl_sequence_from_id

        self.get_ensembl_protein_id_from_mapping = get_ensembl_protein_id_from_mapping
        self.get_uniprot_id_from_mapping = get_uniprot_id_from_mapping
        self.get_preferred_uniprot_id_from_mapping = get_preferred_uniprot_id_from_mapping
        self.get_preferred_ensembl_id_from_mapping = get_preferred_ensembl_id_from_mapping
        self.get_ensembl_species_from_uniprot = get_ensembl_species_from_uniprot

        self.fetch_icgc_variants = _fetch_icgc_variants
        self.fetch_ebi_variants = _fetch_ebi_variants
        self.fetch_ensembl_variants = _fetch_ensembl_variants

        logging.disable(logging.DEBUG)

    def tearDown(self):
        """Remove testing framework."""

        self.uniprotid = None
        self.uniprotid2 = None
        self.pdbid = None
        self.cathid = None
        self.pfamid = None
        self.ensemblid = None
        self.varid = None
        self.Fetcher = None
        self.fetch_best_structures_pdbe = None
        self.fetch_summary_properties_pdbe = None
        self.fetch_preferred_assembly_id = None
        self.fetch_uniprot_variants_ebi = None
        self.fetch_uniprot_id_from_name = None
        self.fetch_uniprot_species_from_id = None
        self.fetch_ensembl_uniprot_ensembl_mapping = None
        self.fetch_ensembl_ensembl_uniprot_mapping = None
        self.fetch_ensembl_transcript_variants = None
        self.fetch_ensembl_somatic_variants = None
        self.fetch_ensembl_variants_by_id = None
        self.fetch_ensembl_sequence_from_id = None

        self.get_ensembl_protein_id_from_mapping = None
        self.get_uniprot_id_from_mapping = None
        self.get_preferred_uniprot_id_from_mapping = None
        self.get_preferred_ensembl_id_from_mapping = None
        self.get_ensembl_species_from_uniprot = None

        self.fetch_icgc_variants = None
        self.fetch_ebi_variants = None
        self.fetch_ensembl_variants = None

        logging.disable(logging.NOTSET)

    def test_fetcher_generic(self):
        url_root = config.http_uniprot
        url_endpoint = self.uniprotid + ".fasta"
        url = url_root + url_endpoint
        b = self.Fetcher(url, cached=False, json=False)
        r = b.response
        self.assertTrue(r.ok)

    def test_best_structure_pdbe(self):
        r = self.fetch_best_structures_pdbe(self.uniprotid)
        self.assertTrue(r.ok)

    def test_summary_properties_pdbe(self):
        r = self.fetch_summary_properties_pdbe(self.pdbid)
        self.assertTrue(r.ok)

    def test_preferred_assembly_pdbe(self):
        r = self.fetch_preferred_assembly_id(self.pdbid)
        self.assertEqual("1", r)

    def test_uniprot_variants_ebi(self):
        r = self.fetch_uniprot_variants_ebi(self.uniprotid)
        self.assertTrue(r.ok)

    def test_fetch_uniprot_id_from_name(self):
        r = self.fetch_uniprot_id_from_name("PH4H_HUMAN")
        self.assertTrue(r.ok)
        self.assertEqual("P00439", str(r.content, encoding='utf-8').strip())

    def test_fetch_uniprot_species_from_id(self):
        r = self.fetch_uniprot_species_from_id(self.uniprotid)
        self.assertTrue(r.ok)
        organism = str(r.content, encoding='utf-8').split('\n')[1]
        species = '_'.join(organism.split()[0:2]).lower()
        self.assertEqual(species, "homo_sapiens")

    def test_fetch_ensembl_uniprot_ensembl_mapping(self):
        r = self.fetch_ensembl_uniprot_ensembl_mapping(self.uniprotid)
        self.assertTrue(r.ok)
        ensps = self.get_ensembl_protein_id_from_mapping(r.json())
        self.assertEqual(ensps, [self.ensemblid])

    def test_fetch_ensembl_ensembl_uniprot_mapping(self):
        r = self.fetch_ensembl_ensembl_uniprot_mapping(self.ensemblid)
        self.assertTrue(r.ok)
        uniprots = self.get_uniprot_id_from_mapping(r.json())
        self.assertEqual(uniprots, ['A0A024RBG4', self.uniprotid])

    def test_fetch_ensembl_transcript_variants(self):
        r = self.fetch_ensembl_transcript_variants(self.ensemblid)
        self.assertTrue(r.ok)

    def fetch_fetch_ensembl_somatic_variants(self):
        r = self.fetch_ensembl_somatic_variants(self.ensemblid)
        self.assertTrue(r.ok)

    def fetch_fetch_ensembl_variants_by_id(self):
        r = self.fetch_ensembl_variants_by_id(self.varid)
        self.assertTrue(r.ok)

    def fetch_fetch_ensembl_sequence_from_id(self):
        r = self.fetch_ensembl_sequence_from_id(self.ensemblid)
        self.assertTrue(r.ok)

    def test_get_ensembl_species_from_uniprot(self):
        data = self.fetch_uniprot_species_from_id(self.uniprotid)
        species = self.get_ensembl_species_from_uniprot(data)
        self.assertEqual(species, 'homo_sapiens')

    def test_get_uniprot_id_from_mapping(self):
        data = self.fetch_ensembl_ensembl_uniprot_mapping(self.ensemblid)

        r = self.get_uniprot_id_from_mapping(data.json(), full_entry=False,
                                             uniprot_id=None)
        self.assertEqual(r, ['A0A024RBG4', self.uniprotid2])

        r = self.get_uniprot_id_from_mapping(data.json(), full_entry=True,
                                             uniprot_id=None)
        self.assertIn('dbname', r[0])
        self.assertIn('xref_start', r[0])
        self.assertIn('ensembl_start', r[0])

        r = self.get_uniprot_id_from_mapping(data.json(), full_entry=False,
                                             uniprot_id=self.uniprotid2)
        self.assertEqual(r, [self.uniprotid2])

    def test_get_ensembl_protein_id_from_mapping(self):
        data = self.fetch_ensembl_uniprot_ensembl_mapping(self.uniprotid2)

        r = self.get_ensembl_protein_id_from_mapping(data.json())
        self.assertEqual(r, [self.ensemblid])

    def test_preferred_uniprot_id_from_mapping(self):
        info = self.fetch_ensembl_ensembl_uniprot_mapping(self.ensemblid)
        data = self.get_uniprot_id_from_mapping(info.json(), full_entry=True)
        best_match = self.get_preferred_uniprot_id_from_mapping(data)
        self.assertEqual(best_match, 'P00439')

    def test_preferred_ensembl_id_from_mapping(self):
        data = self.fetch_uniprot_species_from_id(self.uniprotid)
        species = self.get_ensembl_species_from_uniprot(data)
        info = self.fetch_ensembl_uniprot_ensembl_mapping(self.uniprotid2,
                                                          species=species)
        r = self.get_ensembl_protein_id_from_mapping(info.json())
        best_match = self.get_preferred_ensembl_id_from_mapping(r)
        self.assertEqual(best_match, 'ENSP00000448059')

    def test_fetch_ebi_variants_parsing(self):
        raw_response = """{"accession":"P17612", "entryName":"KAPCA_HUMAN", "sequence":
        "MGNAAAAKKGSEQESVKEFLAKAKEDFLKKWESPAQNTAHLDQFERIKTLGTGSFGRVMLVKHKETGNHYAMKILDKQKVVKLKQIEHTLNEKRILQAVNFPFLVKLEFSFKDNSNLYMVMEYVPGGEMFSHLRRIGRFSEPHARFYAAQIVLTFEYLHSLDLIYRDLKPENLLIDQQGYIQVTDFGFAKRVKGRTWTLCGTPEYLAPEIILSKGYNKAVDWWALGVLIYEMAAGYPPFFADQPIQIYEKIVSGKVRFPSHFSSDLKDLLRNLLQVDLTKRFGNLKNGVNDIKNHKWFATTDWIAIYQRKVEAPFIPKFKGPGDTSNFDDYEEEEIRVSINEKCGKEFSEF",                    "sequenceChecksum":"13793750284533818795", "taxid":9606,                   "features":[{"type":"VARIANT","ftId":"VAR_040591","alternativeSequence":"V","begin":"41","end":"41","xrefs":[{"name":"dbSNP","id":"rs56029020","url":"http://www.ncbi.nlm.nih.gov/SNP/snp_ref.cgi?type=rs&rs=rs56029020"},{"name":"Ensembl","id":"rs56029020","url":"http://www.ensembl.org/id/rs56029020"}],"wildType":"L","somaticStatus":0,"consequenceType":"missense","sourceType":"uniprot"},{"type":"VARIANT","alternativeSequence":"I","begin":"252","end":"252","xrefs":[{"name":"ExAC","id":"rs760535486","url":"http://exac.broadinstitute.org/awesome?query=rs760535486"}],"wildType":"V","polyphenPrediction":"benign","polyphenScore":0.025,"siftPrediction":"tolerated","siftScore":0.21,"somaticStatus":0,"cytogeneticBand":"19p13.12","consequenceType":"missense", "genomicLocation":"NC_000019.10:g.14097372C>T","sourceType":"large_scale_study"}]} """
        mock_response = mock.Mock()
        mock_response.return_value.ok = True
        mock_response.return_value.status = 200
        mock_response.json.return_value = json.loads(raw_response)

        with mock.patch('proteofav.utils.requests.get') as mock_get:
            mock_get.return_value = mock_response
            data = self.fetch_ebi_variants('P17612')

        self.assertEqual(data.shape, (3, 21))
        # deal with duplicated index - two data sources for the first variant
        data = data.reset_index().drop_duplicates(subset='index', keep='last').set_index(
            'index')
        self.assertEqual(data.shape, (2, 21))
        self.assertEqual(data.loc[0, 'begin'], "41")
        # TODO check why it not parsing thois one as float
        self.assertEqual(data.loc[0, 'end'], "41")
        self.assertEqual(data.loc[0, 'consequenceType'], 'missense')
        self.assertEqual(data.loc[0, 'sourceType'], 'uniprot')
        self.assertEqual(data.loc[1, 'sourceType'], 'large_scale_study')
        self.assertEqual(data.loc[0, 'id'], 'rs56029020')
        self.assertEqual(data.loc[1, 'id'], 'rs760535486')

    def test_fetch_ensembl_variants_transcript_variation_parsing(self):
        raw_response = """[{"polyphen":0.908,"sift":0.06,"feature_type":"transcript_variation","clinical_significance":[],"Parent":
"ENST00000288602","codons":"Cca/Gca","end":622,"seq_region_name":"ENSP00000288602","residues":"P/A","minor_allele_frequency":
null,"id":"rs746074624","translation":"ENSP00000288602","allele":"G/C","type":"missense_variant","start":622}
,{"polyphen":0.95,"sift":0,"feature_type":"transcript_variation","clinical_significance":["pathogenic"],"Parent"
:"ENST00000288602","codons":"Gca/Cca","end":246,"seq_region_name":"ENSP00000288602","residues":"A/P",
"minor_allele_frequency":null,"id":"rs180177034","translation":"ENSP00000288602","allele":"C/G","type":"missense_variant","start":246}]"""
        mock_response = mock.Mock()
        mock_response.return_value.ok = True
        mock_response.return_value.status = 200
        mock_response.json.return_value = json.loads(raw_response)

        with mock.patch('proteofav.utils.requests.get') as mock_get:
            mock_get.return_value = mock_response
            data = self.fetch_ensembl_variants('XXXXX',  # ENSP00000309591
                                               feature='transcript_variation')

        self.assertEqual(data.shape, (2, 15))
        self.assertEqual(data.loc[0, 'Parent'], 'ENST00000288602')
        self.assertEqual(data.loc[0, 'end'], 622)
        self.assertEqual(data.loc[1, 'sift'], 0)
        self.assertEqual(data.loc[1, 'start'], 246)
        self.assertEqual(data.loc[1, 'type'], 'missense_variant')

    def test_fetch_ensembl_variants_somatic_transcript_variation_parsing(self):
        raw_response = """[{"polyphen":null,"sift":null,
        "feature_type":"somatic_transcript_variation","clinical_significance":[],
        "Parent":"ENST00000288602","codons":"","end":433,
        "seq_region_name":"ENSP00000288602","residues":"","minor_allele_frequency":null,
        "id":"COSM3832072","translation":"ENSP00000288602","allele":"COSMIC_MUTATION",
        "type":"coding_sequence_variant","start":433},{"polyphen":null,"sift":null,
        "feature_type":"somatic_transcript_variation","clinical_significance":[],
        "Parent":"ENST00000288602","codons":"","end":698,
        "seq_region_name":"ENSP00000288602","residues":"","minor_allele_frequency":null,
        "id":"COSM452456","translation":"ENSP00000288602","allele":"COSMIC_MUTATION",
        "type":"coding_sequence_variant","start":698}] """
        mock_response = mock.Mock()
        mock_response.return_value.ok = True
        mock_response.return_value.status = 200
        mock_response.json.return_value = json.loads(raw_response)

        with mock.patch('proteofav.utils.requests.get') as mock_get:
            mock_get.return_value = mock_response
            data = self.fetch_ensembl_variants('XXXXX',  # ENSP00000309591
                                               feature='somatic_transcript_variation')

        self.assertEqual(data.shape, (2, 15))
        self.assertEqual(data.loc[0, 'Parent'], 'ENST00000288602')
        self.assertEqual(data.loc[0, 'allele'], 'COSMIC_MUTATION')
        self.assertEqual(data.loc[0, 'start'], 433)
        self.assertEqual(data.loc[0, 'type'], 'coding_sequence_variant')
        self.assertEqual(data.loc[1, 'clinical_significance'], [])
        self.assertEqual(data.loc[1, 'type'], 'coding_sequence_variant')


if __name__ == '__main__':
    logging.basicConfig(stream=sys.stderr)
    logging.getLogger("proteofav").setLevel(logging.DEBUG)
    suite = unittest.TestLoader().loadTestsFromTestCase(TestFetchers)
    unittest.TextTestRunner(verbosity=2).run(suite)
