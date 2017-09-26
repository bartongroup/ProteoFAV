# -*- coding: utf-8 -*-

import sys
import logging
import unittest
import requests_cache

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
                                get_ensembl_species_from_uniprot)

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


if __name__ == '__main__':
    logging.basicConfig(stream=sys.stderr)
    logging.getLogger("proteofav").setLevel(logging.DEBUG)
    suite = unittest.TestLoader().loadTestsFromTestCase(TestFetchers)
    unittest.TextTestRunner(verbosity=2).run(suite)
