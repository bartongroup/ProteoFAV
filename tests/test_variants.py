# -*- coding: utf-8 -*-

import os
import sys
import json
import logging
import unittest
import numpy as np
import unittest.mock as mock

from proteofav.config import Defaults

from proteofav.fetchers import (_fetch_sequence_from_ensembl_protein,
                                fetch_uniprot_variants_ebi,
                                fetch_ensembl_transcript_variants,
                                fetch_uniprot_species_from_id,
                                fetch_ensembl_ensembl_uniprot_mapping,
                                fetch_ensembl_uniprot_ensembl_mapping)

from proteofav.variants import (_match_uniprot_ensembl_seq,
                                _compare_sequences,
                                _count_mismatches,
                                parse_uniprot_variants,
                                Variants,
                                flatten_uniprot_variants_ebi,
                                get_ensembl_species_from_uniprot,
                                get_uniprot_id_from_mapping,
                                get_ensembl_protein_id_from_mapping,
                                get_preferred_uniprot_id_from_mapping,
                                get_preferred_ensembl_id_from_mapping,
                                flatten_ensembl_variants)

from proteofav.utils import (flatten_nested_structure,
                             refactor_key_val_singletons)

defaults = Defaults(os.path.join(os.path.dirname(__file__), "config.txt"))


example_uniprot_variants = {
    'accession': 'P40227',
    'entryName': 'TCPZ_HUMAN',
    'sequence': 'MAAVKTLNPKAEVARAQAAMKATNILLVDEIMRAGMSSLKG',
    'sequenceChecksum': '4876218983604824961',
    'taxid': 9606,
    'features': [
        {'type': 'VARIANT',
         'alternativeSequence': 'F',
         'begin': '231',
         'end': '231',
         'xrefs': [{'name': '1000Genomes',
                    'id': 'rs148616984',
                    'url': 'http://www.ensembl.org/Homo_sapiens/Variation/Explore?v=rs148616984;vdb=variation'},
                   {'name': 'ESP',
                    'id': 'rs148616984',
                    'url': 'http://evs.gs.washington.edu/EVS/PopStatsServlet?searchBy=rsID&target=rs148616984&x=0&y=0'},
                   {'name': 'ExAC',
                    'id': 'rs148616984',
                    'url': 'http://exac.broadinstitute.org/awesome?query=rs148616984'}
                   ],
         'wildType': 'L',
         'frequency': 0.000399361,
         'polyphenPrediction': 'benign',
         'polyphenScore': 0.41,
         'siftPrediction': 'deleterious',
         'siftScore': 0.01,
         'somaticStatus': 0,
         'cytogeneticBand': '7p11.2',
         'consequenceType': 'missense',
         'genomicLocation': 'NC_000007.14:g.56058069C>T',
         'sourceType': 'large_scale_study'},
        {'type': 'VARIANT',
         'alternativeSequence': 'S',
         'begin': '83',
         'end': '83',
         'xrefs': [{'name': 'ExAC',
                    'id': 'rs779741978',
                    'url': 'http://exac.broadinstitute.org/awesome?query=rs779741978'}
                   ],
         'wildType': 'A',
         'frequency': 0.0,
         'polyphenPrediction': 'benign',
         'polyphenScore': 0.305,
         'siftPrediction': 'tolerated',
         'siftScore': 0.06,
         'somaticStatus': 0,
         'cytogeneticBand': '7p11.2',
         'consequenceType': 'missense',
         'genomicLocation': 'NC_000007.14:g.56054414G>T',
         'sourceType': 'large_scale_study'},
        {'type': 'VARIANT',
         'description': '[LSS_COSMIC]: primary tissue(s): lung',
         'alternativeSequence': 'V',
         'begin': '292',
         'end': '292',
         'xrefs': [{'name': 'cosmic curated',
                    'id': 'COSM1549991',
                    'url': 'http://cancer.sanger.ac.uk/cosmic/mutation/overview?id=1549991'}
                   ],
         'evidences': [{'code': 'ECO:0000313',
                        'source': {'name': 'cosmic_study',
                                   'id': 'COSU:417',
                                   'url': 'http://cancer.sanger.ac.uk/cosmic/study/overview?study_id=417'}
                        }
                       ],
         'wildType': 'I',
         'polyphenPrediction': 'benign',
         'polyphenScore': 0.021,
         'siftPrediction': 'tolerated',
         'siftScore': 0.1, 'somaticStatus': 0,
         'cytogeneticBand': '7p11.2',
         'consequenceType': 'missense',
         'genomicLocation': '7:g.56058510A>G',
         'sourceType': 'large_scale_study'}
    ]
}


@mock.patch("proteofav.config.defaults", defaults)
class TestVariants(unittest.TestCase):
    """Test the VAR parser methods."""

    def setUp(self):
        """Initialize the framework for testing."""

        self.uniprotid = "P40227"
        self.uniprotid2 = "P00439"
        self.ensemblid = "ENSP00000448059"
        self.data = example_uniprot_variants
        self.flatten_uniprot_variants_ebi = flatten_uniprot_variants_ebi
        self.get_ensembl_species_from_uniprot = get_ensembl_species_from_uniprot
        self.get_uniprot_id_from_mapping = get_uniprot_id_from_mapping
        self.get_ensembl_protein_id_from_mapping = get_ensembl_protein_id_from_mapping
        self.get_preferred_uniprot_id_from_mapping = get_preferred_uniprot_id_from_mapping
        self.get_preferred_ensembl_id_from_mapping = get_preferred_ensembl_id_from_mapping
        self.flatten_ensembl_variants = flatten_ensembl_variants
        self.vagg = Variants

        self.sequence_from_ensembl_protein = _fetch_sequence_from_ensembl_protein
        self.match_uniprot_ensembl_seq = _match_uniprot_ensembl_seq
        self.compare_sequences = _compare_sequences
        self.count_mismatches = _count_mismatches
        self.parse_uniprot_variants = parse_uniprot_variants

        logging.disable(logging.DEBUG)

    def tearDown(self):
        """Remove testing framework."""

        self.uniprotid = None
        self.uniprotid2 = None
        self.ensemblid = None
        self.data = None
        self.flatten_uniprot_variants_ebi = None
        self.get_ensembl_species_from_uniprot = None
        self.get_uniprot_id_from_mapping = None
        self.get_ensembl_protein_id_from_mapping = None
        self.get_preferred_uniprot_id_from_mapping = None
        self.get_preferred_ensembl_id_from_mapping = None
        self.flatten_ensembl_variants = None
        self.vagg = None

        self.sequence_from_ensembl_protein = None
        self.match_uniprot_ensembl_seq = None
        self.compare_sequences = None
        self.parse_uniprot_variants = None

        logging.disable(logging.NOTSET)

    def test_icgc_parsing(self):
        mock_response = mock.Mock()
        mock_response.return_value.ok = True
        mock_response.return_value.status = 200

        with open('testdata/VARIATION/icgc_ENST00000308677.json') as open_f:
            response = json.load(open_f)

        mock_response.json.return_value = response

        with mock.patch('proteofav.utils.requests.get') as mock_get:
            mock_get.return_value = mock_response
            data = self.fetch_icgc_variants('ENST00000308677')

        self.assertEqual(data.shape, (49, 9))
        self.assertEqual(data.loc[0, 'id'], 'MU30632')
        self.assertEqual(data.loc[84, 'id'], 'MU5450975')
        self.assertEqual(data.loc[84, 'new'], '*')
        self.assertEqual(data.loc[84, 'aaMutation'], 'R46*')

    def test_sequence_from_ensembl_protein(self):
        response = 'MGNAAAAKKGSEQESVKEFLAKAKEDFLKKWESPAQNTAHLDQFER'

        with mock.patch('proteofav.utils.fetch_from_url_or_retry') as mock_get:
            mock_get.return_value = response
            sequence = self.sequence_from_ensembl_protein('XXXXXX')  # ENSP00000309591
        self.assertEqual(sequence, response)

    @mock.patch('proteofav.variants.fetch_uniprot_formal_specie', return_value='homo_sapiens')
    @mock.patch('proteofav.variants.fetch_uniprot_sequence', return_value='ABC')
    @mock.patch('proteofav.variants._fetch_sequence_from_ensembl_protein', return_value='ABC')
    def test_match_uniprot_ensembl_seq(self, mock_fun, mock_fun1, mock_fun2):
        import pandas as pd

        ensembl_response = [{"type": "gene", "id": "ENSG00000072062"}, {"type": "transcript",
                                                                        "id": "ENST00000308677"},
                            {"type": "transcript", "id": "ENST00000589994"},
                            {"type": "translation", "id": "XXXX"}, {"type": "translation",
                                                                    "id": "ENSP00000466651"}]

        with mock.patch('proteofav.variants._uniprot_to_ensembl_xref') as mock_fun3:
            mock_fun3.return_value = pd.DataFrame(ensembl_response)
            self.match_uniprot_ensembl_seq('XXXX')

    @mock.patch('proteofav.variants.fetch_uniprot_formal_specie', return_value='alien')
    def test_match_uniprot_ensembl_seq_unknown_specie(self, mock_fun):
        with self.assertRaises(ValueError):
            self.match_uniprot_ensembl_seq('XXXX')

    def test_compare_sequences_correct(self):
        self.assertTrue(self.compare_sequences('AAA', 'AAA', permissive=False))
        self.assertTrue(self.compare_sequences('AAA', 'AAA', permissive=True))
        self.assertTrue(self.compare_sequences('AAA', 'AAC', n_mismatches=1))

    def test_compare_sequences_correct_bad(self):
        self.assertFalse(self.compare_sequences('AAC', 'AAA'))

    def test_count_mismatches(self):
        self.assertEqual(self.count_mismatches('', ''), 0)
        self.assertEqual(self.count_mismatches('AAA', 'BBB'), 3)
        self.assertEqual(self.count_mismatches('AAA', 'ABA'), 1)
        self.assertEqual(self.count_mismatches('AAA', 'AAAA'), 0)

    def test_parse_uniprot_variants(self):
        import pandas as pd
        # data for P17612
        mock_data = dict(annotation={41: "Natural variant: ['L->V'] (['VAR_040591'])",
                                     46: "Natural variant: ['R->Q'] (['VAR_040592'])",
                                     206: "Natural variant: ['In PPNAD4; somatic mutation; the mutation results in cAMP-independent basal protein kinase activity and constitutive activation of protein kinase A. L->R'] (['VAR_071707'])",
                                     264: "Natural variant: ['S->C'] (['VAR_040593'])"})

        with mock.patch('proteofav.variants.map_gff_features_to_sequence') as mock_fun:
            mock_fun.return_value = pd.DataFrame(mock_data)
            data = self.parse_uniprot_variants('XXXX')

        self.assertEqual(data.shape, (4, 3))
        self.assertFalse('annotation' in data)
        self.assertIn('PPNAD4', data.loc[206, 'disease'])

    def test_flatten_uniprot_variants_ebi(self):
        r = fetch_uniprot_variants_ebi(self.uniprotid2)
        if r.ok:
            table = self.flatten_uniprot_variants_ebi(r)
            self.assertIn('accession', list(table))
            self.assertIn('sequence', list(table))
            self.assertIn('polyphenScore', list(table))
            self.assertIn('siftScore', list(table))
            self.assertIn('begin', list(table))
            self.assertIn('end', list(table))

    def test_flatten_uniprot_variants_ebi_mock(self):
        r = self.flatten_uniprot_variants_ebi(self.data)

        # flattening the accession
        self.assertEqual('P40227', self.data['accession'])
        self.assertEqual('P40227', r.loc[0, 'accession'])

        # flattening the 'xrefs'
        self.assertEqual('1000Genomes', self.data['features'][0]['xrefs'][0]['name'])
        self.assertEqual('1000Genomes', r.loc[0, 'xrefs_name'][0])
        self.assertEqual(['1000Genomes', 'ESP', 'ExAC'], r.loc[0, 'xrefs_name'])
        self.assertTrue(np.isnan(r.loc[1, 'evidences_source_name']))

        # flattening the 'evidences'
        self.assertEqual('cosmic_study', self.data['features'][2]['evidences'][0]['source']['name'])
        self.assertEqual('cosmic_study', r.loc[2, 'evidences_source_name'])

    def test_flatten_json(self):
        modified = {}
        flatten_nested_structure(self.data['features'][0], modified)
        r = refactor_key_val_singletons(modified)
        self.assertIn('type', r)
        self.assertIn('xrefs_name', r)
        self.assertIn('xrefs_id', r)
        self.assertIn('xrefs_url', r)
        self.assertIn('1000Genomes', r['xrefs_name'])
        self.assertIn('rs148616984', r['xrefs_id'])
        self.assertEqual('VARIANT', r['type'])

    def test_get_ensembl_species_from_uniprot(self):
        data = fetch_uniprot_species_from_id(self.uniprotid)
        species = self.get_ensembl_species_from_uniprot(data)
        self.assertEqual(species, 'homo_sapiens')

    def test_get_uniprot_id_from_mapping(self):
        data = fetch_ensembl_ensembl_uniprot_mapping(self.ensemblid)

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
        data = fetch_ensembl_uniprot_ensembl_mapping(self.uniprotid2)

        r = self.get_ensembl_protein_id_from_mapping(data.json())
        self.assertEqual(r, [self.ensemblid])

    def test_preferred_uniprot_id_from_mapping(self):
        info = fetch_ensembl_ensembl_uniprot_mapping(self.ensemblid)
        data = get_uniprot_id_from_mapping(info.json(), full_entry=True)
        best_match = get_preferred_uniprot_id_from_mapping(data)
        self.assertEqual(best_match, 'P00439')

    def test_preferred_ensembl_id_from_mapping(self):
        data = fetch_uniprot_species_from_id(self.uniprotid)
        species = self.get_ensembl_species_from_uniprot(data)
        info = fetch_ensembl_uniprot_ensembl_mapping(self.uniprotid2,
                                                     species=species)
        r = self.get_ensembl_protein_id_from_mapping(info.json())
        best_match = self.get_preferred_ensembl_id_from_mapping(r)
        self.assertEqual(best_match, 'ENSP00000448059')

    def test_flatten_ensembl_variants(self):
        r = fetch_ensembl_transcript_variants(self.ensemblid)
        if r.ok:
            table = self.flatten_ensembl_variants(r, synonymous=True)
            self.assertIn('polyphenScore', list(table))
            self.assertIn('siftScore', list(table))
            self.assertIn('begin', list(table))
            self.assertIn('end', list(table))

    def test_variants_agreggator_uniprot(self):
        v = self.vagg(self.uniprotid2)
        table = v.run(uniprot_vars=True, synonymous=True,
                      ensembl_transcript_vars=True,
                      ensembl_somatic_vars=True)

        self.assertIn('polyphenScore', list(table))
        self.assertIn('siftScore', list(table))
        self.assertIn('begin', list(table))
        self.assertIn('end', list(table))

    def test_variants_agreggator_ensembl(self):
        v = self.vagg(self.ensemblid, uniprot=False)
        table = v.run(uniprot_vars=True, synonymous=True,
                      ensembl_transcript_vars=False,
                      ensembl_somatic_vars=False)

        self.assertIn('polyphenScore', list(table))
        self.assertIn('siftScore', list(table))
        self.assertIn('begin', list(table))
        self.assertIn('end', list(table))


if __name__ == '__main__':
    logging.basicConfig(stream=sys.stderr)
    logging.getLogger("proteofav").setLevel(logging.DEBUG)
    suite = unittest.TestLoader().loadTestsFromTestCase(TestVariants)
    unittest.TextTestRunner(verbosity=2).run(suite)
