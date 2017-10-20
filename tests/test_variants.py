# -*- coding: utf-8 -*-

import os
import sys
import json
import logging
import unittest
import numpy as np
import pandas as pd

try:
    import mock
except ImportError:
    import unittest.mock as mock

from proteofav.config import defaults
from proteofav.variants import (_fetch_icgc_variants, fetch_uniprot_variants, fetch_ensembl_variants,
                                fetch_ensembl_sequence_from_id, _match_uniprot_ensembl_seq,
                                _compare_sequences, _count_mismatches, parse_uniprot_variants,
                                fetch_uniprot_sequence, fetch_uniprot_id_from_name,
                                fetch_uniprot_ensembl_mapping, fetch_ensembl_uniprot_mapping,
                                get_uniprot_id_from_mapping, get_ensembl_protein_id_from_mapping,
                                fetch_uniprot_pdb_mapping, fetch_pdb_uniprot_mapping,
                                get_preferred_uniprot_id_from_mapping,
                                get_preferred_ensembl_id_from_mapping,
                                fetch_uniprot_species_from_id, get_ensembl_species_from_uniprot,
                                fetch_uniprot_formal_specie, _uniprot_info,
                                flatten_uniprot_variants_ebi, flatten_ensembl_variants,
                                select_variants, fetch_variants, Variants)

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


@mock.patch("proteofav.structures.defaults", defaults)
class VariantsTestCase(unittest.TestCase):
    """Test for the variants.py"""

    def setUp(self):
        """Initialize the framework for testing."""
        self.uniprot_id = 'O96013'
        self.pdbid = '2pah'
        self.uniprot_id2 = 'P00439'
        self.ensembl_id = 'ENSP00000326864'
        self.ensembl_id2 = "ENSP00000448059"
        self.variant_id = 'rs557625940'
        self.fetch_icgc_variants = _fetch_icgc_variants
        self.fetch_uniprot_variants = fetch_uniprot_variants
        self.fetch_ensembl_variants = fetch_ensembl_variants
        self.fetch_ensembl_sequence_from_id = fetch_ensembl_sequence_from_id
        self.match_uniprot_ensembl_seq = _match_uniprot_ensembl_seq
        self.compare_sequences = _compare_sequences
        self.count_mismatches = _count_mismatches
        self.parse_uniprot_variants = parse_uniprot_variants
        self.uniprot_info = _uniprot_info
        self.fetch_uniprot_ensembl_mapping = fetch_uniprot_ensembl_mapping
        self.fetch_ensembl_uniprot_mapping = fetch_ensembl_uniprot_mapping
        self.get_preferred_uniprot_id_from_mapping = get_preferred_uniprot_id_from_mapping
        self.get_preferred_ensembl_id_from_mapping = get_preferred_ensembl_id_from_mapping
        self.get_uniprot_id_from_mapping = get_uniprot_id_from_mapping
        self.get_ensembl_protein_id_from_mapping = get_ensembl_protein_id_from_mapping
        self.fetch_uniprot_species_from_id = fetch_uniprot_species_from_id
        self.fetch_uniprot_id_from_name = fetch_uniprot_id_from_name
        self.get_ensembl_species_from_uniprot = get_ensembl_species_from_uniprot
        self.get_uniprot_sequence = fetch_uniprot_sequence
        self.get_uniprot_organism = fetch_uniprot_formal_specie
        self.uniprot_info = _uniprot_info
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
        self.fetch_pdb_uniprot_mapping = fetch_pdb_uniprot_mapping
        self.fetch_uniprot_pdb_mapping = fetch_uniprot_pdb_mapping
        self.flatten_uniprot_variants_ebi = flatten_uniprot_variants_ebi
        self.flatten_ensembl_variants = flatten_ensembl_variants
        self.fetch_variants = fetch_variants
        self.select_variants = select_variants
        self.Variants = Variants

    def tearDown(self):
        """Remove testing framework by cleaning the namespace."""

        self.uniprot_id = None
        self.pdbid = None
        self.uniprot_id2 = None
        self.ensembl_id = None
        self.ensembl_id2 = None
        self.variant_id = None
        self.fetch_icgc_variants = None
        self.fetch_uniprot_variants = None
        self.fetch_ensembl_variants = None
        self.fetch_ensembl_sequence_from_id = None
        self.match_uniprot_ensembl_seq = None
        self.compare_sequences = None
        self.parse_uniprot_variants = None
        self.uniprot_info = None
        self.fetch_uniprot_ensembl_mapping = None
        self.fetch_ensembl_uniprot_mapping = None
        self.get_uniprot_id_from_mapping = None
        self.get_ensembl_protein_id_from_mapping = None
        self.get_preferred_uniprot_id_from_mapping = None
        self.get_preferred_ensembl_id_from_mapping = None
        self.fetch_uniprot_species_from_id = None
        self.fetch_uniprot_id_from_name = None
        self.get_ensembl_species_from_uniprot = None
        self.get_uniprot_sequence = None
        self.get_uniprot_organism = None
        self.uniprot_info = None
        self.ccc2_sequence = None
        self.fetch_pdb_uniprot_mapping = None
        self.fetch_uniprot_pdb_mapping = None
        self.flatten_uniprot_variants_ebi = None
        self.flatten_ensembl_variants = None
        self.fetch_variants = None
        self.select_variants = None
        self.Variants = None

    def test_icgc_parsing(self):
        mock_response = mock.Mock()
        mock_response.return_value.ok = True
        mock_response.return_value.status = 200

        with open(os.path.join(os.path.dirname(__file__), "testdata",
                               "variation", "icgc_ENST00000308677.json")) as open_f:
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

    def test_ebi_variants_parsing(self):
        raw_response = """{"accession":"P17612", "entryName":"KAPCA_HUMAN", "sequence":
        "MGNAAAAKKGSEQESVKEFLAKAKEDFLKKWESPAQNTAHLDQFERIKTLGTGSFGRVMLVKHKETGNHYAMKILDKQKVVKLKQIEHTLNEKRILQAVNFPFLVKLEFSFKDNSNLYMVMEYVPGGEMFSHLRRIGRFSEPHARFYAAQIVLTFEYLHSLDLIYRDLKPENLLIDQQGYIQVTDFGFAKRVKGRTWTLCGTPEYLAPEIILSKGYNKAVDWWALGVLIYEMAAGYPPFFADQPIQIYEKIVSGKVRFPSHFSSDLKDLLRNLLQVDLTKRFGNLKNGVNDIKNHKWFATTDWIAIYQRKVEAPFIPKFKGPGDTSNFDDYEEEEIRVSINEKCGKEFSEF",                    "sequenceChecksum":"13793750284533818795", "taxid":9606,                   "features":[{"type":"VARIANT","ftId":"VAR_040591","alternativeSequence":"V","begin":"41","end":"41","xrefs":[{"name":"dbSNP","id":"rs56029020","url":"http://www.ncbi.nlm.nih.gov/SNP/snp_ref.cgi?type=rs&rs=rs56029020"},{"name":"Ensembl","id":"rs56029020","url":"http://www.ensembl.org/id/rs56029020"}],"wildType":"L","somaticStatus":0,"consequenceType":"missense","sourceType":"uniprot"},{"type":"VARIANT","alternativeSequence":"I","begin":"252","end":"252","xrefs":[{"name":"ExAC","id":"rs760535486","url":"http://exac.broadinstitute.org/awesome?query=rs760535486"}],"wildType":"V","polyphenPrediction":"benign","polyphenScore":0.025,"siftPrediction":"tolerated","siftScore":0.21,"somaticStatus":0,"cytogeneticBand":"19p13.12","consequenceType":"missense", "genomicLocation":"NC_000019.10:g.14097372C>T","sourceType":"large_scale_study"}]} """
        mock_response = mock.Mock()
        mock_response.return_value.ok = True
        mock_response.return_value.status = 200
        mock_response.json.return_value = json.loads(raw_response)

        with mock.patch('proteofav.utils.requests.get') as mock_get:
            mock_get.return_value = mock_response
            data = self.fetch_uniprot_variants('P17612')

            from pandas.io.json import json_normalize
            data = json_normalize(data.json(), ['features'],
                                  meta=['accession', 'entryName'])

            # flatten the xref field, which has the id column.
            # ideally this could be normalised with:
            # table = json_normalize(data, ['features', 'xref'] ...
            # but this field is not present in all entries,
            flat_xref = data['xrefs'].apply(pd.Series).stack().apply(pd.Series)
            flat_xref.reset_index(level=1, drop=True, inplace=True)
            data = data.join(flat_xref)

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
            data = pd.DataFrame(data.json())

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
            data = pd.DataFrame(data.json())

        self.assertEqual(data.shape, (2, 15))
        self.assertEqual(data.loc[0, 'Parent'], 'ENST00000288602')
        self.assertEqual(data.loc[0, 'allele'], 'COSMIC_MUTATION')
        self.assertEqual(data.loc[0, 'start'], 433)
        self.assertEqual(data.loc[0, 'type'], 'coding_sequence_variant')
        self.assertEqual(data.loc[1, 'clinical_significance'], [])
        self.assertEqual(data.loc[1, 'type'], 'coding_sequence_variant')

    def test_sequence_from_ensembl_protein(self):
        response = 'MGNAAAAKKGSEQESVKEFLAKAKEDFLKKWESPAQNTAHLDQFER'

        with mock.patch('proteofav.variants.fetch_from_url_or_retry') as mock_get:
            mock_get.return_value = response
            sequence = self.fetch_ensembl_sequence_from_id('XXXXXX')  # ENSP00000309591
        self.assertEqual(sequence, response)

    @mock.patch('proteofav.variants.fetch_uniprot_formal_specie', return_value='homo_sapiens')
    @mock.patch('proteofav.variants.fetch_uniprot_sequence', return_value='ABC')
    @mock.patch('proteofav.variants.fetch_ensembl_sequence_from_id', return_value='ABC')
    def test_match_uniprot_ensembl_seq(self, mock_fun, mock_fun1, mock_fun2):
        import pandas as pd

        ensembl_response = """[{"type": "gene", "id": "ENSG00000072062"}, {"type": "transcript",
                                                                            "id": "ENST00000308677"},
                                {"type": "transcript", "id": "ENST00000589994"},
                                {"type": "translation", "id": "XXXX"}, {"type": "translation",
                                                                        "id": "ENSP00000466651"}]"""
        mock_response = mock.Mock()
        mock_response.return_value.ok = True
        mock_response.return_value.status = 200
        mock_response.json.return_value = json.loads(ensembl_response)
        with mock.patch('proteofav.variants.fetch_uniprot_ensembl_mapping') as mock_fun3:
            mock_fun3.return_value = mock_response
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

        with mock.patch('proteofav.variants.select_annotation') as mock_fun:
            mock_fun.return_value = pd.DataFrame(mock_data)
            data = self.parse_uniprot_variants('XXXX')

        self.assertEqual(data.shape, (4, 3))
        self.assertFalse('annotation' in data)
        self.assertIn('PPNAD4', data.loc[206, 'disease'])

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

    def test_obsolete_uniprot_accession(self):
        data = self.uniprot_info('Q91887')
        self.assertEqual(len(data), 1)

        # number of keys (or columns)
        self.assertEqual(len(data.columns.values), 8)

        # check the values for particular entries
        self.assertTrue(np.isnan(data.iloc[0, -1]))  # 'Length'
        self.assertTrue(np.isnan(data.iloc[0, -1]))  # 'Status'

    def test_uniprot_variants_ebi(self):
        r = self.fetch_uniprot_variants(self.uniprot_id)
        self.assertTrue(r.ok)

    def test_fetch_ensembl_transcript_variants(self):
        r = self.fetch_ensembl_variants(self.ensembl_id,
                                        feature="transcript_variation")
        self.assertTrue(r.ok)

    def test_fetch_ensembl_somatic_variants(self):
        r = self.fetch_ensembl_variants(self.ensembl_id,
                                        feature="somatic_transcript_variation")
        self.assertTrue(r.ok)

    def fetch_fetch_ensembl_sequence_from_id(self):
        r = self.fetch_ensembl_sequence_from_id(self.ensembl_id)
        self.assertTrue(r.ok)

    def test_fetch_uniprot_species_from_id(self):
        r = self.fetch_uniprot_species_from_id(self.uniprot_id2)
        self.assertTrue(r.ok)
        organism = str(r.content, encoding='utf-8').split('\n')[1]
        species = '_'.join(organism.split()[0:2]).lower()
        self.assertEqual(species, "homo_sapiens")

    def test_fetch_uniprot_id_from_name(self):
        r = self.fetch_uniprot_id_from_name("PH4H_HUMAN")
        self.assertTrue(r.ok)
        self.assertEqual("P00439", str(r.content, encoding='utf-8').strip())

    def test_fetch_fetch_uniprot_pdb_mapping(self):
        r = self.fetch_uniprot_pdb_mapping(self.uniprot_id)
        self.assertTrue(r.ok)

    def test_fetch_fetch_pdb_uniprot_mapping(self):
        r = self.fetch_pdb_uniprot_mapping(self.pdbid)
        self.assertTrue(r.ok)

    def test_get_ensembl_species_from_uniprot(self):
        data = self.fetch_uniprot_species_from_id(self.uniprot_id2)
        species = self.get_ensembl_species_from_uniprot(data)
        self.assertEqual(species, 'homo_sapiens')

    def test_fetch_uniprot_ensembl_mapping(self):
        r = self.fetch_uniprot_ensembl_mapping(self.uniprot_id2)
        self.assertTrue(r.ok)
        ensps = self.get_ensembl_protein_id_from_mapping(r.json())
        self.assertEqual(ensps, [self.ensembl_id2])

    def test_fetch_ensembl_uniprot_mapping(self):
        r = self.fetch_ensembl_uniprot_mapping(self.ensembl_id2)
        self.assertTrue(r.ok)
        uniprots = self.get_uniprot_id_from_mapping(r.json())
        self.assertEqual(uniprots, ['A0A024RBG4', self.uniprot_id2])

    def test_to_table_uniprot_ensembl_mapping_full(self):
        """
        Tests the fetching and parsing real UniProt ids.
        This test focuses on the method that parses the residue entries.

        Some checks are made to whether the parsed keys and values
        are the ones we are expecting.
        """

        data = self.fetch_uniprot_ensembl_mapping(self.uniprot_id,
                                                  species='homo_sapiens')

        information = {}
        rows = []
        for entry in data.json():
            typ = entry['type'].upper()
            eid = entry['id']
            try:
                if eid in information[typ]:
                    continue
                information[typ].append(eid)
            except KeyError:
                information[typ] = eid
            except AttributeError:
                information[typ] = [information[typ]]
                information[typ].append(eid)

        rows.append(information)
        data = pd.DataFrame(rows)

        # number of values per column (or rows)
        self.assertEqual(len(data), 1)

        # number of keys (or columns)
        self.assertEqual(len(data.columns.values), 3)

        # check whether there are particular keys
        self.assertIn('TRANSCRIPT', data.columns.values)

        # check the values for particular entries
        self.assertEqual(data['GENE'][0], 'ENSG00000130669')

    def test_get_uniprot_id_from_mapping(self):
        data = self.fetch_ensembl_uniprot_mapping(self.ensembl_id2)

        r = self.get_uniprot_id_from_mapping(data.json(), full_entry=False,
                                             uniprot_id=None)
        self.assertEqual(r, ['A0A024RBG4', self.uniprot_id2])

        r = self.get_uniprot_id_from_mapping(data.json(), full_entry=True,
                                             uniprot_id=None)
        self.assertIn('dbname', r[0])
        self.assertIn('xref_start', r[0])
        self.assertIn('ensembl_start', r[0])

        r = self.get_uniprot_id_from_mapping(data.json(), full_entry=False,
                                             uniprot_id=self.uniprot_id2)
        self.assertEqual(r, [self.uniprot_id2])

    def test_get_ensembl_protein_id_from_mapping(self):
        data = self.fetch_uniprot_ensembl_mapping(self.uniprot_id2)
        r = self.get_ensembl_protein_id_from_mapping(data.json())
        self.assertEqual(r, [self.ensembl_id2])

    def test_preferred_uniprot_id_from_mapping(self):
        info = self.fetch_ensembl_uniprot_mapping(self.ensembl_id2)
        data = self.get_uniprot_id_from_mapping(info.json(), full_entry=True)
        best_match = self.get_preferred_uniprot_id_from_mapping(data)
        self.assertEqual(best_match, 'P00439')

    def test_preferred_ensembl_id_from_mapping(self):
        data = self.fetch_uniprot_species_from_id(self.uniprot_id2)
        species = self.get_ensembl_species_from_uniprot(data)
        info = self.fetch_uniprot_ensembl_mapping(self.uniprot_id2,
                                                  species=species)
        r = self.get_ensembl_protein_id_from_mapping(info.json())
        best_match = self.get_preferred_ensembl_id_from_mapping(r)
        self.assertEqual(best_match, 'ENSP00000448059')

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

    def test_uniprot_to_ensembl_xref(self):
        self.assertEqual(self.fetch_uniprot_ensembl_mapping('P38995').json(), [])
        table = self.fetch_uniprot_ensembl_mapping('P38995', 'saccharomyces_cerevisiae')
        table = pd.DataFrame(table.json())
        self.assertEqual(set(table.columns), {'id', 'type'})

    def test_to_table_pdb_uniprot_sifts_mapping(self):
        """
        Testing the PDBe API for mapping between PDB and UniProt
        accession identifiers.
        """

        information = self.fetch_pdb_uniprot_mapping(self.pdbid)
        rows = []
        for uniprot in information.json()[self.pdbid]['UniProt']:
            uniprots = {'uniprot_id': uniprot}
            rows.append(uniprots)
        data = pd.DataFrame(rows)

        # number of values per column (or rows)
        self.assertEqual(len(data), 1)

        # number of keys (or columns)
        self.assertEqual(len(data.columns.values), 1)

        # check whether there are particular keys
        self.assertIn('uniprot_id', data.columns.values)

        # check the values of particular entries
        self.assertTrue(data['uniprot_id'][0] == 'P00439')

    def test_to_table_uniprot_pdb_sifts_mapping(self):
        """
        Testing the PDBe API for mapping between UniProt and PDB
        accession identifiers.
        """

        information = self.fetch_uniprot_pdb_mapping(self.uniprot_id2)
        rows = []
        for entry in information.json()[self.uniprot_id2]:
            rows.append(entry)
        data = pd.DataFrame(rows)

        # check whether there are particular keys
        self.assertIn('pdb_id', data.columns.values)

        # check the values of particular entries
        self.assertEqual(data['pdb_id'].unique()[0], '2pah')
        self.assertIn('A', data['chain_id'].unique())
        self.assertIn('X-ray diffraction', data['experimental_method'].unique())
        self.assertTrue(type(data['coverage'][0]), float)
        self.assertTrue(type(data['resolution'][0]), float)
        self.assertTrue(type(data['tax_id'][0]), int)

    def test_to_table_transcript_variants_ensembl(self):
        """
        Tests the fetching and parsing real Ensembl Protein ids.

        Some checks are made to whether the parsed keys and values
        are the ones we are expecting.

        Using mocking because variant data varies over time.
        For simplicity here loading the resulting pandas table from a
        dumped csv file.
        """

        # data = self.ensembl_trascript(self.ensembl_id, species='human')
        # mocking the resulting pandas dataframe
        variants = os.path.join(os.path.dirname(__file__), "testdata",
                                "variation", "transcript_variation_{}.csv".format(self.ensembl_id))
        data = pd.read_csv(variants)

        # number of values per column (or rows)
        self.assertEqual(len(data), 206)

        # number of keys (or columns)
        self.assertEqual(len(data.columns.values), 16)

        # check whether there are particular keys
        self.assertIn('Parent', data.columns.values)

        # check the values for particular entries
        self.assertTrue(data['Parent'][0] == 'ENST00000321944')
        self.assertTrue(data['feature_type'][0] == "transcript_variation")
        self.assertEqual(data['codons'][0], 'Cat/Tat')
        self.assertEqual(data['start'][0], 135)

    def test_to_table_somatic_variants_ensembl(self):
        """
        Tests the fetching and parsing real Ensembl Protein ids.

        Some checks are made to whether the parsed keys and values
        are the ones we are expecting.

        Using mocking because variant data varies over time.
        For simplicity here loading the resulting pandas table from a
        dumped csv file.
        """

        # data = self.ensembl_somatic(self.ensembl_id, species='human')
        # mocking the resulting pandas dataframe
        variants = os.path.join(os.path.dirname(__file__), "testdata", "variation",
                                "somatic_variation_{}.csv".format(self.ensembl_id))
        data = pd.read_csv(variants)

        # number of values per column (or rows)
        self.assertEqual(len(data), 40)

        # number of keys (or columns)
        self.assertEqual(len(data.columns.values), 16)

        # check whether there are particular keys
        self.assertIn('Parent', data.columns.values)

        # check the values for particular entries
        self.assertTrue(data['Parent'][0] == 'ENST00000321944')
        self.assertTrue(data['feature_type'][0] == "somatic_transcript_variation")
        self.assertEqual(data['codons'][0], 'Gag/Cag')
        self.assertEqual(data['start'][0], 119)

    def test_to_table_ensembl_variant(self):
        """
        Tests the fetching and parsing real Ensembl Variant ids.

        Some checks are made to whether the parsed keys and values
        are the ones we are expecting.

        Using mocking because variant data varies over time.
        For simplicity here loading the resulting pandas table from a
        dumped csv file.
        """

        # data = self.ensembl_variant(self.variant_id, species='human')
        # mocking the resulting pandas dataframe
        variants = os.path.join(os.path.dirname(__file__), "testdata", "variation",
                                "ensembl_variation_{}.csv".format(self.variant_id))
        data = pd.read_csv(variants)

        # number of values per column (or rows)
        self.assertEqual(len(data), 1)

        # number of keys (or columns)
        self.assertEqual(len(data.columns.values), 19)

        # check whether there are particular keys
        self.assertIn('location', data.columns.values)

        # check the values for particular entries
        self.assertEqual(data['location'][0], '19:39175331-39175331')
        self.assertTrue(data['name'][0] == "rs557625940")
        self.assertEqual(data['seq_region_name'][0], 19)
        self.assertEqual(data['most_severe_consequence'][0], 'missense_variant')

    def test_to_table_uniprot_ensembl_variants(self):
        """
        Tests the wrapper method that goes from a uniprot to a list
        of ensembl protein transcripts and to variants.

        Using mocking because variant data varies over time.
        For simplicity here loading the resulting pandas table from a
        dumped csv file.
        """

        # querying the various ensembl endpoints
        # data = self.uniprot_variants(self.uniprot_id)
        # mocking the resulting pandas dataframe
        variants = os.path.join(os.path.dirname(__file__), "testdata", "variation",
                                "uniprot_variants_{}.csv".format(self.uniprot_id))
        data = pd.read_csv(variants)

        # number of values per column (or rows)
        self.assertEqual(len(data), 903)

        # number of keys (or columns)
        self.assertEqual(len(data.columns.values), 5)

        # check whether there are particular keys
        self.assertIn('translation', data.columns.values)
        self.assertIn('id', data.columns.values)
        self.assertIn('start', data.columns.values)
        self.assertIn('residues', data.columns.values)

        # check the values for particular entries
        self.assertEqual(data['translation'][0], 'ENSP00000351049')
        self.assertTrue(data['id'][0] == "rs769098772")
        self.assertTrue(data['start'][0] == 295)
        self.assertTrue(data['residues'][0] == "E/A")

    def test_flatten_uniprot_variants_ebi(self):
        r = self.fetch_uniprot_variants(self.uniprot_id2)
        if r.ok:
            table = self.flatten_uniprot_variants_ebi(r)
            self.assertIn('accession', list(table))
            self.assertIn('sequence', list(table))
            self.assertIn('polyphenScore', list(table))
            self.assertIn('siftScore', list(table))
            self.assertIn('begin', list(table))
            self.assertIn('end', list(table))

    def test_flatten_uniprot_variants_ebi_mock(self):
        data = example_uniprot_variants
        r = self.flatten_uniprot_variants_ebi(data)

        # flattening the accession
        self.assertEqual('P40227', data['accession'])
        self.assertEqual('P40227', r.loc[0, 'accession'])

        # flattening the 'xrefs'
        self.assertEqual('1000Genomes', data['features'][0]['xrefs'][0]['name'])
        self.assertEqual('1000Genomes', r.loc[0, 'xrefs_name'][0])
        self.assertEqual(['1000Genomes', 'ESP', 'ExAC'], r.loc[0, 'xrefs_name'])
        self.assertTrue(np.isnan(r.loc[1, 'evidences_source_name']))

        # flattening the 'evidences'
        self.assertEqual('cosmic_study', data['features'][2]['evidences'][0]['source']['name'])
        self.assertEqual('cosmic_study', r.loc[2, 'evidences_source_name'])

    def test_flatten_ensembl_variants(self):
        r = self.fetch_ensembl_variants(self.ensembl_id2, feature='transcript_variation')
        if r.ok:
            table = self.flatten_ensembl_variants(r, synonymous=True)
            self.assertIn('polyphenScore', list(table))
            self.assertIn('siftScore', list(table))
            self.assertIn('begin', list(table))
            self.assertIn('end', list(table))

    def test_fetch_variants_all(self):
        uniprot, germline, somatic = self.fetch_variants(self.uniprot_id2,
                                                         id_source='uniprot',
                                                         synonymous=True,
                                                         uniprot_vars=True,
                                                         ensembl_germline_vars=True,
                                                         ensembl_somatic_vars=True)

        self.assertEqual(uniprot.loc[0, 'accession'], self.uniprot_id2)
        self.assertEqual(germline.loc[0, 'translation'], self.ensembl_id2)
        self.assertEqual(somatic.loc[0, 'translation'], self.ensembl_id2)
        self.assertIn('polyphenScore', list(uniprot))
        self.assertIn('siftScore', list(uniprot))
        self.assertIn('begin', list(germline))
        self.assertIn('begin', list(somatic))
        self.assertIn('end', list(germline))
        self.assertIn('end', list(somatic))

    def test_variants_fetch_from_uniprot_id(self):
        uniprot, germline, somatic = self.Variants.fetch(self.uniprot_id2,
                                                         id_source='uniprot',
                                                         synonymous=True,
                                                         uniprot_vars=True,
                                                         ensembl_germline_vars=True,
                                                         ensembl_somatic_vars=True)
        self.assertEqual(uniprot.loc[0, 'accession'], self.uniprot_id2)
        self.assertEqual(germline.loc[0, 'translation'], self.ensembl_id2)
        self.assertEqual(somatic.loc[0, 'translation'], self.ensembl_id2)
        self.assertIn('polyphenScore', list(uniprot))
        self.assertIn('siftScore', list(uniprot))
        self.assertIn('begin', list(germline))
        self.assertIn('begin', list(somatic))
        self.assertIn('end', list(germline))
        self.assertIn('end', list(somatic))

    def test_variants_fetch_from_ensembl_id(self):
        uniprot, germline, somatic = self.Variants.fetch(self.ensembl_id2,
                                                         id_source='ensembl',
                                                         synonymous=True,
                                                         uniprot_vars=True,
                                                         ensembl_germline_vars=True,
                                                         ensembl_somatic_vars=True)
        self.assertEqual(uniprot.loc[0, 'accession'], self.uniprot_id2)
        self.assertEqual(germline.loc[0, 'translation'], self.ensembl_id2)
        self.assertEqual(somatic.loc[0, 'translation'], self.ensembl_id2)
        self.assertIn('polyphenScore', list(uniprot))
        self.assertIn('siftScore', list(uniprot))
        self.assertIn('begin', list(germline))
        self.assertIn('begin', list(somatic))
        self.assertIn('end', list(germline))
        self.assertIn('end', list(somatic))

    def test_variants_select(self):
        uniprot, ensembl = self.select_variants(self.uniprot_id2,
                                                id_source='uniprot',
                                                synonymous=True,
                                                uniprot_vars=True,
                                                ensembl_germline_vars=True,
                                                ensembl_somatic_vars=True)
        self.assertIn('polyphenScore', list(uniprot))
        self.assertIn('siftScore', list(uniprot))
        self.assertIn('begin', list(uniprot))
        self.assertIn('end', list(uniprot))
        self.assertIn('accession', list(uniprot))
        self.assertIn('translation', list(ensembl))
        self.assertIn('begin', list(ensembl))
        self.assertIn('end', list(ensembl))

    def test_variants_select_from_uniprot_id(self):
        uniprot, ensembl = self.Variants.select(self.uniprot_id2,
                                                id_source='uniprot',
                                                synonymous=True,
                                                uniprot_vars=False,
                                                ensembl_germline_vars=True,
                                                ensembl_somatic_vars=True)
        self.assertIn('polyphenScore', list(ensembl))
        self.assertIn('siftScore', list(ensembl))
        self.assertIn('begin', list(ensembl))
        self.assertIn('end', list(ensembl))
        self.assertNotIn('accession', list(ensembl))
        self.assertIn('translation', list(ensembl))
        self.assertIsNone(uniprot)

    def test_variants_select_from_ensembl_id(self):
        uniprot, ensembl = self.Variants.select(self.ensembl_id2,
                                                id_source='ensembl',
                                                synonymous=True,
                                                uniprot_vars=True,
                                                ensembl_germline_vars=False,
                                                ensembl_somatic_vars=False)
        self.assertIn('polyphenScore', list(uniprot))
        self.assertIn('siftScore', list(uniprot))
        self.assertIn('begin', list(uniprot))
        self.assertIn('end', list(uniprot))
        self.assertIn('accession', list(uniprot))
        self.assertNotIn('translation', list(uniprot))
        self.assertIsNone(ensembl)


if __name__ == '__main__':
    logging.basicConfig(stream=sys.stderr)
    logging.getLogger("proteofav.config").setLevel(logging.CRITICAL)
    suite = unittest.TestLoader().loadTestsFromTestCase(VariantsTestCase)
    unittest.TextTestRunner(verbosity=2).run(suite)
