import json
import unittest
from os import path

try:
    import mock
except ImportError:
    import unittest.mock as mock

from proteofav.config import Defaults
from proteofav.variants import (_fetch_icgc_variants, _fetch_ebi_variants, _fetch_ensembl_variants,
                                _sequence_from_ensembl_protein, _match_uniprot_ensembl_seq,
                                _compare_sequences, _count_mismatches, parse_uniprot_variants)

defaults = Defaults(path.join(path.dirname(__file__), "config.txt"))


@mock.patch("proteofav.structures.defaults", defaults)
class VariantsTestCase(unittest.TestCase):
    """Test for the variants.py"""

    def setUp(self):
        """Initialize the framework for testing."""
        self.fetch_icgc_variants = _fetch_icgc_variants
        self.fetch_ebi_variants = _fetch_ebi_variants
        self.fetch_ensembl_variants = _fetch_ensembl_variants
        self.sequence_from_ensembl_protein = _sequence_from_ensembl_protein
        self.match_uniprot_ensembl_seq = _match_uniprot_ensembl_seq
        self.compare_sequences = _compare_sequences
        self.count_mismatches = _count_mismatches
        self.parse_uniprot_variants = parse_uniprot_variants

    def tearDown(self):
        """Remove testing framework by cleaning the namespace."""
        self.fetch_icgc_variants = None
        self.fetch_ebi_variants = None
        self.fetch_ensembl_variants = None
        self.sequence_from_ensembl_protein = None
        self.match_uniprot_ensembl_seq = None
        self.compare_sequences = None
        self.parse_uniprot_variants = None

    def test_icgc_parsing(self):
        mock_response = mock.Mock()
        mock_response.return_value.ok = True
        mock_response.return_value.status = 200

        with open('testdata/variation/icgc_ENST00000308677.json') as open_f:
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

    def test_sequence_from_ensembl_protein(self):
        response = 'MGNAAAAKKGSEQESVKEFLAKAKEDFLKKWESPAQNTAHLDQFER'

        with mock.patch('proteofav.variants.fetch_from_url_or_retry') as mock_get:
            mock_get.return_value = response
            sequence = self.sequence_from_ensembl_protein('XXXXXX')  # ENSP00000309591
        self.assertEqual(sequence, response)

    @mock.patch('proteofav.variants.fetch_uniprot_formal_specie', return_value='homo_sapiens')
    @mock.patch('proteofav.variants.fetch_uniprot_sequence', return_value='ABC')
    @mock.patch('proteofav.variants._sequence_from_ensembl_protein', return_value='ABC')
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

        self.assertEqual(data.shape, (4,3))
        self.assertFalse('annotation' in data)
        self.assertIn('PPNAD4', data.loc[206, 'disease'])




if __name__ == '__main__':
    unittest.main()
