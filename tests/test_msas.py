# -*- coding: utf-8 -*-

import os
import sys
import logging
import unittest

from Bio.Align import MultipleSeqAlignment

try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO
try:
    from mock import patch
except ImportError:
    from unittest.mock import patch

from proteofav.msas import (read_alignments, read_msas,
                            parse_sequence_info_from_description,
                            parse_uniprot_fasta_seq_description,
                            parse_pfam_sth_seq_description,
                            parse_cath_fasta_seq_description,
                            parse_cath_sth_seq_description,
                            parse_generic_seq_description,
                            download_msa_from_cath,
                            download_msa_from_pfam,
                            download_msas, select_msas, MSA)

from proteofav.config import defaults

root = os.path.abspath(os.path.dirname(__file__))
defaults.db_msas = os.path.join(root, "testdata", "msas")


@patch("proteofav.dssp.defaults", defaults)
class TestMSAS(unittest.TestCase):
    """Test the MSAs parser methods."""

    def setUp(self):
        """Initialize the framework for testing."""

        self.cathid = '1.50.10.100_1318'
        self.cathid2 = '1.20.1070.10_7072'
        self.pfamid = 'PF00118'

        self.uniprot_fasta = ("tr|A0A067NRR9|A0A067NRR9_PLEOS Polysaccharide lyase "
                              "family 8 protein OS=Pleurotus ostreatus PC15 "
                              "GN=PLEOSDRAFT_1036297 PE=4 SV=1")
        self.pfam_sto = "A0A067NRR9_PLEOS/22-313"
        self.cath_sto = "3capB00/1-326"
        self.cath_fasta = ("cath|4.1.0|1rwhA01/4-372 CATH_S35=1.50.10.100.1;"
                           "GENE=P84141_ARTAU;GO=GO:0005576,GO:0005975,GO:0016837,"
                           "GO:0030246;MDA=1.50.10.100;ORG=Arthrobacter_aurescens;"
                           "TAXON=43663;UNIPROT=P84141")
        self.cath_fasta2 = ("biomap|4.1.0|08901c114dd3dfa7eb9d39cf55143681/343-713 "
                            "EC=4.2.2.1;GENE=Q8VLQ7_STRSU;GO=GO:0005576,GO:0005975,GO:"
                            "0016020,GO:0030246,GO:0030340;MDA=2.60.120.260_2.60.40.1380_"
                            "1.50.10.100_2.70.98.10_2.60.220.10;ORG=Streptococcus_suis;"
                            "TAXON=1307;UNIPROT=Q8VLQ7")
        self.generic_id = "A0A067NRR9/22-313"

        self.inputcath = os.path.join(defaults.db_msas,
                                      "{}.fasta".format(self.cathid))
        self.inputcath2 = os.path.join(defaults.db_msas,
                                       "{}.sth".format(self.cathid2))
        self.inputpfam = os.path.join(defaults.db_msas,
                                      "{}.sth".format(self.pfamid))
        self.output_fasta = os.path.join(defaults.db_msas, "tmp.fasta")
        self.output_sto = os.path.join(defaults.db_msas, "tmp.sto")
        self.notfound = ""
        self.excluded = ()
        self.entry = {}

        self.read_alignments = read_alignments
        self.read_msas = read_msas

        self.uniprot_fasta_seq_description = parse_uniprot_fasta_seq_description
        self.pfam_sth_seq_description = parse_pfam_sth_seq_description
        self.cath_fasta_seq_description = parse_cath_fasta_seq_description
        self.generic_seq_description = parse_generic_seq_description
        self.sequence_info_from_description = parse_sequence_info_from_description
        self.cath_sth_seq_description = parse_cath_sth_seq_description

        self.download_msa_from_cath = download_msa_from_cath
        self.download_msa_from_pfam = download_msa_from_pfam
        self.download_msas = download_msas
        self.select_msas = select_msas
        self.MSA = MSA

        logging.disable(logging.DEBUG)

    def tearDown(self):
        """Remove testing framework."""

        self.cathid = None
        self.cathid2 = None
        self.pfamid = None

        self.uniprot_fasta = None
        self.pfam_sto = None
        self.cath_sto = None
        self.cath_fasta = None
        self.generic_id = None

        self.inputcath = None
        self.inputcath2 = None
        self.inputpfam = None
        self.output_fasta = None
        self.output_sto = None

        self.notfound = None
        self.excluded = None
        self.entry = None

        self.read_alignments = None
        self.read_msas = None

        self.uniprot_fasta_seq_description = None
        self.pfam_sth_seq_description = None
        self.cath_fasta_seq_description = None
        self.generic_seq_description = None
        self.sequence_info_from_description = None

        self.download_msa_from_cath = None
        self.download_msa_from_pfam = None
        self.download_msas = None
        self.select_msas = None
        self.MSA = None

        logging.disable(logging.NOTSET)

    def test_file_not_found_reader(self):
        with self.assertRaises(IOError):
            self.read_alignments(self.notfound)

    def test_file_not_found_parser(self):
        with self.assertRaises(IOError):
            self.read_msas(self.notfound)

    def test_read_alignment(self):
        align, seq_format = self.read_alignments(self.inputcath)
        self.assertTrue(type(align), MultipleSeqAlignment)
        self.assertEqual(len(align), 42)
        self.assertEqual(seq_format, "fasta")
        align, seq_format = self.read_alignments(self.inputpfam)
        self.assertTrue(type(align), MultipleSeqAlignment)
        self.assertEqual(len(align), 57)
        self.assertEqual(seq_format, "stockholm")
        align, seq_format = self.read_alignments(self.inputcath2)
        self.assertTrue(type(align), MultipleSeqAlignment)
        self.assertEqual(len(align), 59)
        self.assertEqual(seq_format, "stockholm")

    def test_parse_uniprot_fasta_seq_description(self):
        self.uniprot_fasta_seq_description(self.uniprot_fasta,
                                           self.entry)
        self.assertEqual(self.entry['Source'], "UniProt")
        self.assertEqual(self.entry['Collection'], "tr")
        self.assertEqual(self.entry['Accession'], "A0A067NRR9")
        self.assertEqual(self.entry['Name'], "A0A067NRR9_PLEOS")

    def test_parse_pfam_sth_seq_description(self):
        self.pfam_sth_seq_description(self.pfam_sto, self.entry,
                                      get_uniprot_id=True)
        self.assertEqual(self.entry['Source'], "Pfam")
        self.assertEqual(self.entry['Accession'], "A0A067NRR9")
        self.assertEqual(self.entry['Name'], "A0A067NRR9_PLEOS")
        self.assertEqual(self.entry['Start'], 22)
        self.assertEqual(self.entry['End'], 313)

    def test_parse_cath_sth_seq_description(self):
        self.cath_sth_seq_description(self.cath_sto, self.entry,
                                      get_uniprot_id=False)
        self.assertEqual(self.entry['Source'], "CATH")
        self.assertEqual(self.entry['Accession'], "3capB00")
        self.assertEqual(self.entry['pdb_id'], "3cap")
        self.assertEqual(self.entry['chain_id'], "B")
        self.assertEqual(self.entry['domain_id'], "3capB00")
        self.assertEqual(self.entry['Start'], 1)
        self.assertEqual(self.entry['End'], 326)

        self.entry = {}
        self.cath_sth_seq_description(self.cath_sto, self.entry,
                                      get_uniprot_id=True)
        self.assertEqual(self.entry['Source'], "CATH")
        self.assertEqual(self.entry['Accession'], "P02699")
        self.assertEqual(self.entry['pdb_id'], "3cap")
        self.assertEqual(self.entry['chain_id'], "B")
        self.assertEqual(self.entry['domain_id'], "3capB00")
        self.assertEqual(self.entry['Start'], 1)
        self.assertEqual(self.entry['End'], 326)

    def test_parse_cath_fasta_seq_description(self):
        self.cath_fasta_seq_description(self.cath_fasta,
                                        self.entry)
        self.assertEqual(self.entry['Source'], "CATH")
        self.assertEqual(self.entry['Collection'], "cath")
        self.assertEqual(self.entry['Version'], "4.1.0")
        self.assertEqual(self.entry['Accession'], "P84141")
        self.assertEqual(self.entry['Domain'], "1rwhA01")
        self.assertEqual(self.entry['Start'], 4)
        self.assertEqual(self.entry['End'], 372)

    def test_parse_cath_fasta_seq_description_2(self):
        self.cath_fasta_seq_description(self.cath_fasta2,
                                        self.entry)
        self.assertEqual(self.entry['Source'], "CATH")
        self.assertEqual(self.entry['Collection'], "biomap")
        self.assertEqual(self.entry['Version'], "4.1.0")
        self.assertEqual(self.entry['Accession'],
                         "08901c114dd3dfa7eb9d39cf55143681")
        self.assertEqual(self.entry['Start'], 343)
        self.assertEqual(self.entry['End'], 713)

    def test_parse_generic_seq_description(self):
        self.generic_seq_description(self.generic_id,
                                     self.entry)
        self.assertEqual(self.entry['Source'], "GenericParser")
        self.assertEqual(self.entry['Accession'], "A0A067NRR9")
        self.assertEqual(self.entry['Start'], 22)
        self.assertEqual(self.entry['End'], 313)

    def test_parse_sequence_info_from_description(self):
        self.entry = {}
        self.sequence_info_from_description(self.uniprot_fasta,
                                            self.entry, seq_format="fasta")
        self.assertEqual(self.entry['Source'], "UniProt")
        self.assertEqual(self.entry['Collection'], "tr")
        self.assertEqual(self.entry['Accession'], "A0A067NRR9")
        self.assertEqual(self.entry['Name'], "A0A067NRR9_PLEOS")

        self.entry = {}
        self.sequence_info_from_description(self.pfam_sto, self.entry,
                                            seq_format="stockholm",
                                            get_uniprot_id=True)
        self.assertEqual(self.entry['Source'], "Pfam")
        self.assertEqual(self.entry['Accession'], "A0A067NRR9")
        self.assertEqual(self.entry['Name'], "A0A067NRR9_PLEOS")
        self.assertEqual(self.entry['Start'], 22)
        self.assertEqual(self.entry['End'], 313)

        self.entry = {}
        self.sequence_info_from_description(self.cath_fasta,
                                            self.entry,
                                            seq_format="fasta")
        self.assertEqual(self.entry['Source'], "CATH")
        self.assertEqual(self.entry['Collection'], "cath")
        self.assertEqual(self.entry['Version'], "4.1.0")
        self.assertEqual(self.entry['Accession'], "P84141")
        self.assertEqual(self.entry['Domain'], "1rwhA01")
        self.assertEqual(self.entry['Start'], 4)
        self.assertEqual(self.entry['End'], 372)

        self.entry = {}
        self.sequence_info_from_description(self.cath_sto,
                                            self.entry,
                                            seq_format="stockholm")
        self.assertEqual(self.entry['Source'], "CATH")
        self.assertEqual(self.entry['Accession'], "P02699")
        self.assertEqual(self.entry['pdb_id'], "3cap")
        self.assertEqual(self.entry['chain_id'], "B")
        self.assertEqual(self.entry['domain_id'], "3capB00")
        self.assertEqual(self.entry['Start'], 1)
        self.assertEqual(self.entry['End'], 326)

        self.entry = {}
        self.sequence_info_from_description(self.generic_id,
                                            self.entry)
        self.assertEqual(self.entry['Source'], "GenericParser")
        self.assertEqual(self.entry['Accession'], "A0A067NRR9")
        self.assertEqual(self.entry['Start'], 22)
        self.assertEqual(self.entry['End'], 313)

    def test_read_alignment_to_table_cath(self):
        data = self.read_msas(filename=self.inputcath,
                              excluded_cols=self.excluded)
        self.assertIn('Sequence', list(data))
        self.assertIn('Source', list(data))
        self.assertIn('Description', list(data))
        self.assertIn('Collection', list(data))
        self.assertEqual(data.loc[0, 'Accession'], 'Q59288')
        self.assertEqual(data.loc[0, 'Domain'], '1hm3A01')
        self.assertEqual(data.loc[0, 'Collection'], 'cath')
        self.assertEqual(data.loc[0, 'Start'], 27)
        self.assertEqual(data.loc[0, 'End'], 338)

    def test_read_alignment_to_table_pfam(self):
        data = self.read_msas(filename=self.inputpfam,
                              excluded_cols=self.excluded,
                              get_uniprot_id=True)
        self.assertIn('Sequence', list(data))
        self.assertIn('Source', list(data))
        self.assertIn('Name', list(data))
        self.assertIn('Accession', list(data))
        self.assertEqual(data.loc[0, 'Accession'], 'B9LRY6')
        self.assertEqual(data.loc[0, 'Source'], 'Pfam')
        self.assertEqual(data.loc[0, 'Start'], 27)
        self.assertEqual(data.loc[0, 'End'], 514)

    def test_download_msa_from_cath(self):
        self.download_msa_from_cath(self.cathid, self.output_fasta,
                                    seq_format="fasta",
                                    overwrite=True)
        self.assertTrue(os.path.isfile(self.output_fasta))
        os.remove(self.output_fasta)

        self.download_msa_from_cath(self.cathid, self.output_sto,
                                    seq_format="stockholm",
                                    overwrite=True)
        self.assertTrue(os.path.isfile(self.output_sto))
        os.remove(self.output_sto)

    def test_download_msa_from_pfam(self):
        self.download_msa_from_pfam(self.pfamid, self.output_sto,
                                    overwrite=True)
        self.assertTrue(os.path.isfile(self.output_sto))
        os.remove(self.output_sto)

    def test_download_msas(self):
        self.download_msas(self.cathid, self.output_fasta,
                           aln_source="cath", seq_format="fasta",
                           overwrite=True)
        self.assertTrue(os.path.isfile(self.output_fasta))
        os.remove(self.output_fasta)

        self.download_msas(self.cathid, self.output_sto,
                           aln_source="cath", seq_format="stockholm",
                           overwrite=True)
        self.assertTrue(os.path.isfile(self.output_sto))
        os.remove(self.output_sto)

        self.download_msas(self.pfamid, self.output_sto,
                           aln_source="pfam", seq_format="stockholm",
                           overwrite=True)
        self.assertTrue(os.path.isfile(self.output_sto))
        os.remove(self.output_sto)

    def test_select_msas(self):
        data = self.select_msas(superfamily=self.pfamid,
                                family=None, seq_format="stockholm",
                                aln_source='pfam', get_uniprot_id=True,
                                excluded_cols=None, overwrite=False)
        self.assertIn('Sequence', list(data))
        self.assertIn('Source', list(data))
        self.assertIn('Name', list(data))
        self.assertIn('Accession', list(data))
        self.assertEqual(data.loc[0, 'Accession'], 'B9LRY6')
        self.assertEqual(data.loc[0, 'Source'], 'Pfam')
        self.assertEqual(data.loc[0, 'Start'], 27)
        self.assertEqual(data.loc[0, 'End'], 514)

    def test_reader_class_cath(self):
        # read
        data = self.MSA.read(filename=self.inputcath, excluded_cols=self.excluded,
                             get_uniprot_id=True)
        self.assertIn('Sequence', list(data))
        self.assertIn('Source', list(data))
        self.assertIn('Description', list(data))
        self.assertIn('Collection', list(data))
        self.assertEqual(data.loc[0, 'Accession'], 'Q59288')
        self.assertEqual(data.loc[0, 'Domain'], '1hm3A01')
        self.assertEqual(data.loc[0, 'Collection'], 'cath')
        self.assertEqual(data.loc[0, 'Start'], 27)
        self.assertEqual(data.loc[0, 'End'], 338)

        # download
        self.MSA.download(self.pfamid, self.output_sto,
                          aln_source="pfam", seq_format="stockholm",
                          overwrite=True)
        self.assertTrue(os.path.isfile(self.output_sto))
        os.remove(self.output_sto)

        # select
        data = self.MSA.select(superfamily=self.pfamid,
                               family=None, seq_format="stockholm",
                               aln_source='pfam', get_uniprot_id=True,
                               excluded_cols=None, overwrite=False)
        self.assertIn('Sequence', list(data))
        self.assertIn('Source', list(data))
        self.assertIn('Name', list(data))
        self.assertIn('Accession', list(data))
        self.assertEqual(data.loc[0, 'Accession'], 'B9LRY6')
        self.assertEqual(data.loc[0, 'Source'], 'Pfam')
        self.assertEqual(data.loc[0, 'Start'], 27)
        self.assertEqual(data.loc[0, 'End'], 514)


if __name__ == '__main__':
    logging.basicConfig(stream=sys.stderr)
    logging.getLogger("proteofav.config").setLevel(logging.DEBUG)
    suite = unittest.TestLoader().loadTestsFromTestCase(TestMSAS)
    unittest.TextTestRunner(verbosity=2).run(suite)
