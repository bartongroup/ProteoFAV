# -*- coding: utf-8 -*-

import re
import copy
import logging
import pandas as pd
from Bio import AlignIO


from proteofav.variants import (fetch_uniprot_id_from_name, fetch_pdb_uniprot_mapping)
from proteofav.utils import (constrain_column_types, exclude_columns,
                             InputFileHandler)

log = logging.getLogger('proteofav.config')

__all__ = ['read_alignments', 'read_msas',
           'parse_sequence_info_from_description',
           'parse_uniprot_fasta_seq_description',
           'parse_pfam_sth_seq_description',
           'parse_cath_fasta_seq_description',
           'parse_cath_sth_seq_description',
           'parse_generic_seq_description']

SEQ_FORMAT_VALID = ('clustal', 'emboss', 'nexus', 'fasta', 'phylip', 'stockholm')


def read_alignments(filename, seq_format=None):
    """
    Reads an input multiple sequence alignment (MSA).
    Seq Formats recognised by Biopython: see `SEQ_FORMAT_VALID` above.

    :param filename: Input MSA (read by Biopython)
    :param seq_format: (str) or None. Valid formats from Biopython.
    :return: returns the Biopython alignment object
    """

    InputFileHandler(filename)

    if seq_format is not None:
        seq_format = seq_format.lower()
        if seq_format not in SEQ_FORMAT_VALID:
            seq_format = None

    if seq_format is None:
        # guess the format from the file extension
        if filename.endswith('.fasta') or filename.endswith('.fa'):
            seq_format = 'fasta'
        elif filename.endswith('.sto') or filename.endswith('.sth'):
            seq_format = 'stockholm'
        elif filename.endswith('.aln') or filename.endswith('.clw'):
            seq_format = 'clustal'

    if seq_format is None:
        raise ValueError("Alignment format unrecognised...")

    alignment = AlignIO.read(filename, seq_format)
    return alignment, seq_format


def read_msas(filename, excluded_cols=(), seq_format=None, get_uniprot_id=True):
    """
    Reads a Pfam/CATH MSA and returns a pandas table with a
    collection of protein IDs and sequences

    :param filename: path to the MSA file
    :param excluded_cols: option to exclude some columns
    :param seq_format: (str) or None. Valid formats from Biopython.
    :param get_uniprot_id: (boolean)
    :return: returns a pandas DataFrame
    """

    rows = []
    alignment, seq_format = read_alignments(filename, seq_format)
    for record in alignment:
        seq = str(record.seq)
        if seq_format == 'fasta':
            desc = str(record.description)
        else:
            desc = str(record.id)
        # get cross-reference information for each entry
        entry = dict()
        entry['Sequence'] = seq
        entry['Seq_Format'] = seq_format
        # parse the sequence description information
        entry = parse_sequence_info_from_description(desc, entry, seq_format,
                                                     get_uniprot_id)
        rows.append(entry)

    table = pd.DataFrame(rows)

    # excluding columns
    table = exclude_columns(table, excluded_cols)

    # enforce some specific column types
    msa_types = {key: str for key in list(table) if key != 'Start' and key != 'End'}
    table = constrain_column_types(table, msa_types)

    if table.empty:
        raise ValueError('{} resulted in an empty DataFrame...'.format(filename))
    return table


def parse_sequence_info_from_description(description, entry, seq_format=None,
                                         get_uniprot_id=True):
    """
    Parses the Biopython alignment sequence description and tries to guess
    the content. (only works for known formats).

    :param description: Biopython's 'description' field
    :param entry: Mutable dictionary
    :param seq_format: (str) or None
    :param get_uniprot_id: (boolean) Tries to parse UniProt IDs
    :return: (side-effects) updated Dictionary
    """

    prev_entry = copy.deepcopy(entry)

    if seq_format == 'fasta':
        # trying the UniProt fasta seq description
        parse_uniprot_fasta_seq_description(description, entry)
        if entry != prev_entry:
            return entry

        # trying the CATH fasta seq description
        parse_cath_fasta_seq_description(description, entry,
                                         get_uniprot_id=get_uniprot_id)
        if entry != prev_entry:
            return entry

    elif seq_format == 'stockholm':
        # trying the Pfam Stockholm seq description
        parse_pfam_sth_seq_description(description, entry,
                                       get_uniprot_id=get_uniprot_id)
        if entry != prev_entry:
            return entry

        # trying the CATH Stockholm seq description
        parse_cath_sth_seq_description(description, entry,
                                       get_uniprot_id=get_uniprot_id)
        if entry != prev_entry:
            return entry

    else:
        # trying a generic sequence description
        parse_generic_seq_description(description, entry,
                                      get_uniprot_id=get_uniprot_id)
        if entry != prev_entry:
            return entry

    log.debug("Nothing parsed from the MSA sequence description...")
    return entry


def parse_uniprot_fasta_seq_description(description, entry):
    """
    Pattern: <source>|<Accession_ID>|<Accession_Name> ++
    Example: sp|P00439|PH4H_HUMAN Phenylalanine-4-hydroxylase (...)
        OS=Homo sapiens GN=PAH PE=1 SV=1

    :param description: Biopython's 'description' field
    :param entry: Mutable dictionary
    :return: (updated) Dictionary
    """

    # trying the UniProt fasta seq description
    pattern = re.compile("([a-zA-Z])+\|([A-Z0-9])+\|([A-Z0-9])+_([A-Z0-9])+")
    match = re.search(pattern, description, flags=0)
    if match:
        match = match.group()
        # pattern matching
        pattern = re.compile("([a-zA-Z])+\|")
        source = re.search(pattern, match, flags=0)
        if source:
            source = source.group().rstrip('|')
            entry['Collection'] = source

        pattern = re.compile("\|([a-zA-Z0-9])+\|")
        identifier = re.search(pattern, match, flags=0)
        if identifier:
            identifier = identifier.group().lstrip('|').rstrip('|')
            entry['Accession'] = identifier

        pattern = re.compile("\|([A-Z0-9])+_([A-Z0-9])+")
        name = re.search(pattern, match, flags=0)
        if name:
            name = name.group().lstrip('|')
            entry['Name'] = name

        entry['Source'] = 'UniProt'

        # remaining description
        if description != match:
            description = description.replace(match, "")
            entry['Description'] = description.strip()
    return entry


def parse_pfam_sth_seq_description(description, entry, get_uniprot_id=True):
    """
    Pattern: <Accession_Name>/<Start>-<End>
    Example: C7P4T5_HALMD/44-372

    :param description: Biopython's 'description' field
    :param entry: Mutable dictionary
    :param get_uniprot_id: (boolean)
    :return: (updated) Dictionary
    """

    pattern = re.compile("([A-Z0-9])+_([A-Z0-9])+\/[\-]?([0-9])+-[\-]?([0-9])+")
    match = re.search(pattern, description, flags=0)
    if match:
        match = match.group()
        entry = parse_generic_seq_description(match, entry,
                                              get_uniprot_id=get_uniprot_id)

        entry['Source'] = 'Pfam'

        # remaining description
        if description != match:
            description = description.replace(match, "")
            entry['Description'] = description.strip()
    return entry


def parse_cath_sth_seq_description(description, entry, get_uniprot_id=True):
    """
    Pattern: <Accession_Name>/<Start>-<End>
    Example: C7P4T5_HALMD/44-372

    :param description: Biopython's 'description' field
    :param entry: Mutable dictionary
    :param get_uniprot_id: (boolean)
    :return: (updated) Dictionary
    """

    pattern = re.compile("([a-zA-Z0-9])+\/[\-]?([0-9])+-[\-]?([0-9])+")
    match = re.search(pattern, description, flags=0)
    if match:
        match = match.group()
        entry = parse_generic_seq_description(match, entry,
                                              get_uniprot_id=get_uniprot_id)

    entry['Source'] = 'CATH'
    return entry


def parse_cath_fasta_seq_description(description, entry, get_uniprot_id=True):
    """
    Pattern: cath|<version>|<Accession_ID>/<Start>-<End>
    Example:
    # sequence from CATH alignments - with structure domain
    # cath|4.1.0|1rwcA01/4-372 (...)
        CATH_S35=1.50.10.100.1;GENE=P84141_ARTAU;GO=GO:0005576,GO:0005975;
        MDA=1.50.10.100;ORG=Arthrobacter_aurescens;
        TAXON=43663;UNIPROT=P84141

    # sequence from CATH alignments - without structure domain
    # biomap|4.1.0|b7f2808bde19a485bdd0a90545c40d99/29-337 (...)
        EC=4.2.2.5;GENE=CSLA_PEDHD;GO=GO:0005576,GO:0005975;
        MDA=1.50.10.100_2.70.98.10_2.60.220.10;
        ORG=Pedobacter_heparinus;SWISSPROT=1;
        TAXON=485917;UNIPROT=Q59288

    :param description: Biopython's 'description' field
    :param entry: Mutable dictionary
    :param get_uniprot_id: (boolean)
    :return: (updated) Dictionary
    """
    # trying the CATH fasta seq description
    pattern = re.compile("([a-zA-Z])+\|([0-9])(.|-)([0-9])(.|-)([0-9])\|"
                         "([a-zA-Z0-9])+\/[\-]?([0-9])+-[\-]?([0-9])+")
    match = re.search(pattern, description, flags=0)
    if match:
        match = match.group()
        # pattern matching
        pattern = re.compile("([a-zA-Z])+\|")
        source = re.search(pattern, match, flags=0)
        if source:
            source = source.group().rstrip('|')
            entry['Collection'] = source

        pattern = re.compile("\|([0-9])(.|-)([0-9])(.|-)([0-9])\|")
        version = re.search(pattern, match, flags=0)
        if version:
            version = version.group().lstrip('|').rstrip('|')
            entry['Version'] = version

        pattern = re.compile(
            "([A-Z0-9])+[\_]?([a-zA-Z0-9])+\/[\-]?([0-9])+-[\-]?([0-9])+")
        # same as the generic seq description
        nmatch = re.search(pattern, match, flags=0)
        if nmatch:
            entry = parse_generic_seq_description(nmatch.group(), entry,
                                                  get_uniprot_id=get_uniprot_id)

        entry['Source'] = 'CATH'

        # remaining description
        if description != match:
            description = description.replace(match, "")
            entry['Description'] = description.strip()
    return entry


def parse_generic_seq_description(description, entry, get_uniprot_id=True):
    """
    Pattern: <Accession_ID>/<Start>-<End> or <Accession_Name>/<Start>-<End>
    Example: P00439/24-145

    :param description: Biopython's 'description' field
    :param entry: Mutable dictionary
    :param get_uniprot_id: (boolean)
    :return: (updated) Dictionary
    """

    # trying a generic sequence description
    pattern = re.compile(
        "([A-Z0-9])+[\_]?([a-zA-Z0-9])+\/[\-]?([0-9])+-[\-]?([0-9])+")
    match = re.search(pattern, description, flags=0)
    if match:
        match = match.group()
        # pattern matching
        pattern = re.compile("([A-Z0-9])+_([a-zA-Z0-9])+\/")
        name = re.search(pattern, match, flags=0)
        if name:
            name = name.group().rstrip('/')
            entry['Name'] = name

        pattern = re.compile("([a-zA-Z0-9])+\/")
        identifier = re.search(pattern, match, flags=0)
        if identifier:
            identifier = identifier.group().rstrip('/')
            entry['Accession'] = identifier

        pattern = re.compile("^([a-zA-Z0-9]){5}[0-9][0-9]\/")
        domain = re.search(pattern, match, flags=0)
        if domain:
            domain = domain.group()
            entry['pdb_id'] = domain[0:4]
            entry['chain_id'] = domain[4:5]
            entry['domain_id'] = domain.rstrip('/')
            entry['Domain'] = domain.rstrip('/')
            entry['Collection'] = "cath"

        pattern = re.compile("\/[\-]?([0-9])+")
        start = re.search(pattern, match, flags=0)
        if start:
            start = int(start.group().lstrip('/'))
            entry['Start'] = start

        pattern = re.compile("-[\-]?([0-9])+")
        end = re.search(pattern, match, flags=0)
        if end:
            end = int((end.group())[1:])
            entry['End'] = end

        # UniProt ID missing
        if name and get_uniprot_id:
            r = fetch_uniprot_id_from_name(name)
            try:
                identifier = r.json()[0]["id"]
            except Exception:
                # JSONDecodeError
                # new - to match the update in the UniProt API endpoint
                identifier = str(r.content, encoding='utf-8').strip()
            entry['Accession'] = identifier
        elif domain and get_uniprot_id:
            r = fetch_pdb_uniprot_mapping(domain[0:4])
            try:
                identifier = [k for k in r.json()[domain[0:4]]['UniProt']][0]
                entry['Accession'] = identifier
            except Exception:
                pass

        entry['Source'] = 'GenericParser'

        # remaining description
        if description != match:
            description = description.replace(match, "")
            entry['Description'] = description.strip()
    return entry
