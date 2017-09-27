# -*- coding: utf-8

import logging
import pandas as pd
from lxml import etree
from io import StringIO
from string import ascii_uppercase
from collections import OrderedDict

from proteofav.utils import InputFileHandler

log = logging.getLogger('proteofav')


def parse_mmcif_atoms_from_file(filename):
    """
    Parse mmCIF ATOM and HETATM lines.

    :param filename: path to the mmCIF file
    :return: returns a pandas DataFrame
    """

    log.info("Parsing mmCIF atoms from lines...")

    # example lines with some problems
    """
    _atom_site.pdbx_PDB_model_num
    _atom_site.pdbe_label_seq_id
    _atom_site.orig_label_asym_id
    _atom_site.orig_auth_asym_id
    ATOM 1 N N . VAL A 1 1 ? -7.069 21.943 18.770 1.0 56.51 ? ? ? ? ? ? 118 VAL A N 1 1 A A
    ATOM 2 C CA . VAL A 1 1 ? -7.077 21.688 20.244 1.0 59.09 ? ? ? ? ? ? 118 VAL A CA 1 1 A A
    ATOM 3 C C . VAL A 1 1 ? -5.756 21.077 20.700 1.0 44.63 ? ? ? ? ? ? 118 VAL A C 1 1 A A
    ATOM 4 O O . VAL A 1 1 ? -5.346 20.029 20.204 1.0 59.84 ? ? ? ? ? ? 118 VAL A O 1 1 A A
    """

    InputFileHandler(filename)

    # parsing atom lines
    header = []
    lines = []
    with open(filename) as inlines:
        for line in inlines:
            if line.startswith("_atom_site."):
                header.append(line.split('.')[1].rstrip())
            elif line.startswith("ATOM") or "ATOM" in line[0:6]:
                lines.append(line)
            elif line.startswith("HETATM"):
                lines.append(line)
    lines = "".join(lines)

    all_str = {key: str for key in header}
    table = pd.read_table(StringIO(lines), delim_whitespace=True, low_memory=False,
                          names=header, compression=None, converters=all_str,
                          keep_default_na=False)
    return table


def parse_pdb_atoms_from_file(filename):
    """
    Parse PDB ATOM and HETATM lines. The ATOM lines are imported
    to the dictionary names used in the mmCIF format.

    :param filename: path to the PDB file
    :return: returns a pandas DataFrame
    """

    log.info("Parsing PDB atoms from lines...")

    # example lines
    """
    MODEL        1
    ATOM      0  N   SER A  -1     104.083  78.916  -1.349  1.00 61.47           N
    ATOM      1  N   VAL A 118      -7.069  21.943  18.770  1.00 56.51           N  
    ATOM      2  CA  VAL A 118      -7.077  21.688  20.244  1.00 59.09           C  
    ATOM      3  C   VAL A 118      -5.756  21.077  20.700  1.00 44.63           C  
    ATOM      4  O   VAL A 118      -5.346  20.029  20.204  1.00 59.84           O
    """

    InputFileHandler(filename)

    # parsing atom lines, converting it to mmcif-style headers
    lines = []
    modelnumb = '1'
    with open(filename) as inlines:
        for line in inlines:
            line = line.rstrip()
            line = line[0:78]
            if line.startswith("MODEL"):
                modelnumb = line.split()[1]
            elif line.startswith("ATOM"):
                lines.append(line + "%s" % modelnumb)
            elif line.startswith("HETATM"):
                lines.append(line + "%s" % modelnumb)
    lines = "\n".join(lines)

    header = ('group_PDB', 'id', 'label_atom_id', 'label_alt_id', 'label_comp_id',
              'label_asym_id', 'label_seq_id_full', 'label_seq_id',
              'pdbx_PDB_ins_code', 'Cartn_x', 'Cartn_y', 'Cartn_z',
              'occupancy', 'B_iso_or_equiv', 'type_symbol', 'auth_atom_id', 'auth_comp_id',
              'auth_asym_id', 'auth_seq_id_full', 'auth_seq_id', 'pdbx_PDB_model_num')

    # https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html
    widths = ((0, 6), (6, 11), (12, 16), (16, 17), (17, 20), (21, 22), (22, 27), (22, 26), (26, 27),
              (30, 38), (38, 46), (46, 54), (54, 60), (60, 66), (76, 78),  # (72, 76), ('seg_id')
              (12, 16), (17, 20), (21, 22), (22, 27), (22, 26), (78, 79))

    all_str = {key: str for key in header}
    table = pd.read_fwf(StringIO(lines), names=header, colspecs=widths,
                        compression=None, converters=all_str, keep_default_na=False)

    # fixes the 'pdbx_PDB_ins_code'
    table = _fix_pdb_ins_code(table)
    # fixes the 'label_alt_id
    table = _fix_label_alt_id(table)
    # fixes 'type_symbol' if missing
    table = _fix_type_symbol(table)

    return table


# TODO fix for speed
def _fix_pdb_ins_code(table):
    """
    Utility that fixes the 'pdbx_PDB_ins_code' column to match is expected
    in the mmcif format.

    :param table: pandas DataFrame object
    :return: returns a modified pandas DataFrame
    """

    ins_codes = []
    for i in table.index:
        value = table.loc[i, "pdbx_PDB_ins_code"]
        if value == '' or value == ' ' or value == '?':
            value = '?'
        ins_codes.append(value)
    table["pdbx_PDB_ins_code"] = ins_codes
    table['pdbx_PDB_ins_code'] = table['pdbx_PDB_ins_code'].fillna('?').astype(str)
    return table


# TODO fix for speed
def _fix_label_alt_id(table):
    """
    Utility that fixes the 'label_alt_id' column to match what is
    expected in the mmCIF format.

    :param table: pandas DataFrame object
    :return: returns a modified pandas DataFrame
    """

    alt_locs = []
    for i in table.index:
        value = table.loc[i, "label_alt_id"]
        if value == '' or value == ' ' or value == '?':
            value = '.'
        alt_locs.append(value)
    table["label_alt_id"] = alt_locs
    table['label_alt_id'] = table['label_alt_id'].fillna('.').astype(str)
    return table


# TODO fix for speed
def _fix_type_symbol(table):
    """
    Utility that fixes the 'type_symbol' column to match what is
    expected in the mmCIF format - when missing in the Structure.

    :param table: pandas DataFrame object
    :return: returns a modified pandas DataFrame
    """

    def get_type_symbol(table, key, key_fix):
        # this maybe a bit crude way of assigning this value
        if table[key] != " " and table[key] != "" and len(table[key]):
            return table[key]
        else:
            return ''.join([x for x in table[key_fix] if x in ascii_uppercase])[0]

    table.is_copy = False
    table['type_symbol'] = table.apply(get_type_symbol, axis=1,
                                       args=('type_symbol', 'label_atom_id'))
    return table


def parse_dssp_from_file(filename):
    """
    Parse lines of the DSSP file to get entries for every Residue
    in each CHAIN. The hierarchy is maintained. CHAIN->RESIDUE->[...].

    :param filename: path to the DSSP file
    :return: returns a pandas DataFrame
    """

    log.info("Parsing DSSP from lines...")

    # example lines with some problems
    """
      #  RESIDUE AA STRUCTURE BP1 BP2  ACC     N-H-->O    O-->H-N    N-H-->O    O-->H-N    TCO  KAPPA ALPHA  PHI   PSI    X-CA   Y-CA   Z-CA
        1    1 A M              0   0  127      0, 0.0   345,-0.1     0, 0.0     3,-0.1   0.000 360.0 360.0 360.0 162.0  -18.7   21.6  -55.4
        2    2 A R        +     0   0  117      1,-0.1    28,-0.4   343,-0.1     2,-0.3   0.455 360.0  81.5-136.8 -28.7  -17.0   22.3  -52.1

      381  394 A K              0   0  125     -2,-0.4   -21,-0.1   -21,-0.2    -2,-0.0  -0.421 360.0 360.0  64.1 360.0  -22.5   44.2  -25.4
      382        !*             0   0    0      0, 0.0     0, 0.0     0, 0.0     0, 0.0   0.000 360.0 360.0 360.0 360.0    0.0    0.0    0.0
      383    1 A M              0   0  127      0, 0.0   345,-0.1     0, 0.0     3,-0.1   0.000 360.0 360.0 360.0 162.0  -10.0   71.4  -55.4

    10278  103 H H  E     -XZ1023010269W  69     -9,-2.3    -9,-2.2    -2,-0.3     2,-1.0  -0.884  22.6-128.4-108.1 141.6  -97.0   28.7  112.2
    10279  104 H I  E     +XZ1022910268W   0    -50,-2.2   -50,-0.6    -2,-0.4   -11,-0.3  -0.801  30.6 175.4 -90.4  95.6  -98.5   32.0  111.3
    10280  105 H L  E     +     0   0   21    -13,-1.7   -55,-2.5    -2,-1.0     2,-0.3   0.812  62.6   4.9 -70.5 -35.5  -96.3   34.5  113.1

    # missing segment break
      145        !              0   0    0      0, 0.0     0, 0.0     0, 0.0     0, 0.0   0.000 360.0 360.0 360.0 360.0    0.0    0.0    0.0

    # chain break
      382        !*             0   0    0      0, 0.0     0, 0.0     0, 0.0     0, 0.0   0.000 360.0 360.0 360.0 360.0    0.0    0.0    0.0
    """

    InputFileHandler(filename)

    lines = []
    parse = False
    with open(filename) as inlines:
        for line in inlines:
            line = line.rstrip()
            if parse:
                lines.append(line)
            if line.startswith("  #"):
                parse = True
    lines = "\n".join(lines)

    # column width descriptors
    header = ("LINE", "RES", "CHAIN", "AA", "SS", "STRUCTURE",
              "BP1", "BP2", "BP2_CHAIN", "ACC",
              "NH_O_1", "NH_O_1_nrg", "O_HN_1", "O_HN_1_nrg",
              "NH_O_2", "NH_O_2_nrg", "O_HN_2", "O_HN_2_nrg",
              "TCO", "KAPPA", "ALPHA", "PHI", "PSI",
              "X-CA", "Y-CA", "Z-CA")

    widths = ((0, 5), (5, 11), (11, 12), (12, 15), (16, 17), (17, 25),
              (25, 29), (29, 33), (33, 34), (34, 38),
              (38, 45), (46, 50), (50, 56), (57, 61),
              (61, 67), (68, 72), (72, 78), (79, 84),
              (85, 91), (91, 97), (97, 103), (103, 109), (109, 115),
              (115, 123), (123, 130), (130, 137))

    all_str = {key: str for key in header}
    table = pd.read_fwf(StringIO(lines), names=header, colspecs=widths,
                        compression=None, converters=all_str, keep_default_na=False)
    return table


def _parse_sifts_dbs_from_file(filename, excluded_cols=()):
    """
    Parses the residue fields of a SIFTS XML file.

    :param filename: path to the SIFTS file
    :param excluded_cols: option to exclude SIFTS dbSources
    :return: returns a nested dictionary
    """

    log.info("Parsing SIFTS dbs from lines...")

    InputFileHandler(filename)

    # parsing sifts segments
    try:
        tree = etree.parse(filename)
    except etree.XMLSyntaxError:
        raise IOError("{} not available or could not be read...".format(filename))
    root = tree.getroot()
    namespace = root.nsmap[None]
    namespace_map = {'ns': namespace}

    # DB entries and versions
    sifts_object = OrderedDict()
    for db_list in root.iterfind('.//ns:listDB', namespaces=namespace_map):
        for db in db_list:
            source = ''
            db_entries = OrderedDict()
            for k, v in db.attrib.items():
                db_entries[k] = v
                if k == 'dbSource' and v not in excluded_cols:
                    source = v
            if source not in excluded_cols and source != '':
                sifts_object[source] = db_entries

    return sifts_object


def _parse_sifts_regions_from_file(filename, excluded_cols=()):
    """
    Parses the residue fields of a SIFTS XML file.

    :param filename: path to the SIFTS file
    :param excluded_cols: option to exclude SIFTS dbSources
    :return: returns a nested dictionary
    """

    log.info("Parsing SIFTS regions from lines...")

    InputFileHandler(filename)

    # parsing sifts segments
    try:
        tree = etree.parse(filename)
    except etree.XMLSyntaxError:
        raise IOError("{} not available or could not be read...".format(filename))
    root = tree.getroot()
    namespace = root.nsmap[None]
    namespace_map = {'ns': namespace}
    db_reference = "{{{}}}db".format(namespace)

    sifts_object = OrderedDict()

    # Entities
    for entity_list in root.iterfind('.//ns:entity[@type="protein"]',
                                     namespaces=namespace_map):
        entity_id = entity_list.attrib['entityId']
        regions_full = OrderedDict()
        # Entities : Segments
        for segment in entity_list:

            # parse the regions found for this segment
            # Entities : Segments : Regions
            for region_list in segment.iterfind('.//ns:listMapRegion',
                                                namespaces=namespace_map):
                for region in region_list:
                    # get region annotations
                    # parse extra annotations for each region
                    for annotation in region:
                        for k, v in annotation.attrib.items():
                            # db entries
                            if annotation.tag == db_reference:
                                if k == 'dbSource' and v not in excluded_cols:
                                    source = v
                                else:
                                    continue
                                if source not in regions_full:
                                    counter = 1
                                    region_object = OrderedDict()
                                else:
                                    region_object = regions_full[source]
                                    keys = regions_full[source]
                                    last_counter = [int(k) for k in keys][-1]
                                    counter = last_counter + 1
                                try:
                                    coord = annotation.attrib['dbCoordSys']
                                except KeyError:
                                    coord = '-'
                                region_annotation = OrderedDict()
                                region_annotation['dbAccessionId'] = annotation.attrib['dbAccessionId']
                                region_annotation['start'] = int(region.attrib['start'])
                                region_annotation['end'] = int(region.attrib['end'])
                                region_annotation['dbCoordSys'] = coord
                                region_object["{}".format(counter)] = region_annotation
                                regions_full[source] = region_object

        sifts_object[entity_id] = regions_full

    return sifts_object


def parse_sifts_residues_from_file(filename, add_regions=True, add_dbs=False,
                                   excluded_cols=()):
    """
    Parses the residue fields of a SIFTS XML file.

    :param filename: path to the SIFTS file
    :param add_regions: boolean
    :param add_dbs: boolean
    :return: returns a pandas DataFrame
    """

    log.info("Parsing SIFTS residues from lines...")

    # example lines with some problems
    """
    <?xml version='1.0' encoding='UTF-8' standalone='yes'?>
    <entry xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xmlns:align="http://www.ebi.ac.uk/pdbe/docs/sifts/alignment.xsd" xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:data="http://www.ebi.ac.uk/pdbe/docs/sifts/dataTypes.xsd" xmlns="http://www.ebi.ac.uk/pdbe/docs/sifts/eFamily.xsd" dbSource="PDBe" dbVersion="1.3" dbCoordSys="PDBe" dbAccessionId="2pah" dbEntryVersion="2011-07-13" date="2014-07-21" xsi:schemaLocation="http://www.ebi.ac.uk/pdbe/docs/sifts/eFamily.xsd http://www.ebi.ac.uk/pdbe/docs/sifts/eFamily.xsd">
      <rdf:RDF>
        <rdf:Description rdf:about="self">
          <dc:rights rdf:resource="http://pdbe.org/sifts">
            Copyright notice: (c) 2004-2013, EMBL-EBI, PDBe-UniProt
            Jose M. Dana, Paul Gane, Jie Luo, Glen van Ginkel, Claire O'Donovan, Maria J. Martin, Sameer Velankar.
            The information included is supplied as-is, under the terms and conditions of the licence agreement.
        </dc:rights>
        </rdf:Description>
      </rdf:RDF>
      <listDB>
        <db dbSource="Pfam" dbCoordSys="UniProt" dbVersion="27.0"/>
        <db dbSource="InterPro" dbCoordSys="UniProt" dbVersion="48.0"/>
        <db dbSource="CATH" dbCoordSys="PDBresnum" dbVersion="3.5.0"/>
        <db dbSource="EC" dbCoordSys="UniProt" dbVersion="30.14"/>
        <db dbSource="UniProt" dbCoordSys="UniProt" dbVersion="2014.08"/>
        <db dbSource="SCOP" dbCoordSys="PDBresnum" dbVersion="1.75"/>
        <db dbSource="GO" dbCoordSys="UniProt" dbVersion="20140708"/>
        <db dbSource="PDB" dbCoordSys="PDBresnum" dbVersion="30.14"/>
      </listDB>
      <entity type="protein" entityId="A">
        <segment segId="2pah_A_1_335" start="1" end="335">
          <listResidue>
            <residue dbSource="PDBe" dbCoordSys="PDBe" dbResNum="1" dbResName="VAL">
              <crossRefDb dbSource="PDB" dbCoordSys="PDBresnum" dbAccessionId="2pah" dbResNum="118" dbResName="VAL" dbChainId="A"/>
              <crossRefDb dbSource="UniProt" dbCoordSys="UniProt" dbAccessionId="P00439" dbResNum="118" dbResName="V"/>
              <crossRefDb dbSource="CATH" dbCoordSys="PDBresnum" dbAccessionId="1.10.800.10" dbResNum="118" dbResName="VAL" dbChainId="A"/>
              <crossRefDb dbSource="SCOP" dbCoordSys="PDBresnum" dbAccessionId="42581" dbResNum="118" dbResName="VAL" dbChainId="A"/>
              <crossRefDb dbSource="NCBI" dbCoordSys="UniProt" dbAccessionId="9606" dbResNum="118" dbResName="V"/>
              <crossRefDb dbSource="InterPro" dbCoordSys="UniProt" dbAccessionId="IPR001273" dbResNum="118" dbResName="V" dbEvidence="G3DSA:1.10.800.10"/>
              <crossRefDb dbSource="InterPro" dbCoordSys="UniProt" dbAccessionId="IPR019774" dbResNum="118" dbResName="V" dbEvidence="PS51410"/>
              <crossRefDb dbSource="InterPro" dbCoordSys="UniProt" dbAccessionId="IPR019774" dbResNum="118" dbResName="V" dbEvidence="SSF56534"/>
              <crossRefDb dbSource="InterPro" dbCoordSys="UniProt" dbAccessionId="IPR019774" dbResNum="118" dbResName="V" dbEvidence="SSF56534"/>
              <crossRefDb dbSource="InterPro" dbCoordSys="UniProt" dbAccessionId="IPR005961" dbResNum="118" dbResName="V" dbEvidence="TIGR01268"/>
              <crossRefDb dbSource="InterPro" dbCoordSys="UniProt" dbAccessionId="IPR001273" dbResNum="118" dbResName="V" dbEvidence="PTHR11473"/>
              <residueDetail dbSource="PDBe" property="codeSecondaryStructure">T</residueDetail>
              <residueDetail dbSource="PDBe" property="nameSecondaryStructure">loop</residueDetail>
            </residue>
            (...)
    """

    InputFileHandler(filename)

    # parse regions first
    if add_regions:
        regions = _parse_sifts_regions_from_file(filename=filename,
                                                 excluded_cols=excluded_cols)
    if add_dbs:
        dbs = _parse_sifts_dbs_from_file(filename=filename,
                                         excluded_cols=excluded_cols)

    # parsing sifts segments
    try:
        tree = etree.parse(filename)
    except etree.XMLSyntaxError:
        raise IOError("{} not available or could not be read...".format(filename))
    root = tree.getroot()
    namespace = root.nsmap[None]
    namespace_map = {'ns': namespace}
    cross_reference = "{{{}}}crossRefDb".format(namespace)
    residue_detail = "{{{}}}residueDetail".format(namespace)

    rows = []
    # Entities
    for entity_list in root.iterfind('.//ns:entity[@type="protein"]',
                                     namespaces=namespace_map):

        entity_id = entity_list.attrib['entityId']
        # Entities : Segments
        for segment in entity_list:

            # Entities : Segments : Residues
            for list_residue in segment.iterfind('.//ns:listResidue',
                                                 namespaces=namespace_map):
                for residue in list_residue:
                    # get residue annotations
                    residue_annotation = OrderedDict()
                    # key, value pairs
                    for k, v in residue.attrib.items():
                        # skipping dbSource
                        if k == 'dbResNum':
                            # adding to the dictionary
                            # residue_annotation[k] = v
                            resnum = int(v)
                        else:
                            continue

                    # parse extra annotations for each residue
                    for annotation in residue:
                        for k, v in annotation.attrib.items():
                            # crossRefDb entries
                            if annotation.tag == cross_reference:
                                if annotation.attrib["dbSource"] not in excluded_cols:
                                    # skipping various fields
                                    if k == 'dbSource' or k == 'dbCoordSys':
                                        continue
                                    elif (annotation.attrib["dbSource"] != "PDB" and
                                                  annotation.attrib["dbSource"] != "UniProt"):
                                        if k == 'dbResName' or k == 'dbResNum' or k == 'dbChainId':
                                            continue
                                    # elif annotation.attrib["dbSource"] == "PDB" and k == "dbAccessionId":
                                    #     continue

                                    # adding a new column with the regionId from the 'regions'
                                    if add_regions:
                                        if k == "dbAccessionId":
                                            source = annotation.attrib["dbSource"]
                                            if source in regions[entity_id]:
                                                keys = regions[entity_id][source]
                                                for key in [k for k in keys]:
                                                    entry = regions[entity_id][source][key]
                                                    if v == entry["dbAccessionId"]:
                                                        start = int(entry["start"])
                                                        end = int(entry["end"])
                                                        if resnum in range(start, end + 1, 1):
                                                            nk = "{}_regionId".format(source)
                                                            residue_annotation[nk] = key
                                                            nk = "{}_regionStart".format(source)
                                                            residue_annotation[nk] = start
                                                            nk = "{}_regionEnd".format(source)
                                                            residue_annotation[nk] = end
                                                            nk = "{}_regionResNum".format(source)
                                                            residue_annotation[nk] = resnum

                                    if add_dbs:
                                        source = annotation.attrib["dbSource"]
                                        if source in dbs:
                                            nk = "{}_dbVersion".format(source)
                                            residue_annotation[nk] = dbs[source]['dbVersion']

                                    # renaming all keys with dbSource prefix
                                    k = "{}_{}".format(
                                        annotation.attrib["dbSource"], k)
                                    if k == "{}_1".format(annotation.attrib["dbSource"]):
                                        continue

                            if annotation.tag == residue_detail:
                                if 'PDB' not in excluded_cols:
                                    k = "PDB_{}".format(annotation.attrib["property"])
                                    # value is the text field in the XML
                                    v = annotation.text

                            # adding to the dictionary
                            if "_" in k:
                                try:
                                    if v in residue_annotation[k]:
                                        continue
                                    if k != "PDB_Annotation":
                                        residue_annotation[k].append(v)
                                    else:
                                        residue_annotation[k] = v
                                except KeyError:
                                    residue_annotation[k] = v
                                except AttributeError:
                                    residue_annotation[k] = [residue_annotation[k]]
                                    residue_annotation[k].append(v)
                                except TypeError:
                                    # bool column for annotation
                                    residue_annotation[k] = v
                        if 'PDB' not in excluded_cols and "PDB_Annotation" not in residue_annotation:
                            residue_annotation["PDB_Annotation"] = "Observed"

                        # adding the entity_id
                        if 'PDB_entityId' not in residue_annotation:
                            residue_annotation['PDB_entityId'] = entity_id

                    rows.append(residue_annotation)

    table = pd.DataFrame(rows)

    # # enforce some specific column types
    # table = constrain_column_types(table, sifts_types)

    for c in list(table):
        if '_regionId' in c:
            table[c] = table[c].fillna('-').astype(str)
        elif '_regionStart' in c or '_regionEnd' in c:
            table[c] = table[c].fillna(0).astype(int)

    return table
