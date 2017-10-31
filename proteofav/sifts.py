# -*- coding: utf-8 -*-

import os
import logging
import pandas as pd
from lxml import etree
from collections import OrderedDict

from proteofav.variants import fetch_uniprot_pdb_mapping
from proteofav.utils import (row_selector, InputFileHandler, Downloader,
                             GenericInputs, constrain_column_types, exclude_columns)
from proteofav.library import sifts_types

from proteofav.config import defaults

log = logging.getLogger('proteofav.config')

__all__ = ['parse_sifts_residues', 'select_sifts', 'sifts_best',
           'download_sifts', 'SIFTS']


def _parse_sifts_dbs_from_file(filename, excluded_cols=None):
    """
    Parses the residue fields of a SIFTS XML file.

    :param filename: path to the SIFTS file
    :param excluded_cols: option to exclude SIFTS dbSources
    :return: returns a nested dictionary
    """

    log.debug("Parsing SIFTS dbs from lines...")

    InputFileHandler(filename)

    # excluding columns
    if excluded_cols is None:
        excluded_cols = ("InterPro", "GO", "EC", "NCBI")

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


def _parse_sifts_regions_from_file(filename, excluded_cols=None):
    """
    Parses the residue fields of a SIFTS XML file.

    :param filename: path to the SIFTS file
    :param excluded_cols: option to exclude SIFTS dbSources
    :return: returns a nested dictionary
    """

    log.debug("Parsing SIFTS regions from lines...")

    InputFileHandler(filename)

    # excluding columns
    if excluded_cols is None:
        excluded_cols = ("InterPro", "GO", "EC", "NCBI")

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


def parse_sifts_residues(filename, add_regions=True, add_dbs=False,
                         excluded_cols=None):
    """
    Parses the residue fields of a SIFTS XML file.

    :param filename: path to the SIFTS file
    :param add_regions: boolean
    :param add_dbs: boolean
    :param excluded_cols: list of columns to be excluded
    :return: returns a pandas DataFrame
    """

    log.debug("Parsing SIFTS residues from lines...")

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

    # excluding columns
    if excluded_cols is None:
        excluded_cols = ("InterPro", "GO", "EC", "NCBI")

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

    # excluding columns
    table = exclude_columns(table, excluded=excluded_cols)

    # enforce some specific column types
    table = constrain_column_types(table, col_type_dict=sifts_types)

    if table.empty:
        raise ValueError('SIFTS file {} resulted in a empty Dataframe'
                         ''.format(filename))
    return table


def select_sifts(identifier, excluded_cols=None, add_regions=True, add_dbs=False,
                 overwrite=False, **kwargs):
    """
    Produce table ready from SIFTS XML file.

    :param identifier: PDB/mmCIF accession ID
    :param excluded_cols: option to exclude mmCIF columns
    :param add_regions: boolean
    :param add_dbs: boolean
    :param overwrite: boolean
    :return: returns a pandas DataFrame
    """

    filename = os.path.join(defaults.db_sifts, "{}.xml".format(identifier))

    download_sifts(identifier=identifier, filename=filename, overwrite=overwrite)

    table = parse_sifts_residues(filename=filename, excluded_cols=excluded_cols,
                                 add_regions=add_regions, add_dbs=add_dbs)

    table = filter_sifts(table, **kwargs)
    table = constrain_column_types(table, col_type_dict=sifts_types)
    return table


def filter_sifts(table, excluded_cols=None, chains=None,
                 chain_auth=None, res=None, uniprot=None, site=None):
    """
    Parses the residue fields of a SIFTS XML file.

    :param table: pandas DataFrame object
    :param excluded_cols: option to exclude SIFTS dbSources
    :param chains: (tuple) chain IDs or None
    :param chain_auth: (tuple) chain IDs or None
    :param res: (tuple) res IDs or None
    :param uniprot: (tuple) UniProt IDs or None
    :param site: (tuple) UniProt (positional) sites or None
    :return: returns a pandas DataFrame
    """

    # selections / filtering
    # excluding columns
    table = exclude_columns(table, excluded=excluded_cols)

    # excluding rows
    if chains is not None:
        table = row_selector(table, 'PDB_entityId', chains)
        log.debug("SIFTS table filtered by PDB_entityId...")

    if chain_auth is not None:
        table = row_selector(table, 'PDB_dbChainId', chain_auth)
        log.debug("SIFTS table filtered by PDB_dbChainId...")

    if res is not None:
        table = row_selector(table, 'PDB_dbResNum', res)
        log.debug("SIFTS table filtered by PDB_dbResNum...")

    if uniprot is not None:
        table = row_selector(table, 'UniProt_dbAccessionId', uniprot)
        log.debug("SIFTS table filtered by UniProt_dbAccessionId...")

    if site is not None:
        table = row_selector(table, 'UniProt_dbResNum', site)
        log.debug("SIFTS table filtered by UniProt_dbResNum...")
        log.debug("DSSP reset residue number...")

    if table.empty:
        raise ValueError("The filters resulted in an empty DataFrame...")
    return table


def download_sifts(identifier, filename, overwrite=False):
    """
     Downloads a SIFTS xml from the EBI FTP to the filesystem.

    :param identifier: (str) PDB accession ID
    :param filename: path to the SIFTS file
    :param overwrite: (boolean)
    :return: (side effects) output file path
    """

    url_root = defaults.sifts_fetch
    url_endpoint = "{}.xml.gz".format(identifier)
    url = url_root + url_endpoint
    Downloader(url=url, filename=filename,
               decompress=True, overwrite=overwrite)


def sifts_best(identifier, first=False):
    """
    Retrieves the best structures from the SIFTS endpoint in the PDBe API.

    :param identifier: UniProt accession ID
    :param first: gets the first entry
    :return: url content or url content in JSON data structure.
    """
    response = fetch_uniprot_pdb_mapping(identifier=identifier)
    if not response.ok or response is None:
        log.error('No SIFTS mapping found for {}'.format(identifier))
        return None
    return response.json() if not first else response.json()[identifier][0]


class SIFTS(GenericInputs):
    def read(self, filename=None, **kwargs):
        filename = self._get_filename(filename)
        self.table = parse_sifts_residues(filename=filename, **kwargs)
        return self.table

    def download(self, identifier=None, filename=None, **kwargs):
        identifier = self._get_identifier(identifier)
        filename = self._get_filename(filename)
        return download_sifts(identifier=identifier, filename=filename, **kwargs)

    def select(self, identifier=None, add_regions=True, add_dbs=False, **kwargs):
        identifier = self._get_identifier(identifier)
        self.table = select_sifts(identifier=identifier, add_regions=add_regions,
                                  add_dbs=add_dbs, **kwargs)
        return self.table


SIFTS = SIFTS()
