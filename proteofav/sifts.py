# -*- coding: utf-8 -*-

from os import path
import logging
import pandas as pd
from lxml import etree
from requests import HTTPError

from proteofav.utils import fetch_files, row_selector, fetch_from_url_or_retry

from proteofav.config import defaults

log = logging.getLogger('proteofav.config')

__all__ = ['_sifts_residues_regions', 'select_sifts', 'sifts_best']


def _sifts_residues_regions(filename, cols=None,
                            sources=('CATH', 'SCOP', 'Pfam', 'InterPro')):
    """
    Parses the residue fields of a SIFTS XML file to a pandas DataFrame.

    :param filename: input SIFTS xml file
    :param cols: option to select a columns (post-parsing)
    :param sources: option to select SIFTS dbSources
    :return: pandas table dataframe
    """

    # parsing sifts segments
    tree = etree.parse(filename)
    root = tree.getroot()
    namespace = root.nsmap[None]
    namespace_map = {'ns': namespace}
    cross_reference = "{{{}}}crossRefDb".format(namespace)
    residue_detail = "{{{}}}residueDetail".format(namespace)
    db_reference = "{{{}}}db".format(namespace)
    # db_detail = "{{{}}}dbDetail".format(namespace)
    rows = []

    sources += ('PDB', 'UniProt')

    for segment_list in root.iterfind('.//ns:entity[@type="protein"]',
                                      namespaces=namespace_map):

        entity_id = segment_list.attrib['entityId']
        for segment in segment_list:
            # 1st parse the regions found for this segment
            regions_full = {}
            for region_list in segment.iterfind('.//ns:listMapRegion',
                                                namespaces=namespace_map):
                region_source = {}
                for region in region_list:
                    # get region annotations
                    # parse extra annotations for each region
                    for annotation in region:
                        for k, v in annotation.attrib.items():
                            # db entries
                            if annotation.tag == db_reference:
                                if k == 'dbSource' and v in sources:
                                    source = v
                                else:
                                    continue

                                if source not in region_source:
                                    region_source[source] = []
                                try:
                                    coord = annotation.attrib['dbCoordSys']
                                except KeyError:
                                    coord = '-'
                                region_source[source].append([annotation.attrib['dbAccessionId'],
                                                              region.attrib['start'],
                                                              region.attrib['end'],
                                                              coord])
                regions_full[entity_id] = region_source

            # 2nd parse each residue and
            for list_residue in segment.iterfind('.//ns:listResidue',
                                                 namespaces=namespace_map):
                for residue in list_residue:
                    # get residue annotations
                    residue_annotation = {}
                    # key, value pairs
                    for k, v in residue.attrib.items():
                        # skipping dbSource
                        if k == 'dbSource' or k == 'dbCoordSys' or k == 'dbResName':
                            continue

                        k = "PDB_index"
                        # adding to the dictionary
                        residue_annotation[k] = v
                        resnum = int(v)

                    # parse extra annotations for each residue
                    for annotation in residue:
                        for k, v in annotation.attrib.items():
                            # crossRefDb entries
                            if annotation.tag == cross_reference:
                                if annotation.attrib["dbSource"] in sources:
                                    # skipping dbSource
                                    if k == 'dbSource' or k == 'dbCoordSys':
                                        continue
                                    if (annotation.attrib["dbSource"] != "PDB" and
                                                annotation.attrib["dbSource"] != "UniProt"):
                                        if k == 'dbResName' or k == 'dbResNum' or k == 'dbChainId':
                                            continue
                                    if annotation.attrib["dbSource"] == "PDB" and k == "dbAccessionId":
                                        continue

                                    # adding a new column with the regionId from the 'regions'
                                    if k == "dbAccessionId":
                                        if annotation.attrib["dbSource"] in regions_full[entity_id]:
                                            for c, entry in enumerate(
                                                    regions_full[entity_id][annotation.attrib["dbSource"]]):
                                                if v == entry[0]:
                                                    start = int(entry[1])
                                                    end = int(entry[2])
                                                    if resnum in range(start, end + 1, 1):
                                                        nk = "{}_regionId".format(annotation.attrib["dbSource"])
                                                        nv = str(c + 1)
                                                        residue_annotation[nk] = nv

                                    # renaming all keys with dbSource prefix
                                    k = "{}_{}".format(
                                        annotation.attrib["dbSource"], k)

                            if annotation.tag == residue_detail:
                                k = "PDB_{}".format(annotation.attrib["property"])
                                # value is the text field in the XML
                                v = annotation.text

                            # adding to the dictionary
                            if "_" in k:
                                try:
                                    if v in residue_annotation[k]:
                                        continue
                                    residue_annotation[k].append(v)
                                except KeyError:
                                    residue_annotation[k] = v
                                except AttributeError:
                                    residue_annotation[k] = [residue_annotation[k]]
                                    residue_annotation[k].append(v)
                                except TypeError:
                                    # bool column for annotation
                                    residue_annotation[k] = v

                    rows.append(residue_annotation)
    if cols:
        data = pd.DataFrame(rows, columns=cols)
    else:
        data = pd.DataFrame(rows)
    return data


def select_sifts(pdb_id, chains=None):
    """
    Produce table ready from SIFTS XML file.

    :param pdb_id: PDB identifier
    :param chains: Protein structure chain
    :return: table read to be merged
    """
    sifts_path = path.join(defaults.db_sifts, pdb_id + '.xml')

    try:
        sift_table = _sifts_residues_regions(sifts_path,
                                             sources=('CATH', 'SCOP', 'Pfam'))
    except IOError:
        sifts_path = fetch_files(pdb_id, sources='sifts',
                                 directory=defaults.db_sifts)[0]
        sift_table = _sifts_residues_regions(sifts_path)
        # standardise column types
    for col in sift_table:
        #  bool columns
        if col.startswith('is'):
            sift_table[col].fillna(False)
    if chains is None:
        return sift_table
    else:
        return row_selector(sift_table, 'PDB_dbChainId', chains)


def sifts_best(uniprot_id, first=False):
    """
    Retrieves the best structures from the SIFTS endpoint in the PDBe api.

    :param uniprot_id: Uniprot ID
    :param first: gets the first entry
    :return: url content or url content in json data structure.
    """
    sifts_endpoint = "mappings/best_structures/"
    url = defaults.api_pdbe + sifts_endpoint + str(uniprot_id)
    try:
        response = fetch_from_url_or_retry(url, json=True).json()
    except HTTPError as e:
        if e.response.status_code == 404:
            logging.error('No SIFTS mapping found for {}'.format(uniprot_id))
            return None
        else:
            raise
    return response if not first else response[uniprot_id][0]
