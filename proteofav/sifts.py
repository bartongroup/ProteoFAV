# -*- coding: utf-8

import logging

from proteofav.parsers import parse_sifts_residues_from_file
from proteofav.downloaders import Downloader

from proteofav.utils import InputFileHandler
from proteofav.utils import OutputFileHandler
from proteofav.utils import GenericInput
from proteofav.utils import exclude_columns
from proteofav.utils import constrain_column_types
from proteofav.utils import row_selector
from proteofav.library import sifts_types

from proteofav.config import defaults as config

log = logging.getLogger('proteofav')


class SIFTS(GenericInput):
    def read(self, filename=None, excluded_cols=None,
             chains=None, chain_auth=None, res=None, uniprot=None, site=None,
             add_regions=True, add_dbs=False):
        """
        Parses the residue fields of a SIFTS XML file.

        :param filename: path to the SIFTS file
        :param excluded_cols: option to exclude SIFTS dbSources
        :param chains: (tuple) chain IDs or None
        :param chain_auth: (tuple) chain IDs or None
        :param res: (tuple) res IDs or None
        :param uniprot: (tuple) UniProt IDs or None
        :param site: (tuple) UniProt (positional) sites or None
        :param add_regions: boolean
        :param add_dbs: boolean
        :return: returns a pandas DataFrame
        """

        filename = self._get_filename(filename)
        InputFileHandler(filename)

        if excluded_cols is None:
            excluded_cols = ("InterPro", "GO", "EC", "NCBI")

        table = parse_sifts_residues_from_file(filename, excluded_cols=excluded_cols,
                                               add_regions=add_regions, add_dbs=add_dbs)

        # excluding columns
        table = exclude_columns(table, excluded=excluded_cols)

        # enforce some specific column types
        table = constrain_column_types(table, sifts_types)

        # excluding rows
        if chains is not None:
            table = row_selector(table, 'PDB_entityId', chains)
            log.info("SIFTS table filtered by PDB_entityId...")

        if chain_auth is not None:
            table = row_selector(table, 'PDB_dbChainId', chain_auth)
            log.info("SIFTS table filtered by PDB_dbChainId...")

        if res is not None:
            table = row_selector(table, 'PDB_dbResNum', res)
            log.info("SIFTS table filtered by PDB_dbResNum...")

        if uniprot is not None:
            table = row_selector(table, 'UniProt_dbAccessionId', uniprot)
            log.info("SIFTS table filtered by UniProt_dbAccessionId...")

        if site is not None:
            table = row_selector(table, 'UniProt_dbResNum', site)
            log.info("SIFTS table filtered by UniProt_dbResNum...")

        if table.empty:
            raise ValueError('{} resulted in an empty DataFrame...'.format(filename))

        return table

    def download(self, identifier=None, filename=None,
                 overwrite=False, decompress=True):
        """
        Downloads a SIFTS xml from the EBI FTP to the filesystem.

        :param identifier: (str) PDB accession ID
        :param filename: path to the SIFTS file
        :param overwrite: (boolean)
        :param decompress: (boolean) Decompresses the file
        :return: (side effects) output file path
        """

        identifier = self._get_identifier(identifier)
        filename = self._get_filename(filename)

        OutputFileHandler(filename, overwrite=overwrite)

        url_root = config.ftp_sifts
        url_endpoint = "{}.xml.gz".format(identifier)
        url = url_root + url_endpoint
        Downloader(url=url, filename=filename,
                   decompress=decompress, overwrite=overwrite)
