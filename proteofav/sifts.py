# -*- coding: utf-8

import logging

from proteofav.parsers import parse_sifts_residues_from_file
from proteofav.downloaders import Downloader

from proteofav.utils import InputFileHandler
from proteofav.utils import OutputFileHandler
from proteofav.utils import exclude_columns
from proteofav.utils import constrain_column_types
from proteofav.library import sifts_types

from proteofav.config import defaults as config

log = logging.getLogger('proteofav')


class SIFTS(object):
    def read(self, filename, excluded_cols=None, add_regions=True, add_dbs=False):
        """
        Parses the residue fields of a SIFTS XML file.

        :param filename: path to the SIFTS file
        :param excluded_cols: option to exclude SIFTS dbSources
        :param add_regions: boolean
        :param add_dbs: boolean
        :return: returns a pandas DataFrame
        """
        InputFileHandler(filename)

        if excluded_cols is None:
            excluded_cols = ("InterPro", "GO", "EC", "NCBI")

        table = parse_sifts_residues_from_file(filename, excluded_cols=excluded_cols,
                                               add_regions=add_regions, add_dbs=add_dbs)

        # excluding columns
        table = exclude_columns(table, excluded=excluded_cols)

        # enforce some specific column types
        table = constrain_column_types(table, sifts_types)

        if table.empty:
            raise ValueError('{} resulted in an empty DataFrame...'.format(filename))

        return table

    def download(self, identifier, filename, overwrite=False, decompress=True):
        """
        Downloads a SIFTS xml from the EBI FTP to the filesystem.

        :param identifier: (str) PDB accession ID
        :param filename: path to the SIFTS file
        :param overwrite: (boolean)
        :param decompress: (boolean) Decompresses the file
        :return: (side effects) output file path
        """

        OutputFileHandler(filename, overwrite=overwrite)

        url_root = config.ftp_sifts
        url_endpoint = "{}.xml.gz".format(identifier)
        url = url_root + url_endpoint
        Downloader(url=url, filename=filename,
                   decompress=decompress, overwrite=overwrite)
