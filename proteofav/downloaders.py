# -*- coding: utf-8

import os
import gzip
import shutil
import logging
import tempfile

from proteofav.utils import InputFileHandler
from proteofav.utils import OutputFileHandler

from proteofav.config import defaults as config

log = logging.getLogger('proteofav')


class Downloader(object):
    def __init__(self, url, filename, decompress=True, overwrite=False):
        """
        :param filename: (str) Output _filename
        :param decompress: (boolean) Decompresses the file
        :param overwrite: (boolean) Overrides any existing file, if available
        """

        self._url = url
        self._filename = filename
        # self._tempfile = tempfile.NamedTemporaryFile(dir=config.db_tmp).name
        self._tempfile = tempfile.NamedTemporaryFile().name
        self._decompress = decompress
        self._override = overwrite

        OutputFileHandler(self._filename, overwrite=overwrite)
        OutputFileHandler(self._tempfile, overwrite=overwrite)

        self._download()
        if self._decompress:
            self._uncompress()

    def _download(self):
        try:
            try:
                import urllib.request
                from urllib.error import URLError, HTTPError
                with urllib.request.urlopen(self._url) as response, \
                        open(self._tempfile, 'wb') as outfile:
                    shutil.copyfileobj(response, outfile)
            except (AttributeError, ImportError):
                import urllib
                urllib.urlretrieve(self._url, self._tempfile)
            InputFileHandler(self._tempfile)
            with open(self._tempfile, 'rb') as infile, \
                    open(self._filename, 'wb') as outfile:
                shutil.copyfileobj(infile, outfile)
                os.remove(self._tempfile)
        except (URLError, HTTPError, IOError, Exception) as e:
            log.debug("Unable to retrieve %s for %s", self._url, e)

    def _uncompress(self):
        InputFileHandler(self._filename)
        with open(self._filename, 'rb') as infile, \
                open(self._tempfile, 'wb') as outfile:
            shutil.copyfileobj(infile, outfile)
        InputFileHandler(self._tempfile)
        with gzip.open(self._tempfile, 'rb') as infile, \
                open(self._filename, 'wb') as outfile:
            shutil.copyfileobj(infile, outfile)
            os.remove(self._tempfile)
            log.info("Decompressed %s", self._filename)
