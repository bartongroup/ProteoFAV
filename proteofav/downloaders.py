# -*- coding: utf-8

import os
import gzip
import shutil
import logging

from proteofav.fetchers import get_preferred_assembly_id

from proteofav.config import defaults as config

log = logging.getLogger('proteofav')


class Downloader(object):
    def __init__(self, url, filename, decompress=True, override=False):
        """
        :param url: (str) Full web-address
        :param filename: (str) Output filename
        :param decompress: (boolean) Decompresses the file
        :param override: (boolean) Overrides any existing file, if available
        """

        self.url = url
        self.outputfile = filename
        self.outputfile_origin = filename
        self.decompress = decompress
        self.override = override

        if self.decompress:
            if self.outputfile_origin.endswith('.gz'):
                self.outputfile = self.outputfile_origin.rstrip('.gz')

        if not os.path.exists(self.outputfile) or self.override:
            self._download()
            if self.outputfile_origin.endswith('.gz') and self.decompress:
                self._decompress()
        else:
            log.info("%s already available...", self.outputfile)

    def _download(self):
        try:
            try:
                import urllib.request
                from urllib.error import URLError, HTTPError
                with urllib.request.urlopen(self.url) as response, \
                        open(self.outputfile_origin, 'wb') as outfile:
                    shutil.copyfileobj(response, outfile)
            except (AttributeError, ImportError):
                import urllib
                urllib.urlretrieve(self.url, self.outputfile_origin)
        except (URLError, HTTPError, IOError, Exception) as e:
            log.debug("Unable to retrieve %s for %s", self.url, e)

    def _decompress(self):
        with gzip.open(self.outputfile_origin, 'rb') as infile, \
                open(self.outputfile, 'wb') as outfile:
            shutil.copyfileobj(infile, outfile)
            os.remove(self.outputfile_origin)
            log.info("Decompressed %s to %s",
                     self.outputfile_origin, self.outputfile)


def download_structure_from_pdbe(identifier, pdb=False, bio=False, override=False):
    """
    Downloads a structure from the PDBe to the filesystem.

    :param identifier: (str) PDB ID
    :param pdb: (boolean) PDB formatted if True, otherwise mmCIF format
    :param bio: (boolean) if true downloads the preferred Biological Assembly
    :param override: (boolean)
    :return: (side effects)
    """

    if pdb:
        filename = "{}.pdb".format(identifier)
    else:
        if bio:
            filename = "{}_bio.cif.gz".format(identifier)
        else:
            filename = "{}.cif".format(identifier)

    outputfile = os.path.join(config.db_root, config.db_pdbx, filename)
    os.makedirs(os.path.join(config.db_root, config.db_pdbx), exist_ok=True)

    if pdb:
        url_endpoint = "entry-files/download/pdb{}.ent".format(identifier)
    else:
        if bio:
            # atom lines only?
            # url_endpoint = ("static/entry/download/"
            #                "{}-assembly-{}_atom_site.cif.gz".format(identifier, pref))
            pref = get_preferred_assembly_id(identifier=identifier)
            url_endpoint = ("static/entry/download/"
                            "{}-assembly-{}.cif.gz".format(identifier, pref))
        else:
            # original mmCIF?
            # url_endpoint = "entry-files/download/{}.cif".format(pdbid)
            url_endpoint = "entry-files/download/{}_updated.cif".format(identifier)

    url_root = config.http_pdbe
    url = url_root + url_endpoint
    Downloader(url=url, filename=outputfile,
               decompress=True, override=override)
    return


def download_sifts_from_ebi(identifier, override=False):
    """
    Downloads a SIFTS xml from the EBI FTP to the filesystem.

    :param identifier: (str) PDB ID
    :param override: (boolean)
    :return: (side effects)
    """

    filename = "{}.xml.gz".format(identifier)
    outputfile = os.path.join(config.db_root, config.db_sifts, filename)
    os.makedirs(os.path.join(config.db_root, config.db_sifts), exist_ok=True)

    url_root = config.ftp_sifts
    url_endpoint = "{}.xml.gz".format(identifier)
    url = url_root + url_endpoint
    Downloader(url=url, filename=outputfile,
               decompress=True, override=override)
    return


def download_dssp_from_cmbi(identifier, override=False):
    """
    Downloads a pre-computed DSSP from the CMBI Netherlands FTP
    to the filesystem.

    :param identifier: (str) PDB ID
    :param override: (boolean)
    :return: (side effects)
    """

    filename = "{}.dssp".format(identifier)
    outputfile = os.path.join(config.db_root, config.db_dssp, filename)
    os.makedirs(os.path.join(config.db_root, config.db_dssp), exist_ok=True)

    url_root = config.ftp_dssp
    url_endpoint = "{}.dssp".format(identifier)
    url = url_root + url_endpoint
    Downloader(url=url, filename=outputfile,
               decompress=True, override=override)
    return

