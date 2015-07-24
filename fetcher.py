#!/usr/bin/env python
# -*- coding: utf-8

__author__ = 'tbrittoborges'
"""
Created on 11:41 24/07/15 2015 

"""
import urllib


def fetch_pdb_files(pdb_id, destination, data_sources=("cif", "dssp", "sifts")):
    """

    :param pdb_id:
    :param destination:
    :param data_sources:
    :return:
    """

    if "cif" in data_sources:
        file_name = pdb_id + "_updated.cif"
        cif_url = "http://www.ebi.ac.uk/pdbe/entry-files/"
        urllib.urlretrieve(cif_url + file_name,
                           destination + file_name)
    if "dssp" in data_sources:
        file_name = pdb_id + ".dssp"
        dssp_url = "ftp://ftp.cmbi.ru.nl//pub/molbio/data/dssp/"
        urllib.urlretrieve(dssp_url + file_name,
                           destination + file_name)
    if "sifts" in data_sources:
        file_name = pdb_id + ".xml.gz"
        sifts_url = "ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/xml/"
        urllib.urlretrieve(sifts_url + file_name,
                           destination + file_name)


if __name__ == '__main__':
    fetch_pdb_files("4ibw", "/Users/tbrittoborges/Downloads/")
