# -*- coding: utf-8 -*-

from os import path
import logging
import pandas as pd

from proteofav.structures import select_cif
from proteofav.utils import fetch_files, row_selector
from proteofav.library import scop_3to1


from proteofav.config import defaults

log = logging.getLogger('proteofav.config')

__all__ = ['_dssp', '_import_dssp_chains_ids', 'select_dssp']


def _dssp(filename):
    """
    Parses DSSP file output to a pandas DataFrame.

    :param filename: input SIFTS xml file
    :return: pandas table dataframe
    """
    # column width descriptors
    cols_widths = ((0, 5), (6, 11), (11, 12), (13, 14), (16, 17), (35, 38),
                   (103, 109), (109, 115))
    # simplified headers for the table
    dssp_header = ("dssp_index", "icode", "chain_id", "aa", "ss", "acc", "phi", "psi")
    dssp_table = pd.read_fwf(filename,
                             skiprows=28,
                             names=dssp_header,
                             colspecs=cols_widths,
                             index_col=0,
                             compression=None)
    if dssp_table.empty:
        log.error('DSSP file {} resulted in a empty Dataframe'.format(filename))
        raise ValueError('DSSP file {} resulted in a empty Dataframe'.format(
            filename))
    return dssp_table


def _import_dssp_chains_ids(pdb_id):
    """Imports mmCIF chain identifier to DSSP.

    :param pdb_id:
    :return: DSSP table with corrected chain ids.
    """
    dssp_table = select_dssp(pdb_id)
    cif_table = select_cif(pdb_id)
    cif_seq = cif_table.auth_comp_id.apply(scop_3to1.get)
    dssp_has_seq = dssp_table.aa.isin(scop_3to1.values())
    dssp_seq = dssp_table.aa[dssp_has_seq]
    # Import only if the sequences are identical
    if not (cif_seq == dssp_seq).all():
        err = ('Inconsitent DSSP / mmCIF sequence for {} protein structure cannot be fixed'
               'by import_dssp_chains_ids')
        raise ValueError(err.format(pdb_id))
    dssp_table.loc[dssp_has_seq, 'chain_id'] = cif_table.auth_asym_id
    return dssp_table


def select_dssp(pdb_id, chains=None):
    """
    Produce table from DSSP file output.

    :param pdb_id: PDB identifier
    :param chains: PDB protein chain
    :return: pandas dataframe
    """
    dssp_path = path.join(defaults.db_dssp, pdb_id + '.dssp')
    try:
        dssp_table = _dssp(dssp_path)
    except IOError:
        dssp_path = fetch_files(pdb_id, sources='dssp', directory=defaults.db_dssp)[0]
        dssp_table = _dssp(dssp_path)
    except StopIteration:
        raise IOError('{} is unreadable.'.format(dssp_path))

    if chains:
        try:
            dssp_table = row_selector(dssp_table, 'chain_id', chains)
        except ValueError:
            # TODO:
            # Could not find the correct PDB chain. It happens for protein structures with complex
            # chain identifier, as 4v9d.
            # dssp_table = _import_dssp_chains_ids(pdb_id)
            # dssp_table = row_selector(dssp_table, 'chain_id', chains)
            log.error('Error loading DSSP file: Chain {} not in {}'.format(chains, pdb_id))
            return None
    # remove dssp line of transition between chains
    dssp_table = dssp_table[dssp_table.aa != '!']

    dssp_table.reset_index(inplace=True)
    if dssp_table.duplicated(['icode', 'chain_id']).any():
        log.info('DSSP file for {} has not unique index'.format(pdb_id))
    try:
        dssp_table.loc[:, 'icode'] = dssp_table.loc[:, 'icode'].astype(int)
    except ValueError:
        log.warning("{} insertion code detected in the DSSP file.".format(pdb_id))
        dssp_table.loc[:, 'icode'] = dssp_table.loc[:, 'icode'].astype(str)

    return dssp_table
