# -*- coding: utf-8 -*-

from os import path
import logging
import pandas as pd
try:
    # python 2.7
    from StringIO import StringIO
except ImportError:
    from io import StringIO

from proteofav.structures import select_cif
from proteofav.utils import (fetch_files, row_selector, InputFileHandler,
                             constrain_column_types, exclude_columns)
from proteofav.library import scop_3to1, dssp_types


from proteofav.config import defaults

log = logging.getLogger('proteofav.config')

__all__ = ['parse_dssp_residues', '_import_dssp_chains_ids', 'select_dssp']


def parse_dssp_residues(filename, excluded_cols=None):
    """
    Parse lines of the DSSP file to get entries for every Residue
    in each CHAIN. The hierarchy is maintained. CHAIN->RESIDUE->[...].

    :param filename: path to the DSSP file
    :param excluded_cols: list of columns to be excluded
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
    header = ("LINE", "RES", "RES_FULL", "INSCODE", "CHAIN", "AA", "SS", "STRUCTURE",
              "BP1", "BP2", "BP2_CHAIN", "ACC",
              "NH_O_1", "NH_O_1_nrg", "O_HN_1", "O_HN_1_nrg",
              "NH_O_2", "NH_O_2_nrg", "O_HN_2", "O_HN_2_nrg",
              "TCO", "KAPPA", "ALPHA", "PHI", "PSI",
              "X-CA", "Y-CA", "Z-CA")

    widths = ((0, 5), (5, 10), (5, 11), (10, 11), (11, 12), (12, 15), (16, 17), (17, 25),
              (25, 29), (29, 33), (33, 34), (34, 38),
              (38, 45), (46, 50), (50, 56), (57, 61),
              (61, 67), (68, 72), (72, 78), (79, 84),
              (85, 91), (91, 97), (97, 103), (103, 109), (109, 115),
              (115, 123), (123, 130), (130, 137))

    all_str = {key: str for key in header}
    table = pd.read_fwf(StringIO(lines), names=header, colspecs=widths,
                        compression=None, converters=all_str, keep_default_na=False)

    # excluding columns
    if excluded_cols is None:
        excluded_cols = ("LINE", "STRUCTURE", "BP1", "BP2", "BP2_CHAIN",
                         "NH_O_1", "NH_O_1_nrg", "O_HN_1", "O_HN_1_nrg",
                         "NH_O_2", "NH_O_2_nrg", "O_HN_2", "O_HN_2_nrg",
                         "X-CA", "Y-CA", "Z-CA")

    table = exclude_columns(table, excluded=excluded_cols)

    # enforce some specific column types
    table = constrain_column_types(table, col_type_dict=dssp_types)

    if table.empty:
        log.error('DSSP file {} resulted in a empty Dataframe'.format(filename))
        raise ValueError('DSSP file {} resulted in a empty Dataframe'.format(
            filename))
    return table


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
    dssp_table.loc[dssp_has_seq, 'CHAIN'] = cif_table.auth_asym_id
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
        dssp_table = parse_dssp_residues(dssp_path)
    except IOError:
        dssp_path = fetch_files(pdb_id, sources='dssp', directory=defaults.db_dssp)[0]
        dssp_table = parse_dssp_residues(dssp_path)
    except StopIteration:
        raise IOError('{} is unreadable.'.format(dssp_path))

    if chains:
        try:
            dssp_table = row_selector(dssp_table, 'CHAIN', chains)
        except ValueError:
            # TODO:
            # Could not find the correct PDB chain. It happens for protein structures with complex
            # chain identifier, as 4v9d.
            # dssp_table = _import_dssp_chains_ids(pdb_id)
            # dssp_table = row_selector(dssp_table, 'CHAIN', chains)
            log.error('Error loading DSSP file: Chain {} not in {}'.format(chains, pdb_id))
            return None
    # remove dssp line of transition between chains
    dssp_table = dssp_table[dssp_table.AA != '!']

    dssp_table.reset_index(inplace=True)
    if dssp_table.duplicated(['RES_FULL', 'CHAIN']).any():
        log.info('DSSP file for {} has not unique index'.format(pdb_id))
    try:
        dssp_table.loc[:, 'RES_FULL'] = dssp_table.loc[:, 'RES_FULL'].astype(int)
    except ValueError:
        log.warning("{} insertion code detected in the DSSP file.".format(pdb_id))
        dssp_table.loc[:, 'RES_FULL'] = dssp_table.loc[:, 'RES_FULL'].astype(str)

    return dssp_table
