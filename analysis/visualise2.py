#!/usr/bin/env python
# -*- coding: utf-8

import logging
import re
import pandas as pd
import pymol  ##TODO: This import kicks of pymol, consider it for in function.

def view_table(table, show=None, show_group_by=None, biological_assembly=True):

    # First figure out what PDB and chains are in the table
    pdb_id = table.PDB_dbAccessionId.unique()[0]

    # Open the PDB as requested with PyMol and apply a few styles
    pymol.finish_launching()
    if biological_assembly:
        pymol.cmd.fetch(pdb_id, type='pdb1')
        pymol.cmd.set('all_states', 'on')
    else:
        pymol.cmd.fetch(pdb_id)

    pymol.cmd.hide("everything", "all")
    pymol.cmd.show("ribbon", "all")
    pymol.util.cbc()

    # Create a pymol selection for any column names in `show` ----------

    # First find PDB_dbResNum (as could be column or index) and chain IDs
    if table.index.name == 'PDB_dbResNum':
        residueIds = pd.Series(table.index, index=table.index)
    else:
        residueIds = table.PDB_dbResNum

    # Get residue IDs in suitable format and drop residues with no PDB reference
    ##TODO: this may be useful elsewhere
    dropped = residueIds.isnull()
    residueIds = residueIds.dropna()
    residueIds = residueIds.astype('object').astype('int').astype('string')
    table = table[dropped == False]

    chain_ids = table.chain_id
    unique_chains = chain_ids.dropna().unique()

    # Now create labelled boolean selections
    if show:
        for column in show:
            ## TODO: either parametrise test (e.g. allow .isnull() or create inverse selection too.
            select_atom = table[column].notnull()
            select_name = 'has_' + column

            # Selection must be built per chain to avoid ambiguous selections
            for chain in unique_chains:
                select_ResNums = residueIds[select_atom & (chain_ids == chain)]  # Get correct residues
                make_selection(chain, select_ResNums, select_name)

    # Create level selections
    if show_group_by:
        for column in show_group_by:
            unique_values = table[column].dropna().unique()
            for value in unique_values:
                select_atom = (table[column] == value)
                select_name = column + '_{}'.format(value)
                for chain in unique_chains:
                    select_ResNums = residueIds[select_atom & (chain_ids == chain)]  # Get correct residues
                    make_selection(chain, select_ResNums, select_name)



def make_selection(chain, select_ResNums, select_name):
    """
    Create a new selection or add residues to an existing selection in PyMol.

    :param chain: A Chain ID
    :param select_ResNums: PDB residue numbers
    :param select_name: The name for the selection
    :return:
    """
    select_name = re.sub("[!@#$%^&*()'\"[\]{}\|~`<>.?/ ]+", "_", select_name.strip())  # Sanitise select_name
    if select_name not in pymol.cmd.get_names('selections'):
        pymol.cmd.select('"{0}"'.format(select_name), 'none')
        logging.debug('Created selection "{0}" in PyMol'.format(select_name))
    selection = 'chain ' + chain + ' and resi ' + '+'.join(
        select_ResNums) + ' or ' + select_name
    pymol.cmd.select(select_name, selection)
    message = 'Added residues "{0}" to PyMol selection "{1}"'.format(selection, select_name)
    logging.debug(message)
    # Apply some styles
    pymol.cmd.show("lines", select_name)
    pymol.util.cnc(select_name)

