__author__ = 'smacgowan'

import pymol
import pandas as pd

from main import merge_tables


def visualise(pdb_id):

    # Open the PDB with pymol
    pymol.finish_launching()
    pymol.cmd.fetch(pdb_id)
    pymol.cmd.hide("everything", "all")
    pymol.cmd.show("ribbon", "all")
    pymol.util.cbc()

    # Get available variants for any chain
    chains = pymol.cmd.get_chains(pdb_id)

    # Get the variants
    for chain in chains:
        residue_mappings = merge_tables(pdb_id=pdb_id, chain=chain, add_variants=True)
        has_variant = residue_mappings.start.notnull()
        start = residue_mappings.PDB_dbResNum[has_variant]
        variant_residues = list(start[pd.notnull(start)].astype(int).astype(str).unique())

        # Create the selection
        selection = 'chain ' + chain + ' and resi ' + '+'.join(variant_residues)
        name = 'chain_' + chain + '_vars'
        pymol.cmd.select(name, selection)

        # Apply some styles
        pymol.cmd.show("sticks", name)
        pymol.util.cnc(name)

