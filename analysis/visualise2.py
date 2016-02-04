#!/usr/bin/env python
# -*- coding: utf-8

import pymol

def view_table(table, show=None, biological_assembly=True):

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

