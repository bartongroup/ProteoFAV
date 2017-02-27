#!/usr/bin/env python
# -*- coding: utf-8

import logging
import re

import pandas as pd
import pymol  ##TODO: This import kicks of pymol, consider it for in function.
from analysis.clustering import clustered_table_to_partition_and_points, add_clusters_to_points
from pymol.cgo import *
from scipy.spatial import ConvexHull

from analysis.pymol_scripts import drawBoundingBox


def view_table(table, show=None, show_group_by=None, biological_assembly=True):
    """
    View the PDB entry associated with the provided structure table in PyMol and optionally create residue selections
    based on annotation columns.

    :param table: A structure table from `main.merge_tables`
    :param show: List of annotation columns to produce a selection passed on a boolean test (e.g., 'has_variant`)
    :param show_group_by: List of annotation columns to produce a selection based on the annotations value
    :param biological_assembly: Show the biological assembly.
    :return:
    """

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
    selections = []
    if show:
        for column in show:
            ## TODO: either parametrise test (e.g. allow .isnull() or create inverse selection too.
            select_atom = table[column].notnull()
            select_name = 'has_' + column

            # Selection must be built per chain to avoid ambiguous selections
            for chain in unique_chains:
                select_ResNums = residueIds[select_atom & (chain_ids == chain)]  # Get correct residues
                selections.append(make_selection(chain, select_ResNums, select_name))

    # Create level selections
    if show_group_by:
        for column in show_group_by:
            unique_values = table[column].dropna().unique()
            for value in unique_values:
                select_atom = (table[column] == value)
                select_name = column + '_{}'.format(value)
                for chain in unique_chains:
                    select_ResNums = residueIds[select_atom & (chain_ids == chain)]  # Get correct residues
                    selections.append(make_selection(chain, select_ResNums, select_name))

    selections = set(selections)

    # Now style them
    if len(selections) > 0:
        for select_name in selections:
            style_selection(select_name)


def make_selection(chain, select_ResNums, select_name):
    """
    Create a new selection or add residues to an existing selection in PyMol.

    :param chain: A Chain ID
    :param select_ResNums: PDB residue numbers
    :param select_name: The name for the selection
    :return:
    """
    # Sanitise select_name
    select_name = re.sub("[!@#$%^&*()'\"[\]{}\|~`<>.?/ ]+", "_", select_name.strip())

    # Create or append to the selection
    if select_name not in pymol.cmd.get_names('selections'):
        pymol.cmd.select('"{0}"'.format(select_name), 'none')
        logging.debug('Created selection "{0}" in PyMol'.format(select_name))
    selection = 'chain ' + chain + ' and resi ' + '+'.join(
        select_ResNums) + ' or ' + select_name
    pymol.cmd.select(select_name, selection)
    message = 'Added residues "{0}" to PyMol selection "{1}"'.format(selection, select_name)
    logging.debug(message)

    return select_name


def style_selection(select_name):
    # Apply some styles
    pymol.cmd.show("lines", select_name)
    pymol.util.cnc(select_name)

    # # Show bounding surface
    # pymol.cmd.flag('ignore', 'not ' + select_name, 'set')
    # pymol.cmd.delete('indicate')
    # pymol.cmd.show('surface')
    # pymol.cmd.rebuild()
    # pymol.cmd.flag('ignore', 'all', 'reset')

    drawBoundingBox(select_name)


def int_to_chain(x):
    if x <= 25:
        return chr(x + ord('A'))
    else:
        y = (x / 26)
        return chr((y - 1) + ord('A')) + int_to_chain(x - (26 * y))


def plot_convex_hulls(clustered_table):
    # Extract points for convex hull calculation
    part, points = clustered_table_to_partition_and_points(clustered_table)
    point_table = add_clusters_to_points(part, points)

    # Don't bother with clusters of less than four points
    counts = point_table.cluster_id.value_counts()
    clusters_to_drop = list(counts[counts < 4].index)
    for cluster_id in clusters_to_drop:
        point_table = point_table[point_table['cluster_id'] != cluster_id]

    grouped_point_tables = point_table.groupby('cluster_id')

    # Calculate convex hulls
    cluster_hulls = []
    for _, cluster_point_table in grouped_point_tables:
        _, cluster_points = clustered_table_to_partition_and_points(cluster_point_table)
        cluster_hulls.append(ConvexHull(cluster_points))

    # Prepare hull CGO lists
    hull_cgo_lists = []
    for hull, (_, cluster_point_table) in zip(cluster_hulls, grouped_point_tables):
        # Calculate facets
        cgo_facets = []
        for simplex in hull.simplices:
            _, cluster_points = clustered_table_to_partition_and_points(cluster_point_table)
            cgo = build_cgo(cluster_points[simplex, ])
            cgo_facets.append(cgo)
        hull_cgo_lists.append(concat_cgo_lists(cgo_facets))

    # Plot them all
    for i, cgo in enumerate(hull_cgo_lists):
        pymol.cmd.load_cgo(cgo, 'hull_' + str(i))


def build_cgo(array):
    cgo = [
        BEGIN, LINES,
        COLOR, 1., 1., 1.
    ]

    # Draw lines between points
    for start, finish in zip(array[:-1], array[1:]):
        cgo.append(VERTEX)
        for component in start:
            cgo.append(component)
        cgo.append(VERTEX)
        for component in finish:
            cgo.append(component)

    # Once more to join ends
    start, finish = array[-1], array[0]
    cgo.append(VERTEX)
    for component in start:
        cgo.append(component)
    cgo.append(VERTEX)
    for component in finish:
        cgo.append(component)

    cgo.append(END)

    return cgo


def concat_cgo_lists(list_of_cgo_lists):
    concat_cgo = []
    for i, cgo in enumerate(list_of_cgo_lists):
        i += 1
        if i == 1:
            mod_cgo = cgo[:-1]  # Remove END element
        elif i < len(list_of_cgo_lists):
            mod_cgo = cgo[6:-1]  # Remove BEGIN, LINES, COLOR and rgb elements and END
        else:
            mod_cgo = cgo[6:]  # Remove BEGIN, LINES, COLOR and rgb elements
        concat_cgo.append(mod_cgo)
    concat_cgo[:] = [e for cgo in concat_cgo for e in cgo]
    return concat_cgo

