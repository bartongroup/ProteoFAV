#!/usr/bin/env python
# -*- coding: utf-8


import re
import pymol
import logging
import pandas as pd

from main.to_table import merge_tables
from variants.to_table import (_fetch_uniprot_variants,
                               _variant_characteristics_from_identifiers)

__author__ = 'smacgowan'


def visualise(pdb_id, assembly=False, use_ensembl=False, use_uniprot=False):
    """

    :param pdb_id:
    :param group_by_trait:
    :param use_ensembl:
    :param use_uniprot:
    :return:
    """

    def build_selection(chain, group, in_group, mapped_variants):
        """

        :param chain:
        :param group:
        :param in_group:
        :param mapped_variants:
        :return:
        """
        start = mapped_variants.PDB_dbResNum[in_group]
        variant_residues = list(
            start[pd.notnull(start)].astype(int).astype(str).unique())
        # sanitisation of the selection name for pymol is important
        # or the name existence test will always fail and we'll over write
        # entries from previous chains!
        name = group + '_vars'
        name = re.sub("[\W\d]+", "_", name.strip())
        if name not in pymol.cmd.get_names('selections'):
            pymol.cmd.select(name, 'none')
        selection = 'chain ' + chain + ' and resi ' + '+'.join(
            variant_residues) + ' or ' + name
        # name = group + '_vars'
        message = "Creating PyMol selection {} from '{}'".format(name,
                                                                 selection)
        logging.debug(message)
        pymol.cmd.select(name, selection)

    # Open the PDB with pymol
    pymol.finish_launching()
    if assembly:
        pymol.cmd.fetch(pdb_id, type='pdb1')
        pymol.cmd.set('all_states', 'on')
    else:
        pymol.cmd.fetch(pdb_id)
    pymol.cmd.hide("everything", "all")
    pymol.cmd.show("ribbon", "all")
    pymol.util.cbc()

    # Get available variants for any chain
    chains = pymol.cmd.get_chains(pdb_id)

    for chain in chains:

        # Get the variants for this chain
        residue_mappings = merge_tables(pdb_id=pdb_id, chain=chain,
                                        add_variants=True)
        has_variant = residue_mappings.start.notnull()
        start = residue_mappings.PDB_dbResNum[has_variant]
        variant_residues = list(
            start[pd.notnull(start)].astype(int).astype(str).unique())

        # Create the selection
        selection = 'chain ' + chain + ' and resi ' + '+'.join(variant_residues)
        name = 'chain_' + chain + '_vars'
        message = "Creating PyMol selection {} from '{}'".format(name,
                                                                 selection)
        logging.debug(message)
        pymol.cmd.select(name, selection)

        # Apply some styles
        pymol.cmd.show("lines", name)
        pymol.util.cnc(name)

        if use_ensembl:
            variant_ids = residue_mappings.id_y[has_variant]

            # For now, need to iterate with GET requests
            # until POST can retrieve phenotypes
            traits = []
            for variant in variant_ids:
                phenos = \
                    _variant_characteristics_from_identifiers(str(variant))[
                        'phenotypes']
                if phenos == []:
                    traits.append([variant, 'No_Available_Phenotype'])
                else:
                    for entry in phenos:
                        trait = entry['trait']
                        traits.append([variant, trait])

            groups = set(zip(*traits)[1])
            for group in groups:
                variants_in_group = []
                for id, trait in traits:
                    if trait == group:
                        variants_in_group.append(id)

                # Get the relevant residues and create
                # selection for current group
                in_group = [True if x in variants_in_group else False for x in
                            residue_mappings.id_y]

                build_selection(chain, group, in_group, residue_mappings)

        if use_uniprot:

            # Extract first uniprot (REFACTOR THIS, also `merge_tables`)
            structure_uniprots = residue_mappings.UniProt_dbAccessionId
            structure_uniprots = structure_uniprots[
                structure_uniprots.notnull()]
            structure_uniprot = structure_uniprots.unique()[0]

            # Get the variants and merge onto residue_mappings
            uniprot_vars = _fetch_uniprot_variants(structure_uniprot)
            # Not ideal as UniProt_dbResNum  is not unique but is fit for
            # purpose

            mapped = pd.merge(residue_mappings, uniprot_vars,
                              on='UniProt_dbResNum')

            # Create a selection for each trait
            groups = mapped.disease.unique()
            for group in groups:
                in_group = mapped.disease == group

                build_selection(chain, group, in_group, mapped)


if __name__ == '__main__':
    visualise('3tnu', assembly=True, use_uniprot=True)
