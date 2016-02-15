#!/usr/bin/env python
# -*- coding: utf-8


import logging
import re

import pandas as pd
import pymol

from proteofav.main import merge_tables
from proteofav.variants import (_fetch_uniprot_variants,
                                _fetch_variant_characteristics_from_identifiers)

__author__ = 'smacgowan'


def visualise(pdb_id, assembly=False, use_ensembl=False, use_uniprot=False):
    """

    :param pdb_id:
    :param group_by_trait:
    :param use_ensembl:
    :param use_uniprot:
    :return:
    """

    def build_selection(chain, select_name, in_group, mapped_variants):
        """

        :param chain:
        :param group:
        :param in_group:
        :param mapped_variants:
        :return:
        """

        # Find the residues we want to highlight
        in_chain = mapped_variants.chain_id == chain
        start = mapped_variants.PDB_dbResNum[in_group & in_chain]
        variant_residues = list(start.dropna().astype(int).astype(str).unique())

        # Construct PyMol select command
        # Sanitisation of the selection name for pymol is important
        # or the name existence test will always fail and we'll over write
        # entries from previous chains!
        select_name = re.sub("[\W\d]+", "_", select_name.strip())
        if select_name not in pymol.cmd.get_names('selections'):
            pymol.cmd.select(select_name, 'none')
        selection = 'chain ' + chain + ' and resi ' + '+'.join(
                variant_residues) + ' or ' + select_name

        # Make the selection
        message = "Creating PyMol selection {} from '{}'".format(select_name,
                                                                 selection)
        logging.debug(message)
        pymol.cmd.select(select_name, selection)

        # Apply some styles
        pymol.cmd.show("lines", select_name)
        pymol.util.cnc(select_name)

    # Open the PDB as requested with PyMol and apply a few styles
    pymol.finish_launching()
    if assembly:
        pymol.cmd.fetch(pdb_id, type='pdb1')
        pymol.cmd.set('all_states', 'on')
    else:
        pymol.cmd.fetch(pdb_id)

    pymol.cmd.hide("everything", "all")
    pymol.cmd.show("ribbon", "all")
    pymol.util.cbc()

    # Get variants for all chains
    residue_mappings = merge_tables(pdb_id=pdb_id, chain='all',
                                    add_variants=True)
    has_variant = residue_mappings.start.notnull()

    # If we're going to make selections from ensembl traits get that data now
    if use_ensembl:
        variant_ids = residue_mappings.id_y[has_variant]
        traits = ensembl_traits(variant_ids)

    # Now create a PyMol selection for the variants on each chain
    chains = pymol.cmd.get_chains(pdb_id)
    for chain in chains:
        group = 'chain_' + chain
        in_chain = residue_mappings.chain_id == chain
        build_selection(chain, select_name=group, in_group=in_chain & has_variant,
                        mapped_variants=residue_mappings)

        if use_ensembl:
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
            structure_uniprots = residue_mappings.UniProt_dbAccessionId[in_chain]
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


def ensembl_traits(variant_ids):
    """

    :param variant_ids:
    :return:
    """
    # For now, need to iterate with GET requests
    # until POST can retrieve phenotypes
    traits = []
    for variant in variant_ids:
        phenos = \
            _fetch_variant_characteristics_from_identifiers(str(variant))[
                'phenotypes']
        if phenos == []:
            traits.append([variant, 'No_Available_Phenotype'])
        else:
            for entry in phenos:
                trait = entry['trait']
                traits.append([variant, trait])
    return traits


if __name__ == '__main__':
    visualise('3tnu', assembly=True, use_ensembl=True, use_uniprot=True)
