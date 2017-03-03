#!/usr/bin/env python
# -*- coding: utf-8
"""
Created on 18:07 28/02/2017 2017 
Methods for generating files for third party visualisation solutions.
"""
import logging

from proteofav.config import defaults

log = logging.getLogger('proteofav.config')

def make_chimera_attribute_file(column, recipient='residues', match_mode="1 - to - 1"):
    """
    Write a UFSC Chimera attribute file. It assumes that the pd.Series index
        matches to the residue number in the target protein structure.
    :param pd.Series column: column with attributes
    :param str recipient: level of attribute assignment. Choice from atoms,
        residues and molecules.
    :param str match_mode:
    :return None:

    .. example:
        generate_chimera_attrFile(variants.applymap(counts).tcga)

    .. note: Documentation
        https://www.cgl.ucsf.edu/chimera/docs/ContributedSoftware/defineattrib/defineattrib.html

    """
    if recipient not in ['atoms', 'residues', 'molecules']:
        raise ValueError('The recipient is not supported.')
    if recipient != 'residues':
        raise NotImplementedError

    try:
        name = column.name
    except AttributeError:  # column needs to be Series, not an DataFrame
        raise TypeError('Column parameter needs to be a pd.Series.')

    template = """# Generated with ProteoFAV
    attribute: {name}
    match mode: {match_mode}
    recipient: {recipient}
    """

    column = column.reset_index()
    formatter = ["\t:{}".format, "\t{}".format]
    lines = column.to_string(index=False, header=False, formatters=formatter)
    lines = lines.replace(" ", "")
    return template.format(name=name, match_mode=match_mode, recipient=recipient) + lines


def make_chimera_command_file(table, color_secondary_structure=True):
    """
    Creates a Chimera command file for protein structure visualisation.

    :param table:
    :param color_secondary_structure:
    :return:

    .. TODO:: Select models, chains, residues and select atoms:
        display #{model_n}:{start}-{stop}.{chain}@ca

    """
    defaults.db_mmcif
    line = "open {filename}"


    if color_secondary_structure:
        closing_lines = "color green,a helix\ncolor yellow,a strand\ncolor gray,a coil\n"


def write_file(content, file):
    """


    :param str content:
    :param str file:
    """
    with open(file, 'w') as open_file:
        open_file.write(content)
