# -*- coding: utf-8

"""
Created on 18:07 28/02/2017 2017 
Methods aiming to integrate ProteoFAV table to viewers such as UFSC Chimera and Jalview.
"""

import os
import logging

log = logging.getLogger('proteofav.config')


def make_chimera_attribute_file(column, recipient='residues', match_mode="1 - to - 1"):
    """
    Write a UFSC Chimera attribute file. It assumes that the pd.Series index
        matches to the residue number in the target protein structure.
    :param pd.Series column: column with attributes in data and recipient in the index
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


def make_chimera_command_file(filename, content='', color_secondary_structure=True):
    """
    Creates a Chimera command file for protein structure visualisation.

    :param str filename: path to mmcif file
    :param content: personalised content of the chimera command file
    :type content: str or iterable
    :param bool color_secondary_structure: wheter to append secondary structure coloring
    :return str: lines to be written to file

    .. todo:: Select models, chains, residues and select atoms:
        display #{model_n}:{start}-{stop}.{chain}@ca

    """

    line = "open {filename}\n".format(filename=filename)

    if isinstance(content, str):
        line += content
    elif hasattr(content, '__iter__'):
        line += "\n".join(content)

    if color_secondary_structure:
        line += "color green,r helix\ncolor yellow,r strand\ncolor gray,r coil\n"

    return line


def visualise_chimera(filename, column, **kwargs):
    """
    Creates an attribute file and a command file to observed a ProteoFAV processed data in UFSC
    Chimera.

    :param str filename: path to mmcif file
    :param pd.Series column: column with attributes in data and recipient in the index
    :return None: callback
    """

    attribute_file_content = make_chimera_attribute_file(column, **kwargs)
    name = column.name
    write_file(attribute_file_content, '{}.chimera_attrFile'.format(name))

    attribute_line = 'defattr {}.chimera_attrFile'.format(name)
    command_file_content = make_chimera_command_file(filename, content=attribute_line)
    basename = os.path.basename(filename)
    write_file(command_file_content, "{}.com".format(basename.split('.')[0]))


def write_file(content, filename):  # pragma: no cover
    """
    Write content to a new file.

    :param str content: lines of the file as a string
    :param str filename: path to the file
    """
    with open(filename, 'w') as open_file:
        open_file.write(content)
