#!/usr/bin/env python
# -*- coding: utf-8

__author__ = 'fabiomadeira'
"""
Created on 11:06 28/07/15 2015

"""

import shlex
from os import path
from collections import OrderedDict


def _mmcif_info_to_dict(filename):
    """Parse a mmCIF file and return a ordered dictionaries
    with all the category groups, data categories/items.

    :param filename: input CIF file
    :return: OrderedDict with main and leading keys
    """

    if not path.isfile(filename):
        raise IOError('File {} not found or unavailable.'.format(filename))

    def _mmcif_tokenize(handle):
        for line in handle:
            if line.startswith("#"):
                continue
            elif line.startswith(";"):
                token = line[1:].strip()
                for line in handle:
                    line = line.strip()
                    if line == ';':
                        break
                    token += line
                yield token
            else:
                try:
                    tokens = shlex.split(line)
                except ValueError:
                    # error "No closing quotation"
                    line = line.replace("'", '"')
                    # if odd - add a closing " to that line
                    if not line.count('"') % 2 == 0:
                        line = '{}"'.format(line)
                    tokens = shlex.split(line)
                for token in tokens:
                    yield token

    ordict = OrderedDict()
    with open(filename) as handle:
        loop_flag = False
        key = None
        tokens = _mmcif_tokenize(handle)
        token = next(tokens)
        ordict[token[0:5]] = token[5:]
        i = 0
        n = 0
        for token in tokens:
            if token == "loop_":
                loop_flag = True
                keys = []
                i = 0
                n = 0
                continue
            elif loop_flag:
                if token.startswith("_"):
                    if i > 0:
                        loop_flag = False
                    else:
                        ordict[token] = []
                        keys.append(token)
                        n += 1
                        continue
                else:
                    ordict[keys[i % n]].append(token)
                    i += 1
                    continue
            if key is None:
                key = token
            else:
                ordict[key] = token
                key = None

    full = OrderedDict()
    for full_key in ordict:
        # skipping first mmCIF entry and ATOM/HETATOM lines
        if full_key != "data_" and not full_key.startswith('_atom_site'):
            main_key = (full_key.split(".")[0]).lstrip('_')
            lead_key = full_key.split(".")[1]
            if main_key not in full:
                full[main_key] = {}
            full[main_key][lead_key] = ordict[full_key]
    return full

if __name__ == '__main__':
    X = _mmcif_info_to_dict("tests/CIF/2pah.cif")
    print(X['exptl_crystal'])
    print(X['pdbx_struct_assembly'])
    print(X['pdbx_struct_assembly']['oligomeric_details'])
    print(X['pdbx_struct_assembly_gen']['asym_id_list'])
    print(X['pdbx_struct_oper_list']['symmetry_operation'])
    pass