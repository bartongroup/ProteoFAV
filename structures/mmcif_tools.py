#!/usr/bin/env python
# -*- coding: utf-8

__author__ = 'fabiomadeira'
"""
Created on 11:06 28/07/15 2015

"""

import shlex
from collections import OrderedDict
from os import path

import numpy as np
import pandas as pd

from structures.to_table import _mmcif_atom


def _mmcif_info_to_dict(filename):
    """
    Parse a mmCIF file and return a ordered dictionaries
    with all the category groups, data categories/items.
    Based on Biopython's MMCIFDict class.

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


def _bio_unit_parse_operation_expression(expression):
    """
    Parses symmetry operations needed to generate biological units.

    :param expression:
    :return: returns a list
    """
    operations = []
    stops = [",", "-", ")"]
    i = 1

    # iterate over the operation expression
    while i in range(1, len(expression) - 1):
        pos = i

        # read an operation
        while expression[pos] not in stops and pos < len(expression) - 1:
            pos += 1
        current_op = expression[i: pos]

        # handle single operations
        if expression[pos] != "-":
            operations.append(current_op)
            i = pos

        # handle ranges
        if expression[pos] == "-":
            pos += 1
            i = pos

            # read in the range's end value
            while expression[pos] not in stops:
                pos += 1
            end = int(expression[i: pos])

            # add all the operations in [currentOp, end]
            for val in range((int(current_op)), end + 1):
                operations.append(str(val))
            i = pos
        i += 1

    return operations


def _bio_unit_prepare_operation(operations, index1, index2):
    """
    Prepares the symmetry operations based on the operations list
    and two indexes.

    :param operations: operations list
    :param index1: operation index 1
    :param index2: operation index 2
    :return: returns an operation matrix
    """
    # prepare matrices for operations 1 & 2
    op1 = np.array([[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 1]], dtype=np.float_)
    op2 = np.array([[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 1]], dtype=np.float_)

    # fill the operation matrices for operations 1 & 2
    if type(operations['id']) is list:
        for i in range(3):
            op1[i][3] = float(operations["vector[{}]".format(str(i + 1))][index1])
            if index2 != -1:
                op2[i][3] = float(operations["vector[{}]".format(str(i + 1))][index2])
            for j in range(3):
                op1[i][j] = float(operations["matrix[{}][{}]".format(str(i + 1), str(j + 1))][index1])
                if index2 != -1:
                    op2[i][j] = float(operations["matrix[{}][{}]".format(str(i + 1), str(j + 1))][index2])
    else:
        for i in range(3):
            op1[i][3] = float(operations["vector[{}]".format(str(i + 1))])
            if index2 != -1:
                op2[i][3] = float(operations["vector[{}]".format(str(i + 1))])
            for j in range(3):
                op1[i][j] = float(operations["matrix[{}][{}]".format(str(i + 1), str(j + 1))])
                if index2 != -1:
                    op2[i][j] = float(operations["matrix[{}][{}]".format(str(i + 1), str(j + 1))])

    # handles non-Cartesian product expressions
    if index2 == -1:
        return op1

    operation = np.array([[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 1]], dtype=np.float_)

    # handles Cartesian product expressions (4x4 matrix multiplication)
    for row in range(4):
        for col in range(4):
            fsum = 0.0
            for r in range(4):
                fsum += (op1[row][r] * op2[r][col])
            operation[row][col] = fsum

    return operation


def _get_mmcif_bio_units(filename, most_likely=True):
    """
    Method to get the most likely entry from a set of different
    likely biological assemblies. If most_likely is false gets the
    first available biological unit, if any.

    :param filename: input CIF file
    :param most_likely: apply simple rules to find the most likely bio unit
    :return:
    """

    # get the information dictionary
    info = _mmcif_info_to_dict(filename)

    # describes possible macromolecular assemblies
    try:
        assembly = info["pdbx_struct_assembly"]
    except KeyError:
        # no biological assemblies are available
        return None, None, None
    n_assembly = {}

    # details the generation of each macromolecular assembly
    assembly_gen = info["pdbx_struct_assembly_gen"]
    n_assembly_gen = {}

    # details translation and rotation operations required to generate/transform
    # assembly coordinates
    oper_list = info["pdbx_struct_oper_list"]
    # we always need the full list of operations because these might be needed
    # to generate a particular (selected) most likely biological assembly
    # n_oper_list = {}

    # values in assembly are in stored in lists if there are more than one possible
    # biological unit
    if type(assembly['id']) is list:
        length = len(assembly['id'])
        index = 0
        if length > 1:
            if most_likely:
                # define which is most likely biological assembly
                tmp_indexer = 0
                for i in range(length):
                    # get molecule counts
                    o_count = int(assembly["oligomeric_count"][i])

                    # get details of the biological assembly
                    a_detail = assembly["details"][i]

                    # simple scoring each entry
                    if "author" in a_detail and "software" in a_detail:
                        indexer = o_count * 3
                    elif "author" in a_detail or "software" in a_detail:
                        indexer = o_count * 2
                    else:
                        indexer = o_count * 1

                    if indexer > tmp_indexer:
                        tmp_indexer = indexer
                        index = i
            else:
                index = 0

            for entry in assembly:
                n_assembly[entry] = assembly[entry][index]
            for entry in assembly_gen:
                n_assembly_gen[entry] = assembly_gen[entry][index]
            # for entry in oper_list:
            #     n_oper_list[entry] = oper_list[entry][index]

    else:
        # if only one assembly is available
        n_assembly = assembly
        n_assembly_gen = assembly_gen
        # n_oper_list = oper_list

    return n_assembly, n_assembly_gen, oper_list


def _bio_unit_to_table(filename, most_likely=True, method=1):
    """
    Generates a most likely biological assembly for the provided structure.

    :param filename: input CIF file
    :param most_likely: apply simple rules to find the most likely bio unit
    :param method: defaults to 1 = Chain id goes from A to A.1, A.2
                               2 = Model is incremented from 1 to 1 and 2
                               3 = A new column is added with a serial counter
    :return: new pandas table with biological units
    """

    # get the atom lines from mmCIF file
    atom_site_ref = _mmcif_atom(filename)
    # get attributes
    attributes = atom_site_ref.columns.values
    # add a new column with Biological Unit Counter
    attributes = [v for v in attributes]
    if method is 3:
        attributes.append('bio_unit_counter')

    # get the biological assemblies from mmCIF file
    assembly, assembly_gen, oper_list = _get_mmcif_bio_units(filename,
                                                             most_likely=most_likely)

    if not assembly and not assembly_gen and not oper_list:
        # no biological units available
        # returns unmodified mmcif atom table
        return atom_site_ref

    # lists to hold the individual operations
    oper = []
    oper2 = []

    # keep track of the current atom and model number
    atomNum = 1
    modelNum = 1
    asym_counter = 1
    atom_dict = OrderedDict()

    # get the operation expression for this assembly from the oper_expression attribute
    oper_expression = assembly_gen["oper_expression"]
    # count the number of left parentheses in the operation expression
    if "(" in oper_expression:
        parenCount = oper_expression.count("(")
    else:
        oper_expression = "(%s)" % oper_expression
        parenCount = oper_expression.count("(")

    # handles one operation assemblies (e.g., "1")
    if parenCount == 0:
        oper.append(oper_expression)

    # handles multiple operation assemblies, no Cartesian products (e.g., "(1-5)")
    if parenCount == 1:
        oper.extend(_bio_unit_parse_operation_expression(oper_expression))

    # handles Cartesian product expressions (e.g., "(X0)(1-60)")
    if parenCount == 2:
        # Break the expression into two parenthesized expressions and parse them
        temp = oper_expression.find(")")
        oper.extend(_bio_unit_parse_operation_expression(oper_expression[0:temp + 1]))
        oper2.extend(_bio_unit_parse_operation_expression(oper_expression[temp + 1:]))

    # retrieve the asym_id_list, which indicates which atoms to apply the operations to
    asym_id_list = assembly_gen["asym_id_list"]

    # either 1 or more operations
    temp = (1 > len(oper2)) and 1 or len(oper2)

    # for every operation in the first parenthesized list
    for op1 in oper:
        # Find the index of the current operation in the oper_list category table
        op1index = 0
        for row in range(len(oper_list['id'])):
            if oper_list["id"][row] == op1:
                op1index = row
                break

        # for every operation in the second parenthesized list (if there is one)
        for i in range(temp):
            # Find the index of the second operation in the oper_list category table
            op2index = -1
            if oper2:
                for row in range(len(oper_list['id'])):
                    if oper_list.getValue["id"][row] == oper2[i]:
                        op2index = row
                        break

            # prepare the operation matrix
            operation = _bio_unit_prepare_operation(oper_list, op1index, op2index)
            # iterate over every atom in the atom_site reference table
            for r in range(len(atom_site_ref['id'])):

                # if the chain_id of the current atom is not in the chain_id list, skip to the next atom
                if asym_id_list.find(atom_site_ref["label_asym_id"][r]) == -1:
                    continue

                # retrieve the atom's row from the atom_site reference table
                atom = atom_site_ref.ix[r]

                # add this row to the atom_site table for this assembly
                # for s in range(len(attributes) - 1):
                for s in range(len(attributes)):
                    if attributes[s] not in atom_dict:
                        atom_dict[attributes[s]] = []
                    try:
                        atom_dict[attributes[s]].append(atom[attributes[s]])
                    except KeyError:
                        # only excepts on method == 3, because 'bio_uni_counte' is not
                        # found in the original table
                        pass

                # update the atom number and model number for this row
                atom_dict['id'][atomNum - 1] = str(atomNum)

                # method 1: would be to add an identifier to chain_id: e.g. A1, A2, etc.
                if method is 1:
                    asym_id = "%s.%s" % (atom_site_ref["label_asym_id"][r], asym_counter)
                    atom_dict['label_asym_id'][atomNum - 1] = asym_id

                # method 2: incremental counter on the number of models
                elif method is 2:
                    atom_dict['pdbx_PDB_model_num'][atomNum - 1] = modelNum

                # method 3: add a new column to the table with the chain counter
                elif method is 3:
                    if 'bio_unit_counter' not in atom_dict:
                        atom_dict['bio_unit_counter'] = []
                    atom_dict['bio_unit_counter'].append(asym_counter)

                # determine and set the new coordinates
                coords = [atom['Cartn_x'], atom['Cartn_y'], atom['Cartn_z'], 1.0]
                xyz = ['x', 'y', 'z']
                for a in range(3):
                    summ = 0.0
                    for b in range(4):
                        summ += (operation[a][b] * coords[b])
                    atom_dict["Cartn_" + xyz[a]][atomNum - 1] = "%.3f" % summ
                atomNum += 1
            modelNum += 1
            asym_counter += 1

    atom_site = pd.DataFrame(columns=attributes, data=atom_dict)

    return atom_site


def _bio_unit_to_mmcif():
    """
    Method that outputs an mmCIF formatted file with the
    selected biological assembly.

    :return: outputs a file
    """
    # TODO!
    pass


if __name__ == '__main__':
    # X = _mmcif_info_to_dict("tests/CIF/2pah.cif")
    # print(X['exptl_crystal'])
    # print(X['pdbx_struct_assembly'])
    # print(X['pdbx_struct_assembly']['oligomeric_details'])
    # print(X['pdbx_struct_assembly_gen']['asym_id_list'])
    # print(X['pdbx_struct_oper_list']['symmetry_operation'])

    # X = _mmcif_atom("tests/CIF/2pah.cif")
    # X = _mmcif_atom("tests/CIF/3fad.cif")
    # print(X['label_asym_id'], X['Cartn_x'], X['Cartn_y'], X['Cartn_z'])
    # X.to_csv('3fad.tsv', sep='\t')
    # print(X)

    # X = _bio_unit_to_table("tests/CIF/2w4o.cif")
    # X = _bio_unit_to_table("tests/CIF/3fad.cif")
    X = _bio_unit_to_table("tests/CIF/2pah.cif", most_likely=True, method=3)
    print(X['label_asym_id'], X['Cartn_x'], X['Cartn_y'], X['Cartn_z'], X['bio_unit_counter'])
    # X.to_csv('3fad_bio.tsv', sep='\t')
    # print(X)
    pass


def _generic_mmcif_field_to_table(filename, fieldname='_pdbx_poly_seq_scheme.'):
    header = []
    lines = []
    with open(filename) as handle:
        for line in handle:
            if line.startswith(fieldname):
                break
        while line.startswith(fieldname):
            header.append(line.split('.')[1].rstrip())
            line = handle.next()
        while not line.startswith('#'):
            lines.append(line.split())
            line = handle.next()
    return pd.DataFrame(lines, columns=header)