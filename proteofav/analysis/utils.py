#!/usr/bin/env python
# -*- coding: utf-8 -*-


import colorsys
import os
import numpy as np
from os import path

from ..config import defaults


def is_valid_file(parser, arg):
    """
    Check if arg is a valid file and throw a parsing error if not.

    :param parser: argparse.ArgumentParser
    :param arg: argument
    :return: Open file handle
    !FIXME argparse support file as an type
    https://docs.python.org/2/library/argparse.html#type
    """
    try:
        return open(arg, 'r')
    except:
        parser.error("Not a valid file: %s" % arg)


def delete_file(filename):
    """

    :param filename: File to delete
    :return: None
    """
    os.remove(filename)


def create_directory(directory):
    """
    Creates a directory structure if it does not exist.

    :param directory: directory name (expects full path)
    :return: creates a directory if it does not exist yet
    """
    return os.makedirs(directory)


def _get_colors(num_colors):
    """
    See
    http://stackoverflow.com/questions/470690/how-to-automatically-generate-n-distinct-colors

    :param num_colors: number of color
    :return: a list of colors
    :rtype: list
    """
    colors = []
    for i in np.arange(0., 360., 360. / num_colors):
        hue = i / 360.
        lightness = (50 + np.random.rand() * 10) / 100.
        saturation = (90 + np.random.rand() * 10) / 100.
        colors.append(colorsys.hls_to_rgb(hue, lightness, saturation))
    return colors


def _fractional_to_cartesian(coords, pdb_id, matrix_only=False):
    """
    Converts fractional unit cell coords to cartesian.

    :param coords: Atom coordinates
    :param pdb_id: PDB id
    :param matrix_only: boolean
    :return: cartesian coordinates
    :rtype: np.array
    """

    # retrieve and parse unit cell parameters
    d = _mmcif_unit_cell(pdb_id)
    a2r = np.pi / 180.
    alpha = a2r * float(d['angle_alpha'])
    beta = a2r * float(d['angle_beta'])
    gamma = a2r * float(d['angle_gamma'])
    a = float(d['length_a'])
    b = float(d['length_b'])
    c = float(d['length_c'])

    # unit cell volume
    v = np.sqrt(1 - np.cos(alpha) * np.cos(alpha) - np.cos(beta) * np.cos(beta) -
                np.cos(gamma) * np.cos(gamma) + 2 * np.cos(alpha) * np.cos(beta) * np.cos(gamma))

    # build the transformation matrix
    tr = np.matrix([
        [a, b * np.cos(gamma), c * np.cos(beta)],
        [0, b * np.sin(gamma),
         c * ((np.cos(alpha) - np.cos(beta) * np.cos(gamma)) / np.sin(gamma))],
        [0, 0, c * (v / np.sin(gamma))]
    ])

    if matrix_only:
        return tr

    # now apply the transformation to the coordinates
    # TODO: type check this?
    coords = np.matrix(coords)
    coords = coords * tr

    # return the Nx3 results
    return np.array(coords)


def _autoscale_axes(xyz, margin=5):
    """
    Auto scales the xyz axes by a margin.

    :param xyz: Atom coordinates
    :param margin: spacial margin size/length
    :return: array
    """

    x = xyz[:, 0]
    y = xyz[:, 1]
    z = xyz[:, 2]

    rx = max(x) - min(x)
    ry = max(y) - min(y)
    rz = max(z) - min(z)

    max_range = max([rx, ry, rz]) + margin

    pad = []
    for i in [rx, ry, rz]:
        pad.append((max_range - i) / 2)

    return [[max(x) + pad[0], min(x) - pad[0]],
            [max(y) + pad[1], min(y) - pad[1]],
            [max(z) + pad[2], min(z) - pad[2]]]


def _mmcif_unit_cell(pdb_id):
    """
    Loader of mmCIF unit cell parameters.

    :param pdb_id: PDB id
    :return: pandas table dataframe
    :rtype: dict
    """

    cif_path = path.join(defaults.db_mmcif, pdb_id + '.cif')

    lines = []
    with open(cif_path) as inlines:
        for line in inlines:
            if line.startswith("_cell."):
                lines.append(line.split('.')[1].rstrip())

    l = [i.split() for i in lines]
    d = {k: v for k, v in l}

    return d