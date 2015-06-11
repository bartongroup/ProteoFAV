#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
---------
config.py
---------

This defines the methods that load and validate user defined
parameters.
Usage
defaults = get_config("api_rcsb")
r = request.get(url = defaults.api_rcsb)

:module_version: 1.0
:created_on: 26-02-2015

"""
from ConfigParser import ConfigParser

__all__ = ["get_config"]

config = ConfigParser()
config_default = "../config.txt"
config.read(config_default)


class Defaults(object):
    pass

def get_config(*vars, **default):
    """
    Gets the config values defined locally.

    :param input_config: input config file
    :param var: list of [str, ]
    :param default:
    :return: returns an object
    """

    if not default:
        default = Defaults()
    for var in vars:
        found = False
        for section in config.sections():
            for var_name, var_par in config.items(section):
                if var == var_name:
                    if var_par == "...":
                        raise NotImplementedError(  # TODO offer to fill fields
                            "Fill {} parameter in config.txt".format(var_par))
                    else:
                        setattr(default, var, var_par)
                    found = True
        if not found:
            raise TypeError

    return default
