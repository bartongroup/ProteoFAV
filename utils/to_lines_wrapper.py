#!/usr/bin/env python
# -*- coding: utf-8

"""
Created on 18:18 03/06/15 2015 

"""
import logging

import requests

logger = logging.getLogger(__name__)


def to_lines_wrapper(function, url=None, filename=None, fileobj=None):
    """
    wrong
    :param function:
    :param url:
    :param filename: path to file
    :param fileobj: python file objectio
    :return:
    """

    if filename:
        with open(filename) as open_f:
            return function(open_f)

    elif fileobj:
        if hasattr(fileobj, "closed"):
            if not fileobj.closed:
                with fileobj as open_f:
                    return function(open_f)

    elif url:
        r = requests.get(url, stream=True)
        if r.ok:
            return r.iter_lines()

    else:
        err = "{} uses one from the followin parameters: {}".format
        raise ValueError(err(function.__name__,
                             " ".join("filename", "fileobj", "url")) )

