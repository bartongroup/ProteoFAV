#!/usr/bin/env python
# -*- coding: utf-8

"""
ProteoFAV: protein feature aggregation and variants
--------------------------------------------------

Exploring the power of Pandas to work with protein structures,
sequences and genetic variants.

:copyright: (c) 2015-2017.
:license: TBD, see LICENSE for more details.
"""

from .config import defaults
from .main import merge_tables
# from .structures import *
# from .variants import *
# from .analysis import *
# from . import utils
# from . import library


# Set default logging handler to avoid "No handler found" warnings.
import logging
try:  # Python 2.7+
    from logging import NullHandler
except ImportError:
    class NullHandler(logging.Handler):
        def emit(self, record):
            pass

logging.getLogger(__name__).addHandler(NullHandler())


__title__ = 'proteofav'
__version__ = '0.1.0'
__license__ = 'TBD'
__credits__ = [u'Fábio Madeira', u'Thiago Britto-Borges', u'Stuart MacGowan']
