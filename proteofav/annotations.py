#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import absolute_import

import logging

from proteofav.utils import GenericInput

log = logging.getLogger('proteofav')


class Annotations(GenericInput):
    def __init__(self, identifier=None):
        """
        Aggregates res/site-specific Annotations from the UniProt.

        :param identifier: UniProt accession ID
        """

        self.identifier = identifier
        self.table = None

    def fetch(self, idenfitier=None, **kwargs):

        self.identifier = self._get_identifier(idenfitier)
        # TODO
        return self.table

