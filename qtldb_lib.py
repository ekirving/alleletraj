#!/usr/bin/env python
# -*- coding: utf-8 -*-

from collections import defaultdict


def extract_qtl_fields(dbfile, fields):

    data = defaultdict(list)

    # get all the QTL IDs for this scope (i.e. species)
    with open(dbfile, 'rU') as fin:
        # get the column headers
        header = fin.readline().split('\t')

        # get the index of the target fields
        columns = [(idx, field) for idx, field in enumerate(header) if field in fields]

        for line in fin:
            try:
                line = line.split('\t')
                for idx, field in columns:
                    datum = line[idx].strip()
                    if datum:
                        data[field].append(datum)
            except IndexError:
                # ignore badly formatted lines
                pass

    return data

