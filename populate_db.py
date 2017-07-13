#!/usr/bin/env python
# -*- coding: utf-8 -*-

from db_conn import *
from qtldb_api import *
from qtldb_lib import *

import pprint

VERBOSE = True

dbc = db_conn()
api = qtldb_api()

for species, qtldb_file in QTL_FILES.iteritems():

    # get a list of all the QTLDb IDs and trait names
    data = extract_qtl_fields(qtldb_file, ['QTL_ID', 'Trait_name'])

    # # get all the existing traits, indexed by name
    # traits = dbc.get_records('traits', key='name')
    #
    # # populate all the trait records
    # for name in set(data['Trait_name']):
    #
    #     if name not in traits:
    #
    #         if VERBOSE:
    #             print "INFO: Missing trait '%s'" % name
    #
    #         # get the trait record(s) from the API
    #         for record in api.get_traits(species, name):
    #             trait = {
    #                 'id':    record['@traitID'],
    #                 'class': record['traitClass'],
    #                 'type':  record['traitType'],
    #                 'name':  record['traitName']
    #             }
    #
    #             dbc.save_record('traits', trait)

    # get all the QTLs
    qtls = dbc.get_records('qtls', key='id')

    # convert all the IDs to int
    qtl_ids = [int(qtl_id) for qtl_id in data['QTL_ID'] if qtl_id.isdigit()]

    # get the new IDs
    new_ids = list(set(qtl_ids) - set(qtls.keys()))

    # FIXME...
    schema = defaultdict(set)
    ctr = 0

    # get all the new records
    for record in api.get_qtls(species, new_ids):

        # pprint.pprint(dict(record))

        # insert the nested records
        for field, value in record.iteritems():

            if type(value) is OrderedDict:
                value = record.pop(field)
                # pprint.pprint(dict(value))
                schema[field] |= set(value.keys())

        # pprint.pprint(dict(record))
        schema['qtl'] |= set(record.keys())

        ctr += 1

        if ctr > 500:
            break

    pprint.pprint(dict(schema))
