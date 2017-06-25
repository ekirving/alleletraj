#!/usr/bin/env python
# -*- coding: utf-8 -*-

from db_conn import *
from qtldb_api import *
from qtldb_lib import *

VERBOSE = True

dbc = db_conn()
api = qtldb_api()

for species, qtldb_file in QTL_FILES.iteritems():

    # get a list of all the QTLDb IDs and trait names
    data = extract_qtl_fields(qtldb_file, ['QTL_ID', 'Trait_name'])

    # get all the existing traits, indexed by name
    traits = dbc.get_records('traits', key='name')

    # populate all the trait records
    for name in set(data['Trait_name']):
        if name not in traits:
            if VERBOSE:
                print "INFO: Missing trait '%s'" % name

            # get the trait record(s) from the API
            for record in api.get_traits(species, name):
                trait = {
                    'id':    record['@traitID'],
                    'class': record['traitClass'],
                    'type':  record['traitType'],
                    'name':  record['traitName']
                }
                dbc.save_record('traits', trait)

    # get the records
    for record in api.get_qtls(species, data['QTL_ID']):

        # get the trait record
        trait = record.pop('trait')

        print trait

        quit()

