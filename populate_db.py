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

    if VERBOSE:
        print "INFO: Found %s new QTLs to add" % len(new_ids)

    # get all the new records
    for record in api.get_qtls(species, new_ids):

        # ignore placeholder values
        for key, value in record.iteritems():
            record[key] = None if value == '-' else value

        # setup the trait record
        trait = dict((field.replace('trait', '').lower(), value)
                     for field, value in record.pop('trait').iteritems() if value is not None)
        trait['name'] = record.pop('name')

        try:
            dbc.save_record('traits', trait)

        except Exception as e:
            print "qtlID = %s" % record['id']
            for key, value in trait.iteritems():
                print (key, value)
            raise e

        # set the trait ID
        record['traitID'] = trait['id']

        # flatten the other nested records
        for field, value in record.iteritems():

            if type(value) is OrderedDict:
                nested = record.pop(field)

                for child_name, child_value in nested.iteritems():
                    # for doubly nested fields, use the parent name as a prefix
                    if type(child_value) is OrderedDict:
                        for key in child_value:
                            record[child_name + '_' + key] = child_value[key]
                    elif field in ['gene']:
                        record[field + child_name.title()] = child_value
                    else:
                        record[child_name] = child_value

        # filter out any empty values
        qtl = OrderedDict((key, value) for key, value in record.iteritems() if value is not None)

        try:
            dbc.save_record('qtls', qtl)

        except Exception as e:
            for key, value in qtl.iteritems():
                print (key, value)
            raise e

        print 'Added..'