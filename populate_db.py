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

        # setup the trait record
        trait = dict((field.replace('trait', '').lower(), value)
                     for field, value in record.pop('trait').iteritems() if value is not None)
        trait['name'] = record.pop('name')

        dbc.save_record('traits', trait)

        # setup the pubmed record
        pubmed = api.get_publication(species, record['pubmedID'])
        pubmed['id'] = pubmed.pop('pubmed_ID')
        pubmed['year'] = re.search('\(([0-9]{4})\)', pubmed['authors']).group(1)
        pubmed['journal'] = pubmed['journal']['#text'][:-5]

        dbc.save_record('pubmeds', pubmed)

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

        # drop the source field as the
        record.pop('source', None)

        # filter out any empty values
        qtl = OrderedDict((key, value) for key, value in record.iteritems() if value is not None and value != '-')
        dbc.save_record('qtls', qtl)
