#!/usr/bin/env python
# -*- coding: utf-8 -*-

from db_conn import *
from qtldb_api import *
from qtldb_lib import *

# TODO make global variable
VERBOSE = True


def populate_qtls(species):
    """
    Fetch all the QTLs from the QTLdb API and populate the local database.
    """
    dbc = db_conn()
    api = qtldb_api()

    # get the file containing all the QTL IDs
    qtldb_file = QTL_FILES[species]

    # get a list of all the QTLDb IDs from the tsv dump
    data = extract_qtl_fields(qtldb_file, ['QTL_ID'])

    if VERBOSE:
        print "INFO: Processing %s %s QTLs from '%s'" % (len(data['QTL_ID']), species, qtldb_file)

    # convert all the IDs to int
    qtl_ids = [int(qtl_id) for qtl_id in data['QTL_ID'] if qtl_id.isdigit()]

    # get all the QTLs already in the DB
    qtls = dbc.get_records('qtls', key='id')

    # find the new IDs in the list
    new_ids = list(set(qtl_ids) - set(qtls.keys()))

    if VERBOSE:
        print "INFO: Found %s new %s QTLs to add" % (len(new_ids), species)

    added = 0

    # get all the new records
    for record in api.get_qtls(species, new_ids):

        # TODO when resultset is len() = 1 then this throws an error
        # extract the nest trait record
        trait = record.pop('trait')
        trait['name'] = record.pop('name')

        # set the tait foreign key on the main record
        record['traitID'] = trait['traitID']

        # does the trait exist
        if not dbc.exists_record('traits', {'id': record['traitID']}):

            # setup the trait record
            trait = dict((field.replace('trait', '').lower(), value) for field, value in trait.iteritems())
            trait['type'] = api.get_trait_type(species, trait['id'], trait['name'])

            dbc.save_record('traits', trait)

        # does the publication exist
        if not dbc.exists_record('pubmeds', {'id': record['pubmedID']}):

            # setup the pubmed record
            pubmed = api.get_publication(species, record['pubmedID'])

            if pubmed:
                pubmed['id'] = pubmed.pop('pubmed_ID')
                pubmed['year'] = re.search('\(([0-9]{4})\)', pubmed['authors']).group(1)
                pubmed['journal'] = pubmed['journal']['#text'][:-5]

                dbc.save_record('pubmeds', pubmed)
            else:
                # some records have a bogus pubmedID
                record['pubmedID'] = None

        # flatten the other nested records
        for field, value in record.iteritems():

            if type(value) is OrderedDict:
                nested = record.pop(field)

                for nested_name, nested_value in nested.iteritems():
                    # for doubly nested fields, use the parent name as a prefix
                    if type(nested_value) is OrderedDict:
                        for key in nested_value:
                            record[nested_name + '_' + key] = nested_value[key]
                    elif field in ['gene']:
                        record[field + nested_name.title()] = nested_value
                    else:
                        record[nested_name] = nested_value

        # drop any lingering malformed fields
        record.pop('source', None)
        record.pop('breeds', None)
        record.pop('effects', None)
        record.pop('statTests', None)

        # handle malformed data
        for field in ['linkageLoc_end', 'linkageLoc_peak', 'linkageLoc_start']:
            if field in record and record[field] is not None:
                record[field] = re.sub('[^0-9.]', '', record[field])

        # filter out any empty values
        qtl = OrderedDict((key, value) for key, value in record.iteritems() if value != '-')

        qtl['species'] = species

        dbc.save_record('qtls', qtl)

        added += 1

        if VERBOSE and added % 100 == 0:
            print "INFO: Added %5d new QTLs" % added

    if VERBOSE:
        print "INFO: Finished adding %s new %s QTLs" % (len(new_ids), species)
