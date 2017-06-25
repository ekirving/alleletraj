#!/usr/bin/env python
# -*- coding: utf-8 -*-

from db_conn import db_conn
from qtldb_api import qtldb_api

# make the db connection
dbc = db_conn()

qtldb = qtldb_api()

qtl_dump = {
    'pig': 'mapDwnLd46065CEZR.txt'
}

# TODO remove me
itr = 0


def extract_qtl_ids(dbfile):

    ids = []

    # get all the QTL IDs for this scope (i.e. species)
    with open(dbfile, 'rU') as fin:
        fin.readline()
        for line in fin:
            try:
                ids.append(line.split()[0])
            except IndexError:
                # ignore badly formatted lines
                pass

    return ids


def populate_traits(species):

    # get all the existing traits
    existing = dbc.get_records('traits')
    print existing
    # get all the trait types
    for record in qtldb.get_traits(species):
        traitid = int(record['@traitID'])

        if traitid not in existing:
            trait = {
                'id':    record['@traitID'],
                'class': record['traitClass'],
                'type':  record['traitType'],
                'name':  record['traitName']
            }

            dbc.save_record('traits', trait)


for species, dbfile in qtl_dump.iteritems():

    # get a list of all the QTLDb IDs
    ids = extract_qtl_ids(dbfile)

    populate_traits(species)

    quit()

    # get the records
    for record in qtldb.get_qtls(species, ids):

        # get the trait record
        trait = record.pop('trait')

        print trait

        itr += 1

        if itr > 10:
            quit()
        # http://www.animalgenome.org/cgi-bin/QTLdb/API/ifetch?q=1,15&s=pig


