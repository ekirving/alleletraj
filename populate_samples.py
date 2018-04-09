#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import httplib2

from collections import OrderedDict, defaultdict

from pprint import pprint

# custom module
import google_sheets as gs

from db_conn import db_conn

# TODO make global variable
VERBOSE = True

# Pigs_allTo20042016_shared / old Google sheet
# SHEET_ID = '154wbTEcAUPz4d5v7QLLSwIKTdK8AxXP5US-riCjt2og'
# SHEET_TABS = ['Europe and NE Pigs']  #, 'SE Asian Pigs']

GOOGLE_SHEET = {

    # Pig_Table_Final_05_03_18
    'pig': {
        'id':   '1IWCt8OtTz6USOmN5DO0jcYxZOLnnOVdstTGzRcBZolI',
        'tabs': ['Everything for the paper - updated'],
        'cols': OrderedDict([
                    ('Extract No.',       'accession'),
                    ('Total Reads',       'map_reads'),
                    ('% Mapped',          'map_prcnt'),
                    ('Age',               'age'),
                    ('Period',            'period'),
                    ('Neolithic (Y/N/W)', 'neolithic'),
                    ('Location',          'location'),
                    ('Country',           'country'),
                    ('Wild/Dom Status',   'status'),
                    ('GMM Status',        'gmm_status')
                ])
    }
}

# list of junk input to mask with NULL
SHEET_NA = ['n/a', 'NA', '-', '?', 'NULL', 'None', '...', '']

AGE_MAP = {

    'pig': {
        'id': '1bH5u_qDaFXJdTyybeahqgF7je17td0FdyOMs_tlECdA',
        'tabs': ['Age Map'],
        'cols': OrderedDict([
            ('Age',         'age'),
            ('Confident',   'confident'),
            ('Lower (BP)',  'lower'),
            ('Upper (BP)',  'upper'),
            ('Median (BP)', 'median'),
        ])
    }
}

C14_SHEET = {

    'pig': {
        'id': '1odoL9hQh87bLLe3yipbo-CKKXLvIgb5n_kfoqSALHi8',
        'tabs': ['All Dates'],
        'cols': OrderedDict([
            ('Extract_No',              'accession'),
            ('From Cal BP (Int Cal13)', 'lower'),
            ('To Cal BP',               'upper'),
        ])
    }
#
}


# list of permissible countries in Europe
EUROPE = ['Belgium', 'Bulgaria', 'Croatia', 'Czech Rep.', 'Denmark', 'England', 'Estonia', 'Faroes', 'France',
          'Germany', 'Greece', 'Hungary', 'Iceland', 'Italy', 'Macedonia (FYROM)', 'Moldova', 'Netherlands', 'Poland',
          'Portugal', 'Romania', 'Serbia', 'Slovakia', 'Spain', 'Sweden', 'Switzerland', 'Ukraine']


def fetch_google_sheet(sheet_id, sheet_tabs, sheet_columns):
    """
    Fetch a given list of tabs and columns from a Google Sheet.
    """

    # connect to GoogleSheets
    credentials = gs.get_credentials()
    http = credentials.authorize(httplib2.Http())
    discovery_url = 'https://sheets.googleapis.com/$discovery/rest?version=v4'
    service = gs.discovery.build('sheets', 'v4', http=http, discoveryServiceUrl=discovery_url)

    records = []

    # convert junk values into None
    mask_null = lambda x: x if x not in SHEET_NA else None

    for tab in sheet_tabs:
        # get the tab of data from the GoogleSheet
        values = service.spreadsheets().values().get(spreadsheetId=sheet_id, range=tab).execute().get('values', [])

        fields = values[0]

        # fetch the requested columns
        for row in values[1:]:
            records.append(dict([(sheet_columns[fields[idx]], mask_null(value))
                                 for idx, value in enumerate(row) if fields[idx] in sheet_columns]))

    return records


def sync_c14_dates(species):
    """
    Fetch all the C14 dates
    """

    dbc = db_conn()

    print "INFO: Synchronising {} C14 dates".format(species)

    # get the google sheet details
    sheet = C14_SHEET[species]

    # fetch all the manual age mappings
    records = fetch_google_sheet(sheet['id'], sheet['tabs'], sheet['cols'])

    # update the `sample_dates` table
    for record in records:
        if record['accession'] is not None:
            record['confident'] = 'Yes'
            dbc.save_record('sample_dates_c14', record)


def confirm_age_mapping(species):
    """
    Make sure that all the free-text dates have proper numeric mappings
    """

    dbc = db_conn()

    print "INFO: Confirming {} age mappings".format(species)

    # get the google sheet details
    sheet = AGE_MAP[species]

    # fetch all the manual age mappings
    records = fetch_google_sheet(sheet['id'], sheet['tabs'], sheet['cols'])

    # update the `sample_dates` table
    for record in records:
        if record['age'] is not None:
            dbc.save_record('sample_dates', record)

    # check if any samples have an age which is unmapped
    missing = dbc.get_records_sql("""
        SELECT s.*
          FROM samples s
     LEFT JOIN sample_dates sd
            ON s.age = sd.age
         WHERE s.species = '{species}'
           AND s.age IS NOT NULL
           AND sd.id IS NULL""".format(species=species))

    if missing:
        print "ERROR: Not all sample ages have numeric mappings!"

        for id, sample in missing.iteritems():
            print "{} - '{}'".format(sample['accession'], sample['age'])

        quit()


def mark_valid_samples(species):
    """
    Samples are valid if they are from Europe and have a BAM file.
    """

    dbc = db_conn()

    print "INFO: Marking valid {} samples".format(species)

    # make a list of permissible countries
    europe = "','".join(EUROPE)

    dbc.execute_sql("""
        UPDATE samples
          JOIN sample_files
            ON sample_files.sample_id = samples.id
           SET valid = 1
         WHERE species = '{species}'
           AND country IN ('{europe}')
           AND COALESCE(neolithic, '') != 'Y'""".format(species=species, europe=europe))


def populate_samples(species):

    # open a db connection
    dbc = db_conn()

    if species != 'pig':
        # TODO make this work for all species not just pigs
        raise Exception('Not implemented yet for %' % species)

    # get the google sheet details
    sheet = GOOGLE_SHEET[species]

    # fetch all the samples from the GoogleDoc spreadsheet
    samples = fetch_google_sheet(sheet['id'], sheet['tabs'], sheet['cols'])

    if VERBOSE:
        print "INFO: Updating %s %s samples" % (len(samples), species)

    bam_files = defaultdict(list)

    # load the BAM file paths
    with open('./data/bam_files_{}.txt'.format(species), 'r') as fin:
        for line in fin:
            # extract the accession code
            accession = os.path.basename(line).replace('_rmdup.bam', '').strip()
            bam_files[accession].append(line.strip())

    for sample in samples:
        accession = sample['accession']

        # save the sample record
        sample['species'] = species
        dbc.save_record('samples', sample)

        # fetch the sample ID
        sample = dbc.get_record('samples', {'accession': accession})

        # save the BAM file paths
        for path in bam_files[accession]:
            bam_file = dict()
            bam_file['sample_id'] = sample['id']
            bam_file['path'] = path

            if not dbc.exists_record('sample_files', bam_file):
                dbc.save_record('sample_files', bam_file)

    # fetch all the C14 dates
    sync_c14_dates(species)

    # make sure that all the free-text dates have proper numeric mappings
    confirm_age_mapping(species)

    # samples are valid if they are from Europe and have a BAM file
    mark_valid_samples(species)

    if VERBOSE:
        print "INFO: Finished updating %s %s samples" % (len(samples), species)