#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import httplib2
import pandas as pd

from collections import OrderedDict, defaultdict

# custom module
import google_sheets as gs

from db_conn import db_conn

# TODO make global variable
VERBOSE = True

# Pigs_allTo20042016_shared / old Google sheet
# SHEET_ID = '154wbTEcAUPz4d5v7QLLSwIKTdK8AxXP5US-riCjt2og'
# SHEET_TABS = ['Europe and NE Pigs']  #, 'SE Asian Pigs']

SHEET = {
    # Pig_Table_Final_05_03_18 / new Google sheet
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

SHEET_NA = ['n/a', 'NA', '-', '?']

# list of permissible countries in Europe
EUROPE = ['Belgium', 'Bulgaria', 'Croatia', 'Czech Rep.', 'Denmark', 'England', 'Estonia', 'Faroes', 'France',
          'Germany', 'Greece', 'Hungary', 'Iceland', 'Italy', 'Macedonia (FYROM)', 'Moldova', 'Netherlands', 'Poland',
          'Portugal', 'Romania', 'Serbia', 'Slovakia', 'Spain', 'Sweden', 'Switzerland', 'Ukraine']

# make sure we can properly inspect the data if we want to
pd.set_option('max_colwidth', 1000)


def fetch_metadata(sheet_id, sheet_tabs, sheet_columns):

    # connect to GoogleSheets
    credentials = gs.get_credentials()
    http = credentials.authorize(httplib2.Http())
    discoveryUrl = 'https://sheets.googleapis.com/$discovery/rest?version=v4'
    service = gs.discovery.build('sheets', 'v4', http=http, discoveryServiceUrl=discoveryUrl)

    # fetch the sample metadata
    frames = []

    for tab in sheet_tabs:
        # get the tab of data from the GoogleSheet
        values = service.spreadsheets().values().get(spreadsheetId=sheet_id, range=tab).execute().get('values', [])

        # convert it into a data frame
        df = pd.DataFrame.from_records(values[1:], columns=values[0])

        # add to the list
        frames.append(df[sheet_columns.keys()])

    # merge the data frames
    df = pd.concat(frames)

    # assign row names
    df = df.set_index(df['Extract No.'])

    return df


def confirm_age_mapping(species):
    """
    Make sure that all the free-text dates have proper numeric mappings
    """

    dbc = db_conn()

    print "INFO: Confirming {} age mappings".format(species)

    missing = dbc.get_records_sql("""
        SELECT s.*
          FROM samples s
     LEFT JOIN sample_dates sd
            ON s.age <=> sd.age
         WHERE s.species = '{species}'
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
           AND country IN ('{europe}')""".format(species=species, europe=europe))


def populate_samples(species):

    # open a db connection
    dbc = db_conn()

    if species != 'pig':
        # TODO make this work for all species not just pigs
        raise Exception('Not implemented yet for %' % species)

    # fetch all the samples from the GoogleDoc spreadsheet
    df = fetch_metadata(SHEET[species]['id'], SHEET[species]['tab'], SHEET[species]['cols'])

    if VERBOSE:
        print "INFO: Updating %s %s samples" % (len(df), species)

    bam_files = defaultdict(list)

    # load the BAM file paths
    with open('./data/bam_files_{}.txt'.format(species), 'r') as fin:
        for file_path in fin:
            # extract the accession code
            accession = os.path.basename(file_path).replace('_rmdup.bam', '').strip()
            bam_files[accession].append(file_path.strip())

    for accession, data in df.iterrows():
        sample = dbc.get_record('samples', {'accession': accession}) or dict()
        sample['species'] = species
        sample['accession'] = accession
        for field, value in data.iteritems():
            sample[SHEET[species]['cols'][field]] = value if value not in SHEET_NA else None

        dbc.save_record('samples', sample)

        # save the BAM file paths
        for path in bam_files[accession]:
            bam_file = dict()
            bam_file['sample_id'] = sample['id']  # TODO not set on first pass
            bam_file['path'] = path
            dbc.save_record('sample_files', bam_file)

    # make sure that all the free-text dates have proper numeric mappings
    confirm_age_mapping(species)

    # samples are valid if they are from Europe and have a BAM file
    mark_valid_samples(species)

    if VERBOSE:
        print "INFO: Finished updating %s %s samples" % (len(df), species)