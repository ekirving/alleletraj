#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import httplib2
import pandas as pd

from pprint import pprint

from collections import OrderedDict

# custom module
import google_sheets as gs

from db_conn import db_conn

# TODO make global variable
VERBOSE = True

# details of the GoogleSheet
SHEET_ID = '154wbTEcAUPz4d5v7QLLSwIKTdK8AxXP5US-riCjt2og'
SHEET_TABS = ['Europe and NE Pigs']  #, 'SE Asian Pigs']
SHEET_COLS = OrderedDict([
    ('Extract No.',      'accession'),
    ('Mapped Reads',     'map_reads'),
    ('% Mapped',         'map_prcnt'),
    ('% Mapped-Q30',     'map_prcnt_q30'),
    ('Age',              'age'),
    ('Location',         'location'),
    ('Country',          'country'),
    ('Wild/Dom Status',  'status')
])

# make sure we can properly inspect the data if we want to
pd.set_option('max_colwidth', 1000)


def fetch_metadata(species):

    # connect to GoogleSheets
    credentials = gs.get_credentials()
    http = credentials.authorize(httplib2.Http())
    discoveryUrl = 'https://sheets.googleapis.com/$discovery/rest?version=v4'
    service = gs.discovery.build('sheets', 'v4', http=http, discoveryServiceUrl=discoveryUrl)

    # fetch the sample metadata
    frames = []

    for tab in SHEET_TABS:
        # get the tab of data from the GoogleSheet
        values = service.spreadsheets().values().get(spreadsheetId=SHEET_ID, range=tab).execute().get('values', [])

        # convert it into a data frame
        df = pd.DataFrame.from_records(values[1:], columns=values[0])

        # add to the list
        frames.append(df[SHEET_COLS.keys()])

    # merge the data frames
    df = pd.concat(frames)

    # assign row names
    df = df.set_index(df['Extract No.'])

    return df


def populate_samples(species):

    if species != 'pig':
        # TODO make this work for all species not just pigs
        raise Exception('Not implemented yet for %' % species)

    # open a db connection
    dbc = db_conn()

    # fetch all the samples from the GoogleDoc spreadsheet
    df = fetch_metadata(species)

    if VERBOSE:
        print "INFO: Found %s %s samples to update" % (len(df), species)

    bam_files = {}

    # load the BAM file paths
    with open('./data/bam_files_{}.txt'.format(species), 'r') as fin:
        for file_path in fin:
            # extract the accession code
            code = os.path.basename(file_path).split('_')[0]
            bam_files[code] = file_path.strip()

    for accession, data in df.iterrows():
        sample = dict()
        sample['species'] = species
        sample['accession'] = accession
        for field, value in data.iteritems():
            sample[SHEET_COLS[field]] = value if value != 'NA' else None

        sample['path'] = bam_files[accession] if accession in bam_files else None

        dbc.save_record('samples', sample)

    if VERBOSE:
        print "INFO: Finished updating %s %s samples" % (len(df), species)