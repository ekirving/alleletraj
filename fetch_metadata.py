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

# details of the GoogleSheet
SHEET_ID = '154wbTEcAUPz4d5v7QLLSwIKTdK8AxXP5US-riCjt2og'
SHEET_TABS = ['Europe and NE Pigs', 'SE Asian Pigs']
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


def fetch_metadata():

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


def populate_samples():

    df = fetch_metadata()

    # open a db connection
    dbc = db_conn()

    bam_files = {}

    # load the BAM file paths
    with open('./pathtopigs.txt', 'r') as fin:
        for file_path in fin:
            # extract the accession code
            code = os.path.basename(file_path).split('_')[0]
            bam_files[code] = file_path.strip()

    # check all the samples for coverage in this interval
    for accession, data in df.iterrows():
        sample = dict()
        sample['accession'] = accession
        for field, value in data.iteritems():
            sample[SHEET_COLS[field]] = value if value != 'NA' else None

        if accession in bam_files:
            sample['path'] = bam_files[accession]

        dbc.save_record('samples', sample)

populate_samples()