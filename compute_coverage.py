#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import httplib2
import pandas as pd

# custom module
import google_sheets as gs

# details of the GoogleSheet
SHEET_ID = '154wbTEcAUPz4d5v7QLLSwIKTdK8AxXP5US-riCjt2og'
SHEET_TABS = ['Europe and NE Pigs', 'SE Asian Pigs']
SHEET_COLS = ['Extract No.', 'Mapped Reads', '% Mapped', '% Mapped-Q30', 'Age', 'Location', 'Country',
              'Wild/Dom Status']
SHEET_NUMERIC = ['Mapped Reads', '% Mapped', '% Mapped-Q30']

MIN_READS = 1000
MIN_MAPQ30 = 0.1

# load the BAM file paths
bam_files = {}

with open('./pathtopigs.txt', 'r') as fin:
    for file_path in fin:
        # extract the accession code
        code = os.path.basename(file_path).split('_')[0]
        bam_files[code] = file_path.strip()

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
    frames.append(df[SHEET_COLS])

# merge the data frames
df = pd.concat(frames)

# assign row names
df = df.set_index(SHEET_COLS[0])

print "%s total samples" % df.shape[0]

# get the records we have BAM files for
df = df.loc[bam_files.keys()]

print "%s samples with BAM files" % df.shape[0]

# replace whitespace
df.replace('', pd.np.nan, inplace=True)

# count the missing data
for col in SHEET_COLS[1:]:
    missing = df[pd.isnull(df[col])].shape[0]
    if missing:
        print " %s missing records from column '%s'" % (missing, col)

# drop all rows with missing data
df.dropna(inplace=True)
df = df[df['% Mapped'] != 'NA']

print "%s samples with complete metadata" % df.shape[0]

# convert stings into real numbers
df[SHEET_NUMERIC] = df[SHEET_NUMERIC].apply(pd.to_numeric)

df = df[df['Mapped Reads'] >= MIN_READS]
df = df[df['% Mapped-Q30'] >= MIN_MAPQ30]

print "%s samples pass basic quality filters (reads > %s, mapQ30 > %s)" % (df.shape[0], MIN_READS, MIN_MAPQ30)

