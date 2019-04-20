#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function

import httplib2
import os
import google_sheets as gs
from collections import defaultdict, OrderedDict
from database import Database

# from pipeline_utils import *

GOOGLE_SHEET = {

    # Pig_Table_Final_05_03_18
    # 'pig': {
    #     'id':   '1IWCt8OtTz6USOmN5DO0jcYxZOLnnOVdstTGzRcBZolI',
    #     'tabs': ['Everything for the paper - updated'],
    #     'cols': OrderedDict([
    #                 ('Extract No.',       'accession'),
    #                 ('Total Reads',       'map_reads'),
    #                 ('% Mapped',          'map_prcnt'),
    #                 ('Age',               'age'),
    #                 ('Period',            'period'),
    #                 ('Location',          'location'),
    #                 ('Country',           'country'),
    #                 ('Wild/Dom Status',   'status'),
    #                 ('GMM Status',        'gmm_status'),
    #                 ('Group',             'group'),
    #                 ('Haplogroup',        'haplogroup'),
    #                 ('DNA',               'dna')
    #             ])
    # },

    'pig': {
        'id': '1GBxNiRWAqPdz4MdSpi0ec_K8x4TRns31VgUcICq68qo',
        'tabs': ['final combined'],
        'cols': OrderedDict([
                    ('Extract No. / Lab code', 'accession'),
                    ('Total Reads', 'map_reads'),
                    ('% Mapped', 'map_prcnt'),
                    ('Age', 'age'),
                    ('Age (Mean years BP)', 'age_int'),
                    ('Period', 'period'),
                    ('Location', 'location'),
                    ('Country', 'country'),
                    ('Final status (MC1R+Morpho+Context)', 'status'),
                    ('Genotype MC1R', 'mc1r_snp')
                ])
    },

    # HorseSelection_LO4EIP-TRANSFERED
    'horse': {
        'id':   '1BMvIwYj-d8t3mpf67rzabrEvDoB8hBZbyS6XfGwcwUU',
        'tabs': ['Ancient'],
        'cols': OrderedDict([
                    ('Name',     'accession'),
                    ('Status',   'status'),
                    ('path',     'path'),
                    ('Age BP',   'age'),
                    ('Age',      'period'),
                    ('Site',     'location'),
                ])
    }
}

# list of junk input to mask with NULL
SHEET_NA = ['n/a', 'NA', 'N', '-', '?', 'NULL', 'None', '...', '']

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

RADIOCARBON_SHEET = {

    'pig': {
        'id': '1odoL9hQh87bLLe3yipbo-CKKXLvIgb5n_kfoqSALHi8',
        'tabs': ['All Dates'],
        'cols': OrderedDict([
            ('Extract_No',              'accession'),
            ('From Cal BP (Int Cal13)', 'lower'),
            ('To Cal BP',               'upper'),
        ])
    }
}


# list of permissible countries in Europe
EUROPE = [
    'Austria', 'Belgium', 'Bosnia-Herzegovina', 'Bulgaria', 'Crimea, Ukraine', 'Croatia', 'Czech Rep.', 'Denmark',
    'England', 'Estonia', 'Europe', 'Faroe Islands', 'Faroes', 'Finland', 'France', 'Germany', 'Greece', 'Hungary',
    'Iberia', 'Iceland', 'Ireland', 'Italy', 'Macedonia', 'Macedonia (FYROM)', 'Moldova', 'Netherlands', 'Norway',
    'Poland', 'Portugal', 'Portugal/France', 'Romania', 'Russia', 'Sardinia', 'Scotland', 'Serbia', 'Slovakia',
    'Spain', 'Sweden', 'Switzerland', 'UK', 'Ukraine', 'West Caucasus, north slope'
]

# list of non-permissible countries outside of Europe
NON_EUROPE = ['Africa', 'Armenia', 'Azerbaijan', 'Cyprus', 'Egypt', 'Egyptian', 'EuroAm', 'Georgia', 'Iran', 'Iraq',
              'Israel', 'Morocco', 'Sudan', 'Syria', 'Tunisia', 'Turkey', 'Turkmenistan', 'United Arab Emirates']


# the arbitrary +/- age uncertainty for median age dates
MEDIAN_AGE_UNCERT = 100

BIN_WIDTH = 500
BIN_PERCENT = 0.5  # samples must overlap a bin by >= 50%


def unicode_truncate(s, length, encoding='utf-8'):
    encoded = s.encode(encoding)[:length]
    return encoded.decode(encoding, 'ignore')


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
    def mask_null(val):
        return val if val not in SHEET_NA else None

    for tab in sheet_tabs:
        # get the tab of data from the GoogleSheet
        values = service.spreadsheets().values().get(spreadsheetId=sheet_id, range=tab).execute().get('values', [])

        fields = values[0]

        # fetch the requested columns
        for row in values[1:]:
            records.append(dict([(sheet_columns[fields[idx]], mask_null(value))
                                 for idx, value in enumerate(row) if fields[idx] in sheet_columns]))

    return records


def sync_radiocarbon_dates(species):
    """
    Fetch all the radiocarbon dates
    """

    dbc = Database(species)

    print("INFO: Synchronising radiocarbon dates")

    try:
        # get the google sheet details
        sheet = RADIOCARBON_SHEET[species]

    except KeyError:

        print("WARNING: No radiocarbon spreadsheet for {}".format(species))
        return

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

    dbc = Database(species)

    print("INFO: Confirming age mappings")

    try:
        # get the google sheet details
        sheet = AGE_MAP[species]

    except KeyError:

        print("WARNING: No age mapping for {}".format(species))
        return

    # fetch all the manual age mappings
    records = fetch_google_sheet(sheet['id'], sheet['tabs'], sheet['cols'])

    # update the `sample_dates` table
    for record in records:
        if record['age'] is not None:
            dbc.save_record('sample_dates', record)

    # check if any samples have an age which is unmapped
    missing = dbc.get_records_sql("""
        SELECT s.age
          FROM samples s
     LEFT JOIN sample_dates sd
            ON s.age = sd.age
         WHERE s.age IS NOT NULL
           AND sd.id IS NULL
      GROUP BY s.age""", key=None)

    if missing:
        print("ERROR: Not all sample ages have numeric mappings!")

        for age in missing:
            print(age['age'].encode('utf-8'))

        quit()


def confirm_country_mapping(species):
    """
    Make sure that all the free-text countries have been properly mapped to Europe.
    """

    dbc = Database(species)

    print("INFO: Confirming country mappings")

    # make a list of known countries
    world = "','".join(EUROPE + NON_EUROPE)

    missing = dbc.get_records_sql("""
        SELECT DISTINCT s.country
          FROM samples s
         WHERE s.country NOT IN ('{world}')
           """.format(world=world), key=None)

    if missing:
        print("ERROR: Not all countries have been mapped!")

        for country in missing:
            print(country['country'])

        quit()


def mark_valid_pigs(species):
    """
    Pig samples are valid if they are from Europe and have a BAM file or MC1R genotype.
    """
    db_lock_tables = ['samples']

    dbc = Database(species)

    print("INFO: Marking valid samples")

    # make a list of permissible countries
    europe = "','".join(EUROPE)

    dbc.execute_sql("""
        UPDATE samples s
     LEFT JOIN sample_files sf
            ON sf.sample_id = s.id
           SET s.valid = 1
         WHERE s.country IN ('{europe}')
           AND (sf.id IS NOT NULL 
            OR  s.mc1r_snp IS NOT NULL) """.format(europe=europe))


def bin_samples(species):
    """
    Assign samples to temporal bins
    """

    dbc = Database(species)

    print("INFO: Binning samples")

    # sample_bins
    dbc.execute_sql("TRUNCATE TABLE sample_bins")

    # TODO update to just use the median age
    # SQL fragment to get the most precise dates for each sample
    sql_age = """
        SELECT s.id sample_id,
               COALESCE(c14.lower, sd.lower, sd.median + {uncert}) lower,
               COALESCE(c14.upper, sd.upper, sd.median - {uncert}) upper
          FROM samples s
     LEFT JOIN sample_dates sd
            ON s.age = sd.age
     LEFT JOIN sample_dates_c14 c14
            ON c14.accession = s.accession
         WHERE s.valid = 1
           """.format(uncert=MEDIAN_AGE_UNCERT)

    # get the maximum date
    max_lower = dbc.get_records_sql("""
        SELECT MAX(lower) max_lower
          FROM ({age}) AS age
           """.format(age=sql_age), key=None)[0].pop('max_lower')

    # round to nearest multiple of BIN_WIDTH
    bin_start = int(round(float(max_lower)/BIN_WIDTH) * BIN_WIDTH)
    bin_end = 0

    # iterate over each temporal bin
    for bin_lower in range(bin_start, bin_end, -BIN_WIDTH):
        bin_upper = bin_lower - BIN_WIDTH

        # get all the samples which overlap this bin by >= BIN_OVERLAP
        dbc.execute_sql("""
             INSERT IGNORE
               INTO sample_bins (sample_id, bin, overlap, perct_overlap)
             SELECT sample_id,
                    '{binlower} - {binupper}' AS bin,
                    LEAST(lower, {binlower}) - GREATEST(upper, {binupper}) AS overlap,
                    (LEAST(lower, {binlower}) - GREATEST(upper, {binupper})) / (lower - upper) AS perct_overlap
               FROM ({age}) as age
              WHERE lower >= {binupper}
                AND upper <= {binlower}
             HAVING perct_overlap >= {binpercent}
                """.format(age=sql_age, binpercent=BIN_PERCENT, binlower=bin_lower, binupper=bin_upper))


def populate_pig_samples(species):

    # open a db connection
    dbc = Database(species)

    # get the google sheet details
    sheet = GOOGLE_SHEET[species]

    # fetch all the samples from the GoogleDoc spreadsheet
    samples = fetch_google_sheet(sheet['id'], sheet['tabs'], sheet['cols'])

    print("INFO: Updating {} samples".format(len(samples)))

    bam_files = defaultdict(list)

    # load the BAM file paths
    with open('./data/bam_files_{}.txt'.format(species), 'r') as fin:
        for line in fin:
            # extract the accession code
            accession = os.path.basename(line).replace('_rmdup.bam', '').strip()
            bam_files[accession].append(line.strip())

    for sample in samples:
        accession = sample['accession']

        try:
            sample['map_reads'] = int(sample['map_reads'])
        except (ValueError, TypeError, KeyError):
            # handle missing and malformed map_reads info
            sample['map_reads'] = None

        # limit free-text age to 255 chracters
        if sample['age']:
            sample['age'] = unicode_truncate(sample['age'], 255)

        # save the sample record
        dbc.save_record('samples', sample)

        # fetch the sample record (so we can link the BAM files to the ID)
        sample = dbc.get_record('samples', {'accession': accession})

        # save the BAM file paths
        for path in bam_files[accession]:
            bam_file = dict()
            bam_file['sample_id'] = sample['id']
            bam_file['path'] = path

            if not dbc.exists_record('sample_files', bam_file):
                dbc.save_record('sample_files', bam_file)

    # fetch all the C14 dates
    sync_radiocarbon_dates(species)

    # make sure that all the free-text dates have proper numeric mappings
    confirm_age_mapping(species)

    # make sure that all the free-text countries have been properly mapped to Europe
    confirm_country_mapping(species)

    # samples are valid if they are from Europe and have a BAM file
    mark_valid_pigs(species)

    # assign samples to temporal bins
    bin_samples(species)

    print("INFO: Finished updating {} samples".format(len(samples)))


def populate_horse_samples(species):

    # open a db connection
    dbc = Database(species)

    # get the google sheet details
    sheet = GOOGLE_SHEET[species]

    # fetch all the samples from the GoogleDoc spreadsheet
    samples = fetch_google_sheet(sheet['id'], sheet['tabs'], sheet['cols'])

    print("INFO: Updating {} samples".format(len(samples)))

    for sample in samples:
        accession = sample['accession']

        # get the file path
        path = sample.pop('path', None)

        # all horses are valid
        sample['valid'] = 1

        # save the sample record
        dbc.save_record('samples', sample)

        # fetch the sample record (so we can link the BAM files to the ID)
        sample = dbc.get_record('samples', {'accession': accession})

        bam_file = dict()
        bam_file['sample_id'] = sample['id']
        bam_file['path'] = '/home/ludo/inbox/BAMs/ancient/' + os.path.basename(path)

        if not dbc.exists_record('sample_files', bam_file):
            dbc.save_record('sample_files', bam_file)

    # TODO assign samples to temporal bins
    # bin_samples()

    print("INFO: Finished updating {} samples".format(len(samples)))
