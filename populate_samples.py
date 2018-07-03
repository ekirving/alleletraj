#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function

import os
import httplib2
import google_sheets as gs

from pipeline_utils import *


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


def sync_radiocarbon_dates():
    """
    Fetch all the radiocarbon dates
    """

    dbc = db_conn()

    print("INFO: Synchronising radiocarbon dates")

    try:
        # get the google sheet details
        sheet = RADIOCARBON_SHEET[SPECIES]

    except KeyError:

        print("WARNING: No radiocarbon spreadsheet for {}".format(SPECIES))
        return

    # fetch all the manual age mappings
    records = fetch_google_sheet(sheet['id'], sheet['tabs'], sheet['cols'])

    # update the `sample_dates` table
    for record in records:
        if record['accession'] is not None:
            record['confident'] = 'Yes'
            dbc.save_record('sample_dates_c14', record)


def confirm_age_mapping():
    """
    Make sure that all the free-text dates have proper numeric mappings
    """

    dbc = db_conn()

    print("INFO: Confirming age mappings")

    try:
        # get the google sheet details
        sheet = AGE_MAP[SPECIES]

    except KeyError:

        print("WARNING: No age mapping for {}".format(SPECIES))
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
            print(age['age'])

        quit()


def mark_valid_samples():
    """
    Samples are valid if they are from Europe and have a BAM file.
    """

    if SPECIES != 'pig':
        # TODO make this work for all species not just pigs
        raise Exception('Not implemented yet for {}'.format(SPECIES))

    dbc = db_conn()

    print("INFO: Marking valid samples")

    # make a list of permissible countries
    europe = "','".join(EUROPE)

    dbc.execute_sql("""
        UPDATE samples
          JOIN sample_files
            ON sample_files.sample_id = samples.id
           SET valid = 1
         WHERE country IN ('{europe}')""".format(europe=europe))


def bin_samples():
    """
    Assign samples to temporal bins
    """

    dbc = db_conn()

    # sample_bins
    dbc.execute_sql("TRUNCATE TABLE sample_bins")

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


def populate_pig_samples():

    # open a db connection
    dbc = db_conn()

    # get the google sheet details
    sheet = GOOGLE_SHEET[SPECIES]

    # fetch all the samples from the GoogleDoc spreadsheet
    samples = fetch_google_sheet(sheet['id'], sheet['tabs'], sheet['cols'])

    print("INFO: Updating {} samples".format(len(samples)))

    bam_files = defaultdict(list)

    # load the BAM file paths
    with open('./data/bam_files_{}.txt'.format(SPECIES), 'r') as fin:
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
    sync_radiocarbon_dates()

    # make sure that all the free-text dates have proper numeric mappings
    confirm_age_mapping()

    # samples are valid if they are from Europe and have a BAM file
    mark_valid_samples()

    # assign samples to temporal bins
    bin_samples()

    print("INFO: Finished updating {} samples".format(len(samples)))


def populate_horse_samples():

    # open a db connection
    dbc = db_conn()

    # get the google sheet details
    sheet = GOOGLE_SHEET[SPECIES]

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
