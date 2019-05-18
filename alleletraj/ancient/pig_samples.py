#!/usr/bin/env python
# -*- coding: utf-8 -*-

# standard modules
import os
from collections import defaultdict, OrderedDict

# third party modules
import httplib2
import luigi
from apiclient import discovery

# local modules
from alleletraj import utils, gsheet as gs
from alleletraj.db.load import CreateDatabase

# from pipeline_utils import *

# TODO fix this
ANCIENT_PATH = '/home/ludo/inbox/BAMs/ancient/'

GOOGLE_SHEET = {

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
            ('Final status (MC1R+Morpho+Context)', 'population'),
            ('Genotype MC1R', 'mc1r_snp')
        ])
    },

    # TODO carry over the GPS coordinates
    # HorseSelection_LO4EIP-TRANSFERED
    'horse': {
        'id': '1BMvIwYj-d8t3mpf67rzabrEvDoB8hBZbyS6XfGwcwUU',
        'tabs': ['Ancient'],
        'cols': OrderedDict([
            ('Name', 'accession'),
            ('Status', 'population'),
            ('path', 'path'),
            ('sex', 'sex'),
            ('Age BP', 'age_int'),
            ('Age', 'period'),
            ('Site', 'location'),
        ])
    },

    # Goat Samples for Allele Trajectory
    'gaot': {
        'id': '1HPKfpxKgwz7MDrt3b6U8530zoQYaYte035aOGl9MzgU',
        'tabs': ['Sheet1'],
        'cols': OrderedDict([
            ('Sample ID', 'accession'),
            ('Wild/Dom', 'population'),
            # TODO sex
            ('Approximate Age ', 'age'),
            ('Archaeological Context', 'period'),
            ('Archaeological Site', 'location'),
            ('Country', 'country'),
            # TODO Latitude / Longitude
        ])
    },
}

# list of junk input to mask with NULL
SHEET_NA = ['n/a', 'NA', 'N', '-', '?', 'NULL', 'None', '...', '']

AGE_MAP = {

    # TODO use this for the goats as well
    'pig': {
        'id': '1bH5u_qDaFXJdTyybeahqgF7je17td0FdyOMs_tlECdA',
        'tabs': ['Age Map'],
        'cols': OrderedDict([
            ('Age', 'age'),
            ('Confident', 'confident'),
            ('Lower (BP)', 'lower'),
            ('Upper (BP)', 'upper'),
            ('Median (BP)', 'median'),
        ])
    }
}

RADIOCARBON_SHEET = {

    'pig': {
        'id': '1odoL9hQh87bLLe3yipbo-CKKXLvIgb5n_kfoqSALHi8',
        'tabs': ['All Dates'],
        'cols': OrderedDict([
            ('Extract_No', 'accession'),
            ('From Cal BP (Int Cal13)', 'lower'),
            ('To Cal BP', 'upper'),
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
NON_EUROPE = [
    'Africa', 'Armenia', 'Azerbaijan', 'Cyprus', 'Egypt', 'Egyptian', 'EuroAm', 'Georgia', 'Iran', 'Iraq', 'Israel',
    'Morocco', 'Sudan', 'Syria', 'Tunisia', 'Turkey', 'Turkmenistan', 'United Arab Emirates'
]

# the arbitrary +/- age uncertainty for median age dates
MEDIAN_AGE_UNCERT = 100


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
    service = discovery.build('sheets', 'v4', http=http, discoveryServiceUrl=discovery_url, cache_discovery=False)

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


class PopulatePigSamples(utils.MySQLTask):
    """
    Load all the ancient pig samples into the database.

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):
        return CreateDatabase(self.species)

    def queries(self):
        # get the google sheet details
        sheet = GOOGLE_SHEET[self.species]

        # fetch all the samples from the GoogleDoc spreadsheet
        samples = fetch_google_sheet(sheet['id'], sheet['tabs'], sheet['cols'])

        bam_files = defaultdict(list)

        # load the BAM file paths
        with open('data/bam_files_{}.txt'.format(self.species), 'r') as fin:
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
            sample_id = self.dbc.save_record('samples', sample)

            # save the BAM file paths
            for path in bam_files[accession]:
                bam_file = dict()
                bam_file['ancient_id'] = sample_id
                bam_file['path'] = path

                if not self.dbc.exists_record('ancient_files', bam_file):
                    self.dbc.save_record('ancient_files', bam_file)


class SyncRadiocarbonDates(utils.MySQLTask):
    """
    Fetch all the radiocarbon dates.

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):
        return PopulatePigSamples(self.species)

    def queries(self):
        # get the google sheet details
        sheet = RADIOCARBON_SHEET[self.species]

        # fetch all the manual age mappings
        records = fetch_google_sheet(sheet['id'], sheet['tabs'], sheet['cols'])

        # update the `ancient_dates` table
        for record in records:
            if record['accession'] is not None:
                record['confident'] = 'Yes'
                self.dbc.save_record('ancient_dates_c14', record)


class ConfirmAgeMapping(utils.MySQLTask):
    """
    Make sure that all the free-text dates have proper numeric mappings.

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):
        return PopulatePigSamples(self.species)

    def queries(self):
        # get the google sheet details
        sheet = AGE_MAP[self.species]

        # fetch all the manual age mappings
        records = fetch_google_sheet(sheet['id'], sheet['tabs'], sheet['cols'])

        # update the `ancient_dates` table
        for record in records:
            if record['age'] is not None:
                self.dbc.save_record('ancient_dates', record)

        # check if any samples have an age which is unmapped
        missing = self.dbc.get_records_sql("""
            SELECT DISTINCT s.age
              FROM samples s
         LEFT JOIN sample_dates sd
                ON s.age = sd.age
             WHERE s.age IS NOT NULL
               AND sd.id IS NULL
          GROUP BY s.age""", key=None)

        if missing:
            ages = ', '.join(["'{}'".format(age['age']) for age in missing])
            raise Exception("ERROR: Not all sample ages have numeric mappings - {}".format(ages))


class ConfirmCountryMapping(utils.MySQLTask):
    """
    Make sure that all the free-text countries have been properly mapped to Europe.

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):
        return PopulatePigSamples(self.species)

    def queries(self):
        # make a list of known countries
        world = "','".join(EUROPE + NON_EUROPE)

        missing = self.dbc.get_records_sql("""
            SELECT DISTINCT a.country
              FROM ancient a
             WHERE a.country NOT IN ('{world}')
               """.format(world=world), key=None)

        if missing:
            raise RuntimeError('ERROR: Not all countries have been mapped! {}'.format(missing))


class MarkValidPigs(utils.MySQLTask):
    """
    Pig samples are valid if they are from Europe and have a BAM file or MC1R genotype.

    :type species: str
    """
    species = luigi.Parameter()

    db_lock_tables = ['ancient']

    def requires(self):
        yield SyncRadiocarbonDates(self.species)
        yield ConfirmAgeMapping(self.species)
        yield ConfirmCountryMapping(self.species)

    def queries(self):
        # make a list of permissible countries
        europe = "','".join(EUROPE)

        self.dbc.execute_sql("""
            UPDATE samples s
         LEFT JOIN sample_files sf
                ON sf.sample_id = s.id
               SET s.valid = 1
             WHERE s.country IN ('{europe}')
               AND (sf.id IS NOT NULL 
                OR s.mc1r_snp IS NOT NULL) """.format(europe=europe))


if __name__ == '__main__':
    luigi.run()
