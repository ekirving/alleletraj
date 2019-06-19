#!/usr/bin/env python
# -*- coding: utf-8 -*-

# standard modules
import itertools
import math
import re

# third party modules
import luigi
import unicodecsv as csv

# local modules
from alleletraj import utils
from alleletraj.db.load import CreateDatabase

# number of years in a sample bin
BIN_WIDTH = 500


class ExternalCSV(utils.PipelineExternalTask):
    """
    We require a CSV contianing all the sample details.

    N.B. These have been created outside the workflow of this pipeline.

    :type species: str
    :type ancient: bool
    """
    species = luigi.Parameter()
    ancient = luigi.BoolParameter()

    def output(self):
        period = 'ancient' if self.ancient else 'modern'
        return luigi.LocalTarget('data/{}_samples_{}.csv'.format(period, self.species))


class LoadSamples(utils.MySQLTask):
    """
    Load all the samples into the database.

    :type species: str
    :type ancient: bool
    """
    species = luigi.Parameter()
    ancient = luigi.BoolParameter(default=False)

    def requires(self):
        yield CreateDatabase(self.species)
        yield ExternalCSV(self.species, self.ancient)

    def queries(self):
        _, csv_file = self.input()

        with csv_file.open('r') as fin:
            data = csv.DictReader(fin)

            # strip bad Excel chars and lowercase names
            data.unicode_fieldnames = [re.sub(r'\W+', '', field).lower() for field in data.unicode_fieldnames]

            for row in data:
                sample = {
                    'name':       row['sample'],
                    'population': row['population'],
                    'ancient':    1 if self.ancient else 0,
                    'alias':      row.get('alias'),
                    'group_a':    row.get('group_a') if self.ancient else row.get('breed'),
                    'group_b':    row.get('group_b'),
                    'site':       row.get('site'),
                    'location':   row.get('location'),
                    'bp_max':     row.get('bp_max'),
                    'bp_min':     row.get('bp_min'),
                    'bp_median':  row.get('bp_median'),
                    'age':        row.get('age'),
                    'period':     row.get('period'),
                    'lat':        row.get('lat'),
                    'long':       row.get('lat'),
                    'sex':        row.get('sex')[0].upper() if row.get('sex') else None,  # first letter capital
                    'mtdna':      row.get('mtdna'),
                    'sfs':        row.get('sfs'),
                    'path':       row.get('path')
                }

                sample_id = self.dbc.save_record('samples', sample)

                # get the SRA run accessions
                accessions = [acc.strip() for acc in row['accessions'].split(';') if acc.strip() != '']

                # are the accessions paired end or not
                layout = [1 if lib == 'PAIRED' else 0
                          for lib in row.get('librarylayout').split(';') if lib.strip() != '']

                for accession, paired in itertools.izip(accessions, itertools.cycle(layout)):
                    run = {'sample_id': sample_id, 'accession': accession, 'paired': paired}
                    self.dbc.save_record('sample_runs', run)


class CreateSampleBins(utils.MySQLTask):
    """
    Create the sample bins.

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):
        return LoadSamples(self.species, ancient=True)

    def queries(self):
        self.dbc.execute_sql("TRUNCATE TABLE sample_bins")

        # get the maximum date
        max_age = self.dbc.get_records_sql("""
            SELECT MAX(s.bp_median) max_age
              FROM samples s
               """, key=None)[0].pop('max_age')

        # round to nearest multiple of BIN_WIDTH
        max_ceil = int(math.ceil(float(max_age) / BIN_WIDTH) * BIN_WIDTH)

        # TODO improve binning with a clustering algorithm
        for bin_min in range(0, max_ceil, BIN_WIDTH):
            bin_max = bin_min + BIN_WIDTH

            sample_bin = {
                'name': '{} - {} BP'.format(bin_max, bin_min),
                'max': bin_max,
                'min': bin_min + 1
            }

            self.dbc.save_record('sample_bins', sample_bin)


class BinSamples(utils.MySQLTask):
    """
    Assign samples to temporal bins.

    :type species: str
    """
    species = luigi.Parameter()

    db_lock_tables = ['samples']

    def requires(self):
        return CreateSampleBins(self.species)

    # noinspection SqlWithoutWhere
    def queries(self):
        # reset sample bins (just in case)
        self.dbc.execute_sql("""
            UPDATE samples s
               SET s.bin_id = NULL""")

        # assign sample bins
        self.dbc.execute_sql("""
            UPDATE samples s
              JOIN sample_bins sb
                ON s.bp_median BETWEEN sb.min AND sb.max
               SET s.bin_id = sb.id""")

        # reset sample counts
        self.dbc.execute_sql("""
            UPDATE sample_bins sb
              SET sb.num_samples = NULL""")

        # update the the sample count for the bins
        self.dbc.execute_sql("""
            UPDATE sample_bins sb
              JOIN (  SELECT sb.id, COUNT(*) AS cnt
                        FROM samples s
                        JOIN sample_bins sb
                          ON sb.id = s.bin_id
                    GROUP BY sb.id
                   ) AS num
              ON num.id = sb.id
            SET sb.num_samples = num.cnt""")

        # delete unused bins
        self.dbc.execute_sql("""
            DELETE
              FROM sample_bins
             WHERE num_samples IS NULL""")


class LoadAllSamples(utils.PipelineWrapperTask):
    """
    Load both modern and ancient samples, and group ancient samples into temporal bins.

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):

        # load the modern samples
        yield LoadSamples(self.species)

        # load the ancient samples and assign to temporal bins
        yield BinSamples(self.species)


if __name__ == '__main__':
    luigi.run()
