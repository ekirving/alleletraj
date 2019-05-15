#!/usr/bin/env python
# -*- coding: utf-8 -*-

# standard modules
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


class LoadSamples(utils.PipelineTask):
    """
    Load all the samples into the db.

    :type species: str
    :type ancient: bool
    """
    species = luigi.Parameter()
    ancient = luigi.BoolParameter(default=False)

    def requires(self):
        yield CreateDatabase(self.species)
        yield ExternalCSV(self.species, self.ancient)

    def output(self):
        return luigi.LocalTarget('data/db/{}-{}.log'.format(self.basename, self.classname))

    def run(self):
        _, csv_file = self.input()

        dbc = self.db_conn()

        with csv_file.open('r') as fin:
            data = csv.DictReader(fin)

            # strip bad Excel chars and lowercase names
            data.unicode_fieldnames = [re.sub(r'\W+', '', field).lower() for field in data.unicode_fieldnames]

            for row in data:
                sample = {
                    'name':       row['sample'],
                    'population': row['population'],
                    'ancient':    1 if self.ancient else 0,
                    'breed':      row.get('breed'),
                    'site':       row.get('site'),
                    'age_int':    row.get('bp'),
                    'age':        row.get('age'),
                    'period':     row.get('period'),
                    'lat':        row.get('lat'),
                    'long':       row.get('lat'),
                    'sex':        row['sex'][0].upper() if row.get('sex') else None,  # first letter capital
                    'path':       row.get('path')
                }

                sample_id = dbc.save_record('samples', sample)

                # are the accessions paired end or not
                paired = 1 if row.get('librarylayout') == 'PAIRED' else 0

                # get the SRA run accessions
                accessions = [acc.strip() for acc in row['accessions'].split(';') if acc.strip() != '']

                for accession in accessions:
                    dbc.save_record('sample_runs', {'sample_id': sample_id, 'accession': accession, 'paired': paired})

        with self.output().open('w') as fout:
            fout.write('Done!')


class CreateSampleBins(utils.PipelineTask):
    """
    Create the sample bins.

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):
        return LoadSamples(self.species, ancient=True)

    def output(self):
        return luigi.LocalTarget('data/db/{}-{}.log'.format(self.basename, self.classname))

    def run(self):
        dbc = self.db_conn()

        dbc.execute_sql("TRUNCATE TABLE sample_bins")

        # get the maximum date
        max_age = dbc.get_records_sql("""
            SELECT MAX(s.age_int) max_age
              FROM samples s
               """, key=None)[0].pop('max_age')

        # round to nearest multiple of BIN_WIDTH
        bin_max = int(math.ceil(float(max_age) / BIN_WIDTH) * BIN_WIDTH)

        # iterate over each temporal bin
        for bin_upper in range(0, bin_max, BIN_WIDTH):
            bin_lower = bin_upper + BIN_WIDTH

            sample_bin = {
                'name': '{} - {} BP'.format(bin_lower, bin_upper),
                'lower': bin_lower,
                'upper': bin_upper + 1
            }

            dbc.save_record('sample_bins', sample_bin)

        with self.output().open('w') as fout:
            fout.write('Done!')


class BinSamples(utils.PipelineTask):
    """
    Assign samples to temporal bins.

    :type species: str
    """
    species = luigi.Parameter()

    db_lock_tables = ['samples']

    def requires(self):
        return CreateSampleBins(self.species)

    def output(self):
        return luigi.LocalTarget('data/db/{}-{}.log'.format(self.basename, self.classname))

    # noinspection SqlWithoutWhere
    def run(self):
        dbc = self.db_conn()

        # reset sample bins (just in case)
        dbc.execute_sql("""
            UPDATE samples s
               SET s.bin_id = NULL""")

        # assign sample bins
        dbc.execute_sql("""
            UPDATE samples s
              JOIN sample_bins sb
                ON s.age_int BETWEEN sb.upper AND sb.lower
               SET s.bin_id = sb.id""")

        # reset sample counts
        dbc.execute_sql("""
            UPDATE sample_bins sb
              SET sb.num_samples = NULL""")

        # update the the sample count for the bins
        dbc.execute_sql("""
            UPDATE sample_bins sb
              JOIN (  SELECT sb.id, COUNT(*) AS cnt
                        FROM samples s
                        JOIN sample_bins sb
                          ON sb.id = s.bin_id
                    GROUP BY sb.id
                   ) AS num
              ON num.id = sb.id
            SET sb.num_samples = num.cnt""")

        with self.output().open('w') as fout:
            fout.write('Done!')


class LoadSamplesPipeline(utils.PipelineWrapperTask):
    """
    Load both modern and ancient samples, and group ancient samples into temporal bins

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
