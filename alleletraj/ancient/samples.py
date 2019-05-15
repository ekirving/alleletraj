#!/usr/bin/env python
# -*- coding: utf-8 -*-

# standard modules
import math
import re
import csv

# third party modules
import luigi

# local modules
from alleletraj import utils
from alleletraj.db.load import CreateDatabase

# number of years in a sample bin
BIN_WIDTH = 500


class ExternalCSV(utils.PipelineExternalTask):
    """
    We require an ancient samples CSV.

    N.B. These have been created outside the workflow of this pipeline.

    :type species: str
    """
    species = luigi.Parameter()

    def output(self):
        return luigi.LocalTarget('data/ancient_samples_{}.csv'.format(self.species))


class LoadAncientSamples(utils.PipelineTask):
    """
    Load all the ancient samples into the db.

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):
        yield CreateDatabase(self.species)
        yield ExternalCSV(self.species)

    def output(self):
        return luigi.LocalTarget('data/db/{}-{}.log'.format(self.basename, self.classname))

    def run(self):
        _, csv_file = self.input()

        dbc = self.db_conn()

        with open(csv_file, 'r') as fin:
            data = csv.DictReader(fin)

            # strip bad Excel chars and lowercase names
            data._fieldnames = [re.sub(r'\W+', '', field).lower() for field in data.fieldnames]

            for row in data:
                # normalise sex naming (e.g. male -> M)
                sex = row['sex'][0].upper() if row['sex'] else None

                # are the accessions paired end or not
                paired = True if row['librarylayout'] == 'PAIRED' else False

                record = {
                    'sample':     row['sample'],
                    'population': row['population'],
                    'age_int':    row['bp'],
                    'age':        row.get('age', None),
                    'site':       row.get('site', None),
                    'period':     row.get('period', None),
                    'lat':        row.get('lat', None),
                    'long':       row.get('lat', None),
                    'sex':        sex,
                    'paired':     paired,
                    'path':       row.get('path', None)
                }

                ancient_id = dbc.save_record('ancient', record)

                # get the SRA run accessions
                accessions = row.pop('accessions').split(';')

                for accession in accessions:
                    dbc.save_record('ancient_runs', {'ancient_id': ancient_id, 'accession': accession})

        with self.output().open('w') as fout:
            fout.write('Done!')


class CreateSampleBins(utils.PipelineTask):
    """
    Create the sample bins.

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):
        return LoadAncientSamples(self.species)

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

            dbc.save_record('ancient_bins', sample_bin)

        with self.output().open('w') as fout:
            fout.write('Done!')


class BinSamples(utils.PipelineTask):
    """
    Assign samples to temporal bins.

    :type species: str
    """
    species = luigi.Parameter()

    db_lock_tables = ['ancient']

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


class AncientSamplesPipeline(utils.PipelineWrapperTask):
    """
    Load all the ancient samples and group them in temporal bins

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):
        # assign samples to temporal bins
        return BinSamples(self.species)


if __name__ == '__main__':
    luigi.run()
