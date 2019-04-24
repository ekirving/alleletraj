#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function

from pipeline_sample_reads import *
from pipeline_utils import *


MIN_BASE_QUAL = 30
MIN_MAP_QUAL = 30

# number of bases to soft clip
SOFT_CLIP_DIST = 3


class ResetFlags(PipelineTask):
    """
    Reset all the analysis flags to NULL

    :type species: str
    :type population: str
    :type chrom: str
    """
    species = luigi.Parameter()
    population = luigi.Parameter()
    chrom = luigi.Parameter()

    db_lock_tables = ['sample_reads_{chrom}']

    def requires(self):
        return SampleReadsPipeline(self.species)

    def output(self):
        return luigi.LocalTarget('db/{}-{}.log'.format(self.basename, self.classname))

    def run(self):
        dbc = self.db_conn()

        exec_time = dbc.execute_sql("""
             UPDATE sample_reads sr
               JOIN samples s
                 ON s.id = sr.sample_id
                SET sr.quality = NULL,
                    sr.called = NULL
              WHERE s.population = '{pop}'
                AND sr.chrom = '{chrom}'
              """.format(pop=self.population, chrom=self.chrom))

        with self.output().open('w') as fout:
            fout.write('Execution took {}'.format(exec_time))


class ApplyQualityFilters(PipelineTask):
    """
    Apply MIN_BASE_QUAL, MIN_MAPPING_QUAL and SOFT_CLIP_DIST quality filters.

    :type species: str
    :type population: str
    :type chrom: str
    """
    species = luigi.Parameter()
    population = luigi.Parameter()
    chrom = luigi.Parameter()

    db_lock_tables = ['sample_reads_{chrom}']

    def requires(self):
        return ResetFlags(self.species, self.population, self.chrom)

    def output(self):
        return luigi.LocalTarget('db/{}-{}.log'.format(self.basename, self.classname))

    def run(self):
        dbc = self.db_conn()

        exec_time = dbc.execute_sql("""
            UPDATE sample_reads sr
              JOIN samples s
                ON s.id = sr.sample_id
               SET sr.quality = 1
             WHERE s.population = '{pop}'
               AND sr.chrom = '{chrom}'
               AND sr.baseq >= {baseq}
               AND sr.mapq >= {mapq}
               AND sr.dist > {clip}
               """.format(pop=self.population, chrom=self.chrom, baseq=MIN_BASE_QUAL, mapq=MIN_MAP_QUAL,
                          clip=SOFT_CLIP_DIST))

        with self.output().open('w') as fout:
            fout.write('Execution took {}'.format(exec_time))


class ChooseRandomRead(PipelineTask):
    """
    Choose a random read from those that pass quality filters

    :type species: str
    :type population: str
    :type chrom: str
    """
    species = luigi.Parameter()
    population = luigi.Parameter()
    chrom = luigi.Parameter()

    db_lock_tables = ['sample_reads_{chrom}']

    def requires(self):
        return ApplyQualityFilters(self.species, self.population, self.chrom)

    def output(self):
        return luigi.LocalTarget('db/{}-{}.log'.format(self.basename, self.classname))

    def run(self):
        dbc = self.db_conn()

        # noinspection SqlResolve
        exec_time = dbc.execute_sql("""
        UPDATE sample_reads sr
          JOIN (  
                  SELECT SUBSTRING_INDEX(GROUP_CONCAT(sr.id ORDER BY RAND()), ',',  1) id
                    FROM sample_reads sr
                    JOIN samples s
                      ON s.id = sr.sample_id
                    JOIN modern_snps ms
                      ON ms.population = s.population
                     AND ms.chrom = sr.chrom
                     AND ms.site = sr.site
                   WHERE s.population = '{pop}'
                     AND sr.chrom = '{chrom}'
                     AND sr.quality = 1
                     AND sr.base IN (ms.ancestral, ms.derived)
                GROUP BY sr.chrom, sr.site, sr.sample_id
                
               ) AS rand ON rand.id = sample_reads.id
         WHERE sr.chrom = '{chrom}'
           SET sr.called = 1
           """.format(pop=self.population, chrom=self.chrom))

        with self.output().open('w') as fout:
            fout.write('Execution took {}'.format(exec_time))


class ApplyGenotypeFilters(PipelineTask):
    """
    Apply MIN_GENO_QUAL quality filter.

    :type species: str
    :type population: str
    :type chrom: str
    """
    species = luigi.Parameter()
    population = luigi.Parameter()
    chrom = luigi.Parameter()

    db_lock_tables = ['sample_reads_{chrom}']

    def requires(self):
        return ChooseRandomRead(self.species, self.population, self.chrom)

    def output(self):
        return luigi.LocalTarget('db/{}-{}.log'.format(self.basename, self.classname))

    def run(self):
        dbc = self.db_conn()

        exec_time = dbc.execute_sql("""
        UPDATE sample_reads sr
          JOIN samples s
            ON s.id = sr.sample_id
           SET sr.quality = 1,
               sr.called = 1 
         WHERE s.population = '{pop}'
           AND sr.chrom = '{chrom}'
           AND sr.genoq >= {genoq}
           """.format(pop=self.population, chrom=self.chrom, genoq=MIN_GENO_QUAL))

        with self.output().open('w') as fout:
            fout.write('Execution took {}'.format(exec_time))


class DiscoverSNPsPipeline(PipelineWrapperTask):
    """
    Run queries to mark callable SNPs in the ancient populations.

    Uses mysql table partitioning (w/ MyISAM engine) at the chromosome level to prevent disk swapping caused by massive
    RAM usage for tables with billions of records.

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):
        for pop in self.populations:
            for chrom in self.chromosomes:
                yield ApplyGenotypeFilters(self.species, pop, chrom)
