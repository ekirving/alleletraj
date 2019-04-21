#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function

from time import time
from datetime import timedelta

from pipeline_sample_reads import *
from pipeline_qtls import *

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

    db_lock_tables = ['sample_reads']

    def requires(self):
        return SampleReadsPipeline(self.species)

    def output(self):
        return luigi.LocalTarget('db/{}-{}.log'.format(self.basename, self.classname))

    def run(self):
        dbc = self.db_conn()

        # TODO add population
        exec_time = dbc.execute_sql("""
             UPDATE sample_reads
                SET quality = NULL,
                    called = NULL
              WHERE chrom = '{chrom}'
              """.format(chrom=self.chrom))

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

    db_lock_tables = ['sample_reads']

    def requires(self):
        return ResetFlags(self.species, self.population, self.chrom)

    def output(self):
        return luigi.LocalTarget('db/{}-{}.log'.format(self.basename, self.classname))

    def run(self):
        dbc = self.db_conn()

        # TODO add population
        exec_time = dbc.execute_sql("""
            UPDATE sample_reads sr
               SET sr.quality = 1
             WHERE sr.chrom = '{chrom}'
               AND sr.baseq >= {baseq}
               AND sr.mapq >= {mapq}
               AND sr.dist > {clip}
               """.format(chrom=self.chrom, baseq=MIN_BASE_QUAL, mapq=MIN_MAP_QUAL, clip=SOFT_CLIP_DIST))

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

    db_lock_tables = ['sample_reads']

    def requires(self):
        return ApplyQualityFilters(self.species, self.population, self.chrom)

    def output(self):
        return luigi.LocalTarget('db/{}-{}.log'.format(self.basename, self.classname))

    def run(self):
        dbc = self.db_conn()

        # TODO add population
        exec_time = dbc.execute_sql("""
        UPDATE sample_reads
          JOIN (  
                  SELECT SUBSTRING_INDEX(GROUP_CONCAT(sr.id ORDER BY RAND()), ',',  1) id
                    FROM sample_reads sr
                    JOIN modern_snps ms
                      ON ms.population = '{population}'
                     AND ms.chrom = sr.chrom
                     AND ms.site = sr.site
                   WHERE sr.chrom = '{chrom}'
                     AND sr.quality = 1
                     AND sr.base IN (ms.ancestral, ms.derived)
                GROUP BY sr.chrom, sr.site, sr.sample_id
                
               ) AS rand ON rand.id = sample_reads.id
           SET called = 1
           """.format(population=self.population, chrom=self.chrom))

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

    db_lock_tables = ['sample_reads']

    def requires(self):
        return ChooseRandomRead(self.species, self.population, self.chrom)

    def output(self):
        return luigi.LocalTarget('db/{}-{}.log'.format(self.basename, self.classname))

    def run(self):
        dbc = self.db_conn()

        # TODO add population
        exec_time = dbc.execute_sql("""
        UPDATE sample_reads sr
           SET sr.quality = 1,
               sr.called = 1 
         WHERE sr.chrom = '{chrom}'
           AND sr.genoq >= {genoq}
           """.format(chrom=self.chrom, genoq=MIN_GENO_QUAL))

        with self.output().open('w') as fout:
            fout.write('Execution took {}'.format(exec_time))


class DiscoverSNPs(PipelineWrapperTask):
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
