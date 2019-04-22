#!/usr/bin/env python
# -*- coding: utf-8 -*-

import luigi

from pipeline_discover_snps import DiscoverSNPsPipeline
from pipeline_utils import PipelineTask, PipelineWrapperTask

# the number of SNPs to model per QTL
SNPS_PER_QTL = 3


class CountSNPCoverage(PipelineTask):
    """
    Count the read coverage for each SNP

    :type species: str
    :type population: str
    :type chrom: str
    """
    species = luigi.Parameter()
    population = luigi.Parameter()
    chrom = luigi.Parameter()

    db_lock_tables = ['qtl_snps']

    def requires(self):
        return DiscoverSNPsPipeline(self.species, self.population, self.chrom)

    def output(self):
        return luigi.LocalTarget('db/{}-{}.log'.format(self.basename, self.classname))

    def run(self):
        # open a db connection
        dbc = self.db_conn()

        # TODO add population
        exec_time = dbc.execute_sql("""
            UPDATE qtl_snps
              JOIN (
                      SELECT qs.id,
                             COUNT(sr.id) num_reads
                        FROM qtls q
                        JOIN qtl_snps qs
                          ON qs.qtl_id = q.id
                        JOIN modern_snps ms
                          ON ms.id = qs.modsnp_id
                        JOIN sample_reads sr
                          ON sr.chrom = ms.chrom
                         AND sr.site = ms.site
                         AND sr.called = 1
                        JOIN samples s
                          ON s.id = sr.sample_id
                       WHERE q.chrom = '{chrom}'
                         AND q.valid = 1
                    GROUP BY qs.id
    
                    ) AS num
                      ON num.id = qtl_snps.id
    
               SET qtl_snps.num_reads = num.num_reads""".format(chrom=self.chrom))

        with self.output().open('w') as fout:
            fout.write('Execution took {}'.format(exec_time))


class FindBestSNPs(PipelineTask):
    """
    Choose the best SNPs for each QTL (based on number of reads and closeness to the GWAS peak)

    :type species: str
    :type population: str
    :type chrom: str
    """
    species = luigi.Parameter()
    population = luigi.Parameter()
    chrom = luigi.Parameter()

    db_lock_tables = ['qtl_snps']

    def requires(self):
        return CountSNPCoverage(self.species, self.population, self.chrom)

    def output(self):
        return luigi.LocalTarget('db/{}-{}.log'.format(self.basename, self.classname))

    def run(self):
        # open a db connection
        dbc = self.db_conn()

        # TODO add population
        exec_time = dbc.execute_sql("""
            UPDATE qtl_snps
              JOIN (
                      SELECT qtl_id,
                             SUBSTRING_INDEX(
                                 GROUP_CONCAT(qs.id ORDER BY num_reads DESC, ABS(q.site-ms.site) ASC),
                                 ',', {num_snps}) qtl_snps
                        FROM qtls q
                        JOIN qtl_snps qs
                          ON qs.qtl_id = q.id
                        JOIN modern_snps ms
                          ON ms.id = qs.modsnp_id
                       WHERE q.chrom = '{chrom}'
                         AND q.valid = 1
                         AND qs.num_reads IS NOT NULL
                    GROUP BY qtl_id
    
                    ) AS best
                      ON qtl_snps.qtl_id = best.qtl_id
                     AND FIND_IN_SET(qtl_snps.id, best.qtl_snps)
    
                SET qtl_snps.best = 1""".format(chrom=self.chrom, num_snps=SNPS_PER_QTL))

        with self.output().open('w') as fout:
            fout.write('Execution took {}'.format(exec_time))


class CalculateSummaryStats(PipelineTask):
    """
    Calculate summary stats for each QTL

    :type species: str
    :type population: str
    :type chrom: str
    """
    species = luigi.Parameter()
    population = luigi.Parameter()
    chrom = luigi.Parameter()

    def requires(self):
        return CountSNPCoverage(self.species, self.population, self.chrom)

    def output(self):
        return luigi.LocalTarget('db/{}-{}.log'.format(self.basename, self.classname))

    def run(self):
        # open a db connection
        dbc = self.db_conn()

        # remove any existing stats for this chromosome
        dbc.execute_sql("TRUNCATE qtl_stats")

        # TODO add population
        exec_time = dbc.execute_sql("""
            INSERT INTO qtl_stats (qtl_id, chrom, class, type, name, Pvalue, significance, snps, max_samples, 
                                   avg_samples, max_reads, avg_reads)
                 SELECT qtl_id, chrom, class, type, name, Pvalue, significance,
                        COUNT(qtl_id) AS snps,
                        MAX(num_samples) max_samples,
                        AVG(num_samples) avg_samples,
                        MAX(num_reads) max_reads,
                        AVG(num_reads) avg_reads
                   FROM (
                          SELECT q.id AS qtl_id, q.pubmed_id, q.Pvalue, q.significance,
                                 t.class, t.type, t.name,
                                 sr.chrom, sr.site,
                                 COUNT(DISTINCT sr.sample_id) num_samples,
                                 COUNT(sr.id) num_reads
                            FROM qtls q
                            JOIN traits t
                              ON t.id = q.trait_id
                            JOIN qtl_snps qs
                              ON qs.qtl_id = q.id
                            JOIN modern_snps ms
                              ON ms.id = qs.modsnp_id
                            JOIN sample_reads sr
                              ON sr.chrom = ms.chrom
                             AND sr.site = ms.site
                             AND sr.called = 1
                           WHERE q.valid = 1
                             AND q.chrom = '{chrom}'
                        GROUP BY qs.id

                         ) as snps
               GROUP BY snps.qtl_id""".format(chrom=self.chrom))

        with self.output().open('w') as fout:
            fout.write('Execution took {}'.format(exec_time))


class AnalyseQTLsPipeline(PipelineWrapperTask):
    """
    Run queries to pick the best QTLs and SNPs to run the selection scan on.

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):

        # process SNPs for all populations and all chromosomes
        for pop in self.populations:
            for chrom in self.chromosomes:

                # choose the best SNPs for each QTL
                yield FindBestSNPs(self.species, pop, chrom)

                # calculate summary stats for each QTL
                yield CalculateSummaryStats(self.species, pop, chrom)