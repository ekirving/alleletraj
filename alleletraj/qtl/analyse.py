#!/usr/bin/env python
# -*- coding: utf-8 -*-

# third party modules
import luigi

# local modules
from alleletraj import utils
from alleletraj.ancient.snps import LoadAncientSNPs

# the number of SNPs to model per QTL
SNPS_PER_QTL = 3


class CountSNPCoverage(utils.MySQLTask):
    """
    Count the read coverage for each SNP.

    :type species: str
    :type chrom: str
    """
    species = luigi.Parameter()
    chrom = luigi.Parameter()

    db_lock_tables = ['qtl_snps']

    def requires(self):
        return LoadAncientSNPs(self.species, self.chrom)

    # noinspection SqlWithoutWhere
    def queries(self):
        # clear any existing values
        self.dbc.execute_sql("""
            UPDATE qtl_snps qs
              JOIN qtls q
                ON q.id = qs.qtl_id 
               SET qs.num_reads = NULL
             WHERE q.chrom = '{chrom}'
               AND qs.num_reads IS NOT NULL""".format(chrom=self.chrom))

        # count the SNPs for each QTL
        self.dbc.execute_sql("""
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
                       WHERE q.chrom = '{chrom}'
                         AND q.valid = 1
                    GROUP BY qs.id
    
                    ) AS num
                      ON num.id = qtl_snps.id
    
               SET qtl_snps.num_reads = num.num_reads""".format(chrom=self.chrom))


class FindBestSNPs(utils.MySQLTask):
    """
    Choose the best SNPs for each QTL (based on number of reads and closeness to the GWAS peak).

    :type species: str
    :type chrom: str
    """
    species = luigi.Parameter()
    chrom = luigi.Parameter()

    db_lock_tables = ['qtl_snps']

    def requires(self):
        return CountSNPCoverage(self.species, self.chrom)

    # noinspection SqlResolve, SqlWithoutWhere
    def queries(self):
        # clear any existing values
        self.dbc.execute_sql("""
            UPDATE qtl_snps qs
              JOIN qtls q
                ON q.id = qs.qtl_id 
               SET qs.best = NULL
             WHERE q.chrom = '{chrom}'
               AND qs.best IS NOT NULL""".format(chrom=self.chrom))

        # find the N best SNPs
        self.dbc.execute_sql("""
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


class CalculateSummaryStats(utils.MySQLTask):
    """
    Calculate summary stats for each QTL.

    :type species: str
    :type chrom: str
    """
    species = luigi.Parameter()
    chrom = luigi.Parameter()

    db_lock_tables = ['qtl_stats']

    def requires(self):
        return FindBestSNPs(self.species, self.chrom)

    def queries(self):
        # remove any existing stats for this chromosome
        self.dbc.execute_sql("""
            DELETE 
              FROM qtl_stats qs
             WHERE qs.chrom = '{chrom}'""".format(chrom=self.chrom))

        self.dbc.execute_sql("""
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
                           WHERE q.chrom = '{chrom}'
                             AND q.valid = 1
                        GROUP BY qs.id

                         ) as snps
               GROUP BY snps.qtl_id""".format(chrom=self.chrom))


class AnalyseQTLsPipeline(utils.PipelineWrapperTask):
    """
    Run queries to pick the best QTLs and SNPs to run the selection scan on.

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):
        for chrom in self.chromosomes:
            # Calculate summary stats for each QTL
            yield CalculateSummaryStats(self.species, chrom)
