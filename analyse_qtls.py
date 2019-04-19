#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function

from dbconn import DBConn
from time import time
from datetime import timedelta

from pipeline_qtls import *
from pipeline_sample_reads import *

# the number of SNPs to model per QTL
SNPS_PER_QTL = 3


def calculate_summary_stats():
    """
    Calculate summary stats for each QTL
    """
    dbc = DBConn()

    # remove any existing stats for this chromosome
    dbc.execute_sql("TRUNCATE qtl_stats")

    start = time()

    print("INFO: Calculating QTL summary stats... ", end='')

    # chunk the queries by chrom (to avoid temp tables)
    for chrom in self.chromosomes:
        dbc.execute_sql("""
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
               GROUP BY snps.qtl_id""".format(chrom=chrom))

    print("({}).".format(timedelta(seconds=time() - start)))


def count_snp_coverage():
    """
    Count the read coverage for each SNP
    """
    dbc = DBConn()

    start = time()

    print("INFO: Counting the read coverage for each SNP... ", end='')

    # chunk the queries by chrom (to avoid temp tables)
    for chrom in self.chromosomes:
        dbc.execute_sql("""
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
    
               SET qtl_snps.num_reads = num.num_reads""".format(chrom=chrom))

    print("({}).".format(timedelta(seconds=time() - start)))


def find_best_snps():
    """
    Choose the best SNPs for each QTL (based on number of reads and closeness to the GWAS peak)
    """
    dbc = DBConn()

    start = time()

    print("INFO: Finding the {} best SNPs for each QTL... ".format(SNPS_PER_QTL), end='')

    # chunk the queries by chrom (to avoid temp tables)
    for chrom in self.chromosomes:
        dbc.execute_sql("""
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
    
                SET qtl_snps.best = 1""".format(chrom=chrom, num_snps=SNPS_PER_QTL))

    print("({}).".format(timedelta(seconds=time() - start)))


def analyse_qtls():
    """
    Run queries to pick the best QTLs and SNPs to run the selection scan on.
    """

    print("INFO: Starting QTL analysis")

    # count the read coverage for each SNP
    count_snp_coverage()

    # choose the best SNPs for each QTL
    find_best_snps()

    # calculate summary stats for each QTL
    calculate_summary_stats()

    print("SUCCESS: Finished the SNP discovery")
