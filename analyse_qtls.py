#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function

from db_conn import db_conn
from time import time
from datetime import timedelta

from populate_qtls import *
from populate_coverage import *

# the number of SNPs to model per QTL
SNPS_PER_QTL = 3


def calculate_summary_stats(species, chrom):
    """
    Calculate some summary stats
    """
    dbc = db_conn()

    # remove any existing stats for this chromosome
    dbc.delete_records('qtl_stats', {'species': species, 'chrom': chrom})

    dbc.execute_sql("""
     INSERT INTO qtl_stats (species, qtl_id, chrom, class, type, name, Pvalue, significance, snps, max_samples, 
                            avg_samples, max_reads, avg_reads)
          SELECT species, qtl_id, chrom, class, type, name, Pvalue, significance,
                 COUNT(qtl_id) AS snps,
                 MAX(num_samples) max_samples,
                 AVG(num_samples) avg_samples,
                 MAX(num_reads) max_reads,
                 AVG(num_reads) avg_reads
            FROM (
                      SELECT q.species, q.id AS qtl_id, q.pubmed_id, q.Pvalue, q.significance,
                             t.class, t.type, t.name,
                             sr.chrom, sr.site,
                             COUNT(DISTINCT sr.sample_id) num_samples,
                             COUNT(sr.id) num_reads
                        FROM qtls q
                        JOIN traits t
                          ON t.id = q.trait_id
                        JOIN sample_reads sr
                          ON sr.chrom = q.chrom
                         AND sr.snp = 1
                         AND sr.site BETWEEN q.start and q.end
                       WHERE q.species = '{species}'
                         AND q.valid = 1
                         AND sr.chrom = '{chrom}'
                    GROUP BY q.id, sr.chrom, sr.site

                  ) as snps
        GROUP BY snps.qtl_id""".format(species=species, chrom=chrom))


def analyse_qtl_snps(species, chrom):
    """
    Find the best SNPs for each QTL
    """
    dbc = db_conn()

    # count the number of reads for each SNP
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
                     AND sr.snp = 1
                    JOIN samples s
                      ON s.id = sr.sample_id
                     AND s.species = q.species
                   WHERE q.species = '{species}'
                     AND q.chrom = '{chrom}'
                     AND q.valid = 1
                GROUP BY qs.id

                ) AS num
                  ON num.id = qtl_snps.id

           SET qtl_snps.num_reads = num.num_reads""".format(species=species, chrom=chrom))

    # # choose the best SNPs for each QTL (based on number of reads and closeness to the GWAS peak)
    # dbc.execute_sql("""
    #     UPDATE qtl_snps
    #       JOIN (
    #               SELECT qtl_id,
    #                      SUBSTRING_INDEX(
    #                          GROUP_CONCAT(qs.id ORDER BY num_reads DESC, ABS(v.start-ms.site) ASC),
    #                          ',', {num_snps}) qtl_snps
    #                 FROM qtls q
    #                 JOIN ensembl_variants v
    #                   ON v.rsnumber = q.peak
    #                 JOIN qtl_snps qs
    #                   ON qs.qtl_id = q.id
    #                 JOIN modern_snps ms
    #                   ON ms.id = qs.modsnp_id
    #                WHERE q.species = '{species}'
    #                  AND q.chrom = '{chrom}'
    #                  AND q.valid = 1
    #                  AND qs.num_reads IS NOT NULL
    #             GROUP BY qtl_id
    #
    #             ) AS best
    #               ON qtl_snps.qtl_id = best.qtl_id
    #              AND FIND_IN_SET(qtl_snps.id, best.qtl_snps)
    #
    #         SET qtl_snps.best = 1""".format(species=species, chrom=chrom, num_snps=SNPS_PER_QTL))


def analyse_qtls(species):
    """
    Run queries to pick the best QTLs and SNPs to run the selection scan on.
    """

    start = began = time()

    # chunk all the queries by chrom (otherwise we get massive temp tables as the results can't be held in memory)
    chroms = CHROM_SIZE[species].keys()

    print("INFO: Starting QTL analysis for {}".format(species))

    # # -----------------------------------------------
    # print("INFO: Calculating some summary stats... ", end='')
    # # -----------------------------------------------
    #
    # # TODO refactor to use qtl_snps
    # for chrom in chroms:
    #     calculate_summary_stats(species, chrom)
    #
    # print("({}).".format(timedelta(seconds=time() - start)))
    # start = time()

    print("INFO: Analysing {} QTL SNPs... ".format(species), end='')

    for chrom in chroms:
        analyse_qtl_snps(species, chrom)

    print("({}).".format(timedelta(seconds=time() - start)))

    print("SUCCESS: Finished the {} SNP discovery ({})".format(species, timedelta(seconds=time() - began)))
