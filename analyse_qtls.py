#!/usr/bin/env python
# -*- coding: utf-8 -*-

from db_conn import db_conn
from time import time
from datetime import timedelta
from natsort import natsorted

from populate_qtls import *
from populate_coverage import *

# TODO the minimum number of reads per SNP
MIN_READS_PER_SNP = 12

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


def populate_qtl_snps(species, chrom):
    """
    Now we have ascertained all the modern SNPs, let's find those that intersect with the QTLs.
    """

    dbc = db_conn()

    # insert linking records to make future queries much quicker
    dbc.execute_sql("""
        INSERT INTO qtls_snps (qtl_id, modsnp_id)
             SELECT q.id, ms.id
               FROM qtls q
               JOIN modern_snps ms
                 ON ms.species = q.species
                AND ms.chrom = q.chrom
                AND ms.site BETWEEN q.start AND q.end
              WHERE q.species = '{species}'
                AND q.chrom = '{chrom}'
                AND q.valid = 1
                AND ms.maf >= {maf}""".format(species=species, chrom=chrom, maf=MIN_MAF))


def analyse_qtl_snps(species, chrom):
    """
    Find the best SNPs for each QTL
    """
    dbc = db_conn()

    # TODO remove any existing records for this chromosome

    # count the number of reads for each SNP
    dbc.execute_sql("""
        UPDATE qtls_snps
          JOIN (
                  SELECT qs.id,
                         COUNT(sr.id) num_reads
                    FROM qtls q
                    JOIN qtls_snps qs
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
                  ON num.id = qtls_snps.id

           SET qtls_snps.num_reads = num.num_reads""".format(species=species, chrom=chrom))

    # choose the best SNPs for each QTL (based on number of reads and closeness to the GWAS peak)
    dbc.execute_sql("""
        UPDATE qtls_snps
          JOIN (
                  SELECT qtl_id,
                         SUBSTRING_INDEX(
                             GROUP_CONCAT(qs.id ORDER BY num_reads DESC, ABS(d.site-ms.site) ASC), 
                             ',', {num_snps}) qtl_snps
                    FROM qtls q
                    JOIN dbsnp d
                      ON d.rsnumber = q.peak
                    JOIN qtls_snps qs
                      ON qs.qtl_id = q.id
                    JOIN modern_snps ms
                      ON ms.id = qs.modsnp_id
                   WHERE q.species = '{species}'
                     AND q.chrom = '{chrom}'
                     AND q.valid = 1
                     AND qs.num_reads IS NOT NULL
                GROUP BY qtl_id

                ) AS best
                  ON qtls_snps.qtl_id = best.qtl_id
                 AND FIND_IN_SET(qtls_snps.id, best.qtl_snps)

            SET qtls_snps.best = 1""".format(species=species, chrom=chrom, num_snps=SNPS_PER_QTL))


def analyse_qtls(species):
    """
    Run queries to pick the best QTLs and SNPs to run the selection scan on.
    """

    start = began = time()

    # chunk all the queries by chrom (otherwise we get massive temp tables as the results can't be held in memory)
    chroms = natsorted(CHROM_SIZE[species].keys())

    print "INFO: Starting QTL analysis for %s" % species

    # # -----------------------------------------------
    # print "INFO: Calculating some summary stats... ",
    # # -----------------------------------------------
    #
    # for chrom in chroms:
    #     calculate_summary_stats(species, chrom)
    #
    # print "(%s)." % timedelta(seconds=time() - start)
    # start = time()
    #
    # # -----------------------------------------------
    # print "INFO: Populating {} QTL SNPs... ".format(species),
    # # -----------------------------------------------
    #
    # for chrom in chroms:
    #     populate_qtl_snps(species, chrom)
    # print "(%s)." % timedelta(seconds=time() - start)
    # start = time()

    # -----------------------------------------------
    print "INFO: Analysing {} QTL SNPs... ".format(species),
    # -----------------------------------------------

    for chrom in chroms:
        analyse_qtl_snps(species, chrom)

    print "(%s)." % timedelta(seconds=time() - start)
    start = time()

    # TODO more analysis...

    print "SUCCESS: Finished the %s SNP discovery (%s)" % (species, timedelta(seconds=time() - began))