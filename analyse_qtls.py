#!/usr/bin/env python
# -*- coding: utf-8 -*-

from db_conn import db_conn
from time import time
from datetime import timedelta

from populate_qtls import *

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
                      SELECT q.species, q.id AS qtl_id, q.pubmedID, q.Pvalue, q.significance,
                             t.class, t.type, t.name,
                             sr.chrom, sr.pos,
                             COUNT(DISTINCT sr.sampleID) num_samples,
                             COUNT(sr.id) num_reads
                        FROM qtls q
                        JOIN traits t
                          ON t.id = q.traitID
                        JOIN sample_reads sr
                          ON sr.chrom = q.chromosome
                         AND sr.snp = 1
                         AND sr.pos BETWEEN q.start and q.end
                       WHERE q.species = '{species}'
                         AND q.valid = 1
                         AND sr.chrom = '{chrom}'
                    GROUP BY q.id, sr.chrom, sr.pos

                  ) as snps
        GROUP BY snps.qtl_id""".format(species=species, chrom=chrom))


def analyse_qtl_snps(species, chrom):
    """
    Find the best SNPs for each QTL, and link their sample reads to the QTL.
    """
    dbc = db_conn()

    # TODO remove any existing records for this chromosome
    # dbc.delete_records('qtl_reads', {'species': species, 'chrom': chrom})

    dbc.execute_sql("""
     INSERT INTO qtl_reads (qtl_id, read_id)
          SELECT qtl_id, sr.id read_id
            FROM sample_reads sr
            JOIN (
                    # for each QTL, get the best SNPs (i.e. those with the most reads)
                    SELECT snp.qtl_id, snp.chrom,
                           SUBSTRING_INDEX(GROUP_CONCAT(snp.pos ORDER BY snp.num_reads DESC), ',',  {num_snps}) sites
                      FROM (
                              # for each QTL, get the count of reads for each SNP
                              SELECT q.id AS qtl_id,
                                     sr.chrom, sr.pos,
                                     COUNT(sr.id) num_reads
                                FROM qtls q
                                JOIN sample_reads sr
                                  ON sr.chrom = q.chromosome
                                 AND sr.pos BETWEEN q.start and q.end
                               WHERE q.species = '{species}'
                                 AND q.valid = 1
                                 AND sr.chrom = '{chrom}'
                                 AND sr.snp = 1
                            GROUP BY q.id, sr.chrom, sr.pos
                            
                            ) AS snp
                   GROUP BY qtl_id
                   
                 ) AS best
                   ON sr.chrom = best.chrom
                  AND FIND_IN_SET(sr.pos, best.sites)
           WHERE sr.chrom = '{chrom}'
             AND sr.snp = 1""".format(num_snps=SNPS_PER_QTL, species=species, chrom=chrom))


def analyse_qtls(species):
    """
    Run queries to pick the best QTLs and SNPs to run the selection scan on.
    """

    start = began = time()

    # chunk all the queries by chrom (otherwise we get massive temp tables as the results can't be held in memory)
    chroms = CHROM_SIZE[species].keys()

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

    # -----------------------------------------------
    print "INFO: Analysing QTL SNPs... ",
    # -----------------------------------------------

    for chrom in chroms:
        analyse_qtl_snps(species, chrom)

    print "(%s)." % timedelta(seconds=time() - start)
    start = time()





    analyse_qtl_snps(species, chrom)


    print "SUCCESS: Finished the %s SNP discovery (%s)" % (species, timedelta(seconds=time() - began))