#!/usr/bin/env python
# -*- coding: utf-8 -*-

from db_conn import db_conn
from time import time
from datetime import timedelta

from populate_coverage import *
from populate_qtls import *

# the minimum phred scaled genotype quality (30 = 99.9%)
MIN_BASE_QUAL = 30
MIN_MAPPING_QUAL = 30

# number of bases to soft clip
SOFT_CLIP_DIST = 5


def reset_flags(species, chrom):
    """
    Reset all the analysis flags to NULL
    """
    dbc = db_conn()

    dbc.execute_sql("""
         UPDATE sample_reads
           JOIN samples
             ON samples.id = sample_reads.sampleID
            SET quality = NULL,
                random = NULL,
                snp = NULL
          WHERE species = '%s'
            AND chrom = '%s'
          """ % (species, chrom))


def apply_quality_filters(species, chrom):
    """
    Apply MIN_BASE_QUAL, MIN_MAPPING_QUAL and SOFT_CLIP_DIST quality filters.
    """
    dbc = db_conn()

    dbc.execute_sql("""
        UPDATE sample_reads
          JOIN samples 
            ON samples.id = sample_reads.sampleID
           SET sample_reads.quality = 1
         WHERE samples.species = '%s'
           AND sample_reads.chrom = '%s'
           AND sample_reads.baseq >= %s
           AND sample_reads.mapq >= %s
           AND sample_reads.dist > %s
           """ % (species, chrom, MIN_BASE_QUAL, MIN_MAPPING_QUAL, SOFT_CLIP_DIST))


def choose_random_read(species, chrom):
    """
    Choose a random read from those that pass quality filters
    """
    dbc = db_conn()

    dbc.execute_sql("""
        UPDATE sample_reads
          JOIN (  
                  SELECT substring_index(group_concat(sr.id ORDER BY rand()), ',', 1) id
                    FROM samples s
                    JOIN sample_reads sr
                      ON sr.sampleID = s.id
                   WHERE s.species = '%s'
                     AND sr.chrom = '%s'
                     AND sr.quality = 1
                GROUP BY sr.chrom, sr.pos, sr.sampleID

               ) AS rand ON rand.id = sample_reads.id
           SET random = 1
           """ % (species, chrom))


def call_ancient_snps(species, chrom):
    """
    Mark the sites which contain SNPs
    """
    dbc = db_conn()

    dbc.execute_sql("""
        UPDATE sample_reads
          JOIN samples
            ON samples.id = sample_reads.sampleID
          JOIN (
                  SELECT sr.chrom, sr.pos 
                    FROM samples s
                    JOIN sample_reads sr
                      ON sr.sampleID = s.id
                   WHERE s.species = '%s'
                     AND sr.chrom = '%s'
                     AND sr.random = 1
                GROUP BY sr.chrom, sr.pos
                  HAVING COUNT(sr.sampleID) > 1
                     AND COUNT(DISTINCT sr.base) > 1

                ) AS sub ON sub.chrom = sample_reads.chrom 
                        AND sub.pos = sample_reads.pos
          SET sample_reads.snp = 1
        WHERE samples.species = '%s'
          AND sample_reads.chrom = '%s' 
          AND sample_reads.random = 1
        """ % (species, chrom, species, chrom))


def discover_snps(species):

    dbc = db_conn()

    # TODO try multi-threading the individual chrom queries (will this improve CPU utilisation?)
    # TODO try partitioning tables https://dev.mysql.com/doc/refman/5.7/en/partitioning-overview.html

    start = began = time()

    # chunk all the queries by chrom (otherwise we get massive temp tables as the results can't be held in memory)
    chroms = CHROM_SIZE[species].keys()

    print "INFO: Starting SNP discovery for %s" % species

    print "INFO: Resetting analysis flags... ",

    for chrom in chroms:
        reset_flags(species, chrom)

    print "(%s)." % timedelta(seconds=time() - start)
    start = time()

    print "INFO: Applying quality filters... ",

    for chrom in chroms:
        apply_quality_filters(species, chrom)

    print "(%s)." % timedelta(seconds=time() - start)
    start = time()

    print "INFO: Choosing a random read from those that pass quality filters... ",

    for chrom in chroms:
        choose_random_read(species, chrom)

    print "(%s)." % timedelta(seconds=time() - start)
    start = time()

    print "INFO: Marking the sites which contain SNPs... ",

    for chrom in chroms:
        call_ancient_snps(species, chrom)

    print "(%s)." % timedelta(seconds=time() - start)
    start = time()

    # print "INFO: Calculating some summary stats... ",
    #
    # # calculate some summary stats
    # dbc.execute_sql("""
    #     CREATE TABLE %s
    #           SELECT id, name, Pvalue, significance, length,
    #                  COUNT(id) AS snps,
    #                  MAX(num_samples) max_samples,
    #                  AVG(num_samples) avg_samples
    #             FROM (
    #                       SELECT q.id, t.name,
    #                                q.genomeLoc_end-q.genomeLoc_start AS length,
    #                                q.pubmedID, q.Pvalue, q.significance,
    #                                sr.chrom, sr.pos,
    #                                COUNT(DISTINCT sr.sampleID) num_samples
    #                         FROM qtls q
    #                         JOIN traits t
    #                           ON t.id = q.traitID
    #                         JOIN sample_reads sr
    #                           ON sr.chrom = q.chromosome
    #                          AND sr.snp = 1
    #                          AND sr.pos BETWEEN q.genomeLoc_start and q.genomeLoc_end
    #                         WHERE (genomeLoc_end - genomeLoc_start) <= %s
    #                      GROUP BY q.id, sr.chrom, sr.pos
    #
    #                   ) as sub
    #         GROUP BY sub.id
    #         ORDER BY max_samples DESC, snps DESC""" % (tablename2, max_qtl))
    #
    # print "(%s)." % timedelta(seconds=time() - start)

    print "SUCCESS: Finished the %s SNP discovery (%s)" % (species, timedelta(seconds=time() - began))
