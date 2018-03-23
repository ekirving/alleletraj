#!/usr/bin/env python
# -*- coding: utf-8 -*-

from db_conn import db_conn
from time import time
from datetime import timedelta

from populate_coverage import *

# the minimum phred scaled genotype quality (30 = 99.9%)
MIN_BASE_QUAL = 30
MIN_MAPPING_QUAL = 30

# number of bases to soft clip
SOFT_CLIP_DIST = 5


def discover_snps(species, tablename, min_baseq, min_mapq, min_dist, max_qtl):

    dbc = db_conn()

    # TODO turn on query profiling - see https://dev.mysql.com/doc/refman/5.7/en/show-profile.html
    # TODO add indexes for all the group by conditions - see https://dev.mysql.com/doc/refman/5.7/en/table-scan-avoidance.html
    # TODO set innodb_buffer_pool_size to 100 GB, then decrease when done - see https://dev.mysql.com/doc/refman/5.7/en/innodb-buffer-pool-resize.html
    # TODO try compressing the tables - see
    # TODO try multi-threading the individual chrom queries (will this improve CPU utilisation?) * requires multiple connections
    # TODO try partitioning tables https://dev.mysql.com/doc/refman/5.7/en/partitioning-overview.html

    start = began = time()

    # chunk all the queries by chrom (otherwise we get massive temp tables as the results can't be held in memory)
    chroms = CHROM_SIZE[species].keys()

    print "INFO: Starting SNP discovery for %s (%s, %s, %s)" % (tablename, min_baseq, min_mapq, min_dist)

    print "INFO: Resetting existing flags... ",

    # TODO add species filter
    for chrom in chroms:
        # clear the derived columns
        dbc.cursor.execute("""
            UPDATE sample_reads
               SET quality = NULL,
                   random = NULL,
                   snp = NULL
             WHERE chrom = '%s'
               AND (quality IS NOT NULL
                 OR random IS NOT NULL
                 OR snp IS NOT NULL)""" % chrom)
        dbc.cnx.commit()

    print "(%s)." % timedelta(seconds=time() - start)
    start = time()

    print "INFO: Applying quality filters... ",

    for chrom in chroms:
        # apply quality filters
        dbc.cursor.execute("""
            UPDATE sample_reads
               SET quality = 1
             WHERE chrom = '%s'
               AND baseq >= %s
               AND mapq >= %s
               AND dist > %s""" % (chrom, min_baseq, min_mapq, min_dist))
        dbc.cnx.commit()

    print "(%s)." % timedelta(seconds=time() - start)
    start = time()

    print "INFO: Choosing a random read from those that pass quality filters... ",

    for chrom in chroms:
        # choose a random read from those that pass quality filters
        dbc.cursor.execute("""
            UPDATE sample_reads
              JOIN (  
                      SELECT substring_index(group_concat(id ORDER BY rand()), ',', 1) id
                        FROM sample_reads
                       WHERE chrom = '%s'
                         AND quality = 1
                    GROUP BY chrom, pos, sampleID
    
                   ) AS rand ON rand.id = sample_reads.id
               SET random = 1""" % chrom)
        dbc.cnx.commit()

    print "(%s)." % timedelta(seconds=time() - start)
    start = time()

    print "INFO: Marking the sites which contain SNPs... ",

    for chrom in chroms:
        # mark the sites which contain SNPs
        dbc.cursor.execute("""
            UPDATE sample_reads
              JOIN (
                      SELECT chrom, pos 
                        FROM sample_reads
                       WHERE chrom = '%s'
                         AND random = 1
                    GROUP BY chrom, pos
                      HAVING COUNT(id) > 1
                         AND COUNT(DISTINCT base) > 1
                       
                    ) AS sub ON sub.chrom = sample_reads.chrom 
                            AND sub.pos = sample_reads.pos
              SET snp = 1
            WHERE sample_reads.random = 1""" % chrom)
        dbc.cnx.commit()

    print "(%s)." % timedelta(seconds=time() - start)
    start = time()

    print "INFO: Making a backup of these results... ",

    dbc.cursor.execute("""
        DROP TABLE IF EXISTS %s""" % tablename)
    dbc.cnx.commit()

    # make a backup of these results
    dbc.cursor.execute("""
        CREATE TABLE %s
              SELECT *
                FROM sample_reads
               WHERE snp = 1""" % tablename)

    print "(%s)." % timedelta(seconds=time() - start)
    start = time()

    tablename2 = tablename + '_sum'

    print "INFO: Calculating some summary stats... ",

    dbc.cursor.execute("""
        DROP TABLE IF EXISTS %s""" % tablename2)
    dbc.cnx.commit()

    # # calculate some summary stats
    # dbc.cursor.execute("""
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
    # dbc.cnx.commit()
    #
    # print "(%s)." % timedelta(seconds=time() - start)

    print "SUCCESS: Finished the SNP discovery (%s)" % timedelta(seconds=time() - began)
