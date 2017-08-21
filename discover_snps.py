#!/usr/bin/env python
# -*- coding: utf-8 -*-

from db_conn import db_conn

# the minimum phred scaled genotype quality (30 = 99.9%)
MIN_BASE_QUAL = 30
MIN_MAPPING_QUAL = 30

# number of bases to soft clip
SOFT_CLIP_DIST = 5


def discover_snps(tablename, min_baseq, min_mapq, min_dist, max_qtl, norand=False):

    dbc = db_conn()

    # chunk all the queries by chrom (otherwise we get massive temp tables as the results can't be held in memory)
    chroms = dbc.get_records_sql("""
        SELECT DISTINCT chrom
          FROM intervals""")


    print "INFO: Starting SNP discovery for %s (%s, %s, %s)" % (tablename, min_baseq, min_mapq, min_dist)

    for chrom in chroms:
        # clear the derived columns
        dbc.cursor.execute("""
            UPDATE sample_reads
               SET quality = NULL,
                   random = NULL,
                   snp = NULL
             WHERE chrom = %s""" % chrom)
        dbc.cnx.commit()

    print "INFO: Applying quality filters"

    for chrom in chroms:
        # apply quality filters
        dbc.cursor.execute("""
            UPDATE sample_reads
               SET quality = 1
             WHERE chrom = %s
               AND baseq >= %s
               AND mapq >= %s
               AND dist > %s""" % (chrom, min_baseq, min_mapq, min_dist))
        dbc.cnx.commit()

    print "INFO: Choosing a random read from those that pass quality filters"

    for chrom in chroms:
        # choose a random read from those that pass quality filters
        dbc.cursor.execute("""
            UPDATE sample_reads
              JOIN (  
                      SELECT substring_index(group_concat(id ORDER BY rand()), ',', 1) id
                        FROM sample_reads
                       WHERE chrom = %s
                         AND quality = 1
                    GROUP BY sampleID, chrom, pos
    
                   ) AS rand ON rand.id = sample_reads.id
               SET random = 1""" % chrom)
        dbc.cnx.commit()

    print "INFO: Marking the sites which contain SNPs"

    for chrom in chroms:
        # mark the sites which contain SNPs
        dbc.cursor.execute("""
            UPDATE sample_reads
              JOIN (
                      SELECT chrom, pos 
                        FROM sample_reads
                       WHERE chrom = %s
                         AND random = 1
                    GROUP BY chrom, pos
                      HAVING COUNT(id) > 1
                         AND COUNT(DISTINCT base) > 1
                       
                    ) AS sub ON sub.chrom = sample_reads.chrom 
                            AND sub.pos = sample_reads.pos
              SET snp = 1
            WHERE sample_reads.random = 1""" % chrom)
        dbc.cnx.commit()

    print "INFO: Making a backup of these results"

    # make a backup of these results
    dbc.cursor.execute("""
        CREATE TABLE %s
              SELECT *
                FROM sample_reads
               WHERE snp = 1""" % tablename)

    tablename2 = tablename + '_sum'

    print "INFO: Calculating some summary stats"

    # calculate some summary stats
    dbc.cursor.execute("""
        CREATE TABLE %s
              SELECT id, name, Pvalue, significance, length, 
                     COUNT(id) AS snps, 
                     MAX(num_samples) max_samples, 
                     AVG(num_samples) avg_samples
                FROM (
                          SELECT q.id, t.name, 
                                   q.genomeLoc_end-q.genomeLoc_start AS length, 
                                   q.pubmedID, q.Pvalue, q.significance, 
                                   sr.chrom, sr.pos,
                                   COUNT(DISTINCT sr.sampleID) num_samples
                            FROM qtls q
                            JOIN traits t
                              ON t.id = q.traitID
                            JOIN sample_reads sr
                              ON sr.chrom = q.chromosome
                             AND sr.snp = 1
                             AND sr.pos BETWEEN q.genomeLoc_start and q.genomeLoc_end 
                            WHERE (genomeLoc_end - genomeLoc_start) <= %s
                         GROUP BY q.id, sr.chrom, sr.pos
                
                      ) as sub
            GROUP BY sub.id
            ORDER BY max_samples DESC, snps DESC""" % (tablename, max_qtl))
    dbc.cnx.commit()

    print "SUCCESS: Finished the SNP discovery"
