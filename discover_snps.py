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
MIN_GENO_QUAL = 30

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
                called = NULL,
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
                  SELECT SUBSTRING_INDEX(GROUP_CONCAT(sr.id ORDER BY RAND()), ',',  1) id
                    FROM samples s
                    JOIN sample_reads sr
                      ON sr.sampleID = s.id
                   WHERE s.species = '%s'
                     AND sr.chrom = '%s'
                     AND sr.quality = 1
                GROUP BY sr.chrom, sr.pos, sr.sampleID
                
               ) AS rand ON rand.id = sample_reads.id
           SET called = 1
           """ % (species, chrom))


def apply_genotype_filters(species, chrom):
    """
    Apply MIN_GENO_QUAL quality filter.
    """
    dbc = db_conn()

    dbc.execute_sql("""
        UPDATE sample_reads
          JOIN samples 
            ON samples.id = sample_reads.sampleID
           SET sample_reads.quality = 1,
               sample_reads.called = 1 
         WHERE samples.species = '%s'
           AND sample_reads.chrom = '%s'
           AND sample_reads.genoq >= %s
           """ % (species, chrom, MIN_GENO_QUAL))


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
                # get the good SNPs
                SELECT a.chrom, a.pos
                FROM (
                      # get all the ancient callable biallelic sites
                      SELECT s.species, sr.chrom, sr.pos, 
                             GROUP_CONCAT(DISTINCT sr.base) AS alleles
                        FROM samples s
                        JOIN sample_reads sr
                          ON sr.sampleID = s.id
                       WHERE s.species = '%s'
                         AND sr.chrom = '%s'
                         AND sr.called = 1
                    GROUP BY sr.chrom, sr.pos
                      HAVING COUNT(DISTINCT sr.base) = 2

                     ) AS a
    
                     # make sure they have the same 2 alleles as the modern SNPs
                    JOIN modern_snps as ms 
                      ON ms.species = a.species
                     AND ms.chrom = a.chrom
                     AND ms.site = a.pos 
                     AND FIND_IN_SET(ms.ancestral, a.alleles)
                     AND FIND_IN_SET(ms.derived, a.alleles)

                ) AS snp ON snp.chrom = sample_reads.chrom 
                        AND snp.pos = sample_reads.pos
                        
          SET sample_reads.snp = 1
        WHERE samples.species = '%s'
          AND sample_reads.chrom = '%s' 
          AND sample_reads.called = 1
        """ % (species, chrom, species, chrom))


def discover_snps(species):
    """
    Run queries to mark callable SNPs in the ancient populations.

    Uses mysql table partitioning (w/ MyISAM engine) at the chromosome level to prevent disk swapping caused by massive
    RAM usage for tables with billions of records.
    """

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

    print "INFO: Applying genotype quality filters to diploid calls... ",

    for chrom in chroms:
        apply_genotype_filters(species, chrom)

    print "(%s)." % timedelta(seconds=time() - start)
    start = time()

    print "INFO: Marking the sites which contain SNPs... ",

    for chrom in chroms:
        call_ancient_snps(species, chrom)

    print "(%s)." % timedelta(seconds=time() - start)

    print "SUCCESS: Finished the %s SNP discovery (%s)" % (species, timedelta(seconds=time() - began))
