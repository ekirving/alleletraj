#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function

from time import time
from datetime import timedelta
from natsort import natsorted

from populate_coverage import *
from populate_qtls import *

# the minimum phred scaled genotype quality (30 = 99.9%)
MIN_BASE_QUAL = 30
MIN_MAP_QUAL = 30
MIN_GENO_QUAL = 30

# number of bases to soft clip
SOFT_CLIP_DIST = 3
from pipeline_utils import *

def reset_flags(species, chrom):
    """
    Reset all the analysis flags to NULL
    """
    dbc = db_conn()

    dbc.execute_sql("""
         UPDATE sample_reads
           JOIN samples
             ON samples.id = sample_reads.sample_id
            SET quality = NULL,
                called = NULL,
                snp = NULL
          WHERE species = '{species}'
            AND chrom = '{chrom}'
          """.format(species=species, chrom=chrom))


def apply_quality_filters(species, chrom):
    """
    Apply MIN_BASE_QUAL, MIN_MAPPING_QUAL and SOFT_CLIP_DIST quality filters.
    """
    dbc = db_conn()

    dbc.execute_sql("""
        UPDATE sample_reads
          JOIN samples 
            ON samples.id = sample_reads.sample_id
           SET sample_reads.quality = 1
         WHERE samples.species = '{species}'
           AND sample_reads.chrom = '{chrom}'
           AND sample_reads.baseq >= {baseq}
           AND sample_reads.mapq >= {mapq}
           AND sample_reads.dist > {clip}
           """.format(species=species, chrom=chrom, baseq=MIN_BASE_QUAL, mapq=MIN_MAP_QUAL, clip=SOFT_CLIP_DIST))


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
                      ON sr.sample_id = s.id
                    JOIN modern_snps ms
                      ON ms.chrom = sr.chrom
                     AND ms.site = sr.site
                     AND ms.species = s.species
                   WHERE s.species = '{species}'
                     AND sr.chrom = '{chrom}'
                     AND sr.quality = 1
                     AND sr.base IN (ms.ancestral, ms.derived)
                GROUP BY sr.chrom, sr.site, sr.sample_id
                
               ) AS rand ON rand.id = sample_reads.id
           SET called = 1
           """.format(species=species, chrom=chrom))


def apply_genotype_filters(species, chrom):
    """
    Apply MIN_GENO_QUAL quality filter.
    """
    dbc = db_conn()

    dbc.execute_sql("""
        UPDATE sample_reads
          JOIN samples 
            ON samples.id = sample_reads.sample_id
           SET sample_reads.quality = 1,
               sample_reads.called = 1 
         WHERE samples.species = '{species}'
           AND sample_reads.chrom = '{chrom}'
           AND sample_reads.genoq >= {genoq}
           """.format(species=species, chrom=chrom, genoq=MIN_GENO_QUAL))


def call_ancient_snps(species, chrom):
    """
    Mark the sites which contain SNPs
    """
    dbc = db_conn()

    # TODO this may now be redundant
    dbc.execute_sql("""
        UPDATE sample_reads
          JOIN samples
            ON samples.id = sample_reads.sample_id
          JOIN (
                # get the good SNPs
                SELECT a.chrom, a.site
                FROM (
                      # get all the ancient callable biallelic sites
                      SELECT s.species, sr.chrom, sr.site, 
                             GROUP_CONCAT(DISTINCT sr.base) AS alleles
                        FROM samples s
                        JOIN sample_reads sr
                          ON sr.sample_id = s.id
                       WHERE s.species = '{species}'
                         AND sr.chrom = '{chrom}'
                         AND sr.called = 1
                    GROUP BY sr.chrom, sr.site
                      HAVING COUNT(DISTINCT sr.base) = 2

                     ) AS a
    
                     # make sure they have the same 2 alleles as the modern SNPs
                    JOIN modern_snps as ms 
                      ON ms.species = a.species
                     AND ms.chrom = a.chrom
                     AND ms.site = a.site 
                     AND FIND_IN_SET(ms.ancestral, a.alleles)
                     AND FIND_IN_SET(ms.derived, a.alleles)

                ) AS snp ON snp.chrom = sample_reads.chrom 
                        AND snp.site = sample_reads.site
                        
          SET sample_reads.snp = 1
        WHERE samples.species = '{species}'
          AND sample_reads.chrom = '{chrom}' 
          AND sample_reads.called = 1
        """.format(species=species, chrom=chrom))


def discover_snps(species):
    """
    Run queries to mark callable SNPs in the ancient populations.

    Uses mysql table partitioning (w/ MyISAM engine) at the chromosome level to prevent disk swapping caused by massive
    RAM usage for tables with billions of records.

    :param species: The species name
    :rtype: None
    """

    start = began = time()

    # chunk all the queries by chrom (otherwise we get massive temp tables as the results can't be held in memory)
    chroms = natsorted(CHROM_SIZE[species].keys())

    print("INFO: Starting SNP discovery for {}".format(species))

    print("INFO: Resetting analysis flags... ", end='')

    for chrom in chroms:
        reset_flags(species, chrom)

    print("({}).".format(timedelta(seconds=time() - start)))
    start = time()

    print("INFO: Applying quality filters... ", end='')

    for chrom in chroms:
        apply_quality_filters(species, chrom)

    print("({}).".format(timedelta(seconds=time() - start)))
    start = time()

    print("INFO: Choosing a random read from those that pass quality filters... ", end='')

    for chrom in chroms:
        choose_random_read(species, chrom)

    print("({}).".format(timedelta(seconds=time() - start)))
    start = time()

    print("INFO: Applying genotype quality filters to diploid calls... ", end='')

    for chrom in chroms:
        apply_genotype_filters(species, chrom)

    print("({}).".format(timedelta(seconds=time() - start)))
    start = time()

    print("INFO: Marking the sites which contain SNPs... ", end='')

    for chrom in chroms:
        call_ancient_snps(species, chrom)

    print("({}).".format(timedelta(seconds=time() - start)))

    print("SUCCESS: Finished the {} SNP discovery ({})".format(species, timedelta(seconds=time() - began)))
