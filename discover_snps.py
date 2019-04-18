#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function

from time import time
from datetime import timedelta

from populate_coverage import *
from pipeline_qtls import *

from pipeline_utils import *


def reset_flags(chrom):
    """
    Reset all the analysis flags to NULL
    """
    dbc = db_conn()

    dbc.execute_sql("""
         UPDATE sample_reads
            SET quality = NULL,
                called = NULL
          WHERE chrom = '{chrom}'
          """.format(chrom=chrom))


def apply_quality_filters(chrom):
    """
    Apply MIN_BASE_QUAL, MIN_MAPPING_QUAL and SOFT_CLIP_DIST quality filters.
    """
    dbc = db_conn()

    dbc.execute_sql("""
        UPDATE sample_reads sr
           SET sr.quality = 1
         WHERE sr.chrom = '{chrom}'
           AND sr.baseq >= {baseq}
           AND sr.mapq >= {mapq}
           AND sr.dist > {clip}
           """.format(chrom=chrom, baseq=MIN_BASE_QUAL, mapq=MIN_MAP_QUAL, clip=SOFT_CLIP_DIST))


def choose_random_read(chrom, population):
    """
    Choose a random read from those that pass quality filters
    """
    dbc = db_conn()

    dbc.execute_sql("""
        UPDATE sample_reads
          JOIN (  
                  SELECT SUBSTRING_INDEX(GROUP_CONCAT(sr.id ORDER BY RAND()), ',',  1) id
                    FROM sample_reads sr
                    JOIN modern_snps ms
                      ON ms.population = '{population}'
                     AND ms.chrom = sr.chrom
                     AND ms.site = sr.site
                   WHERE sr.chrom = '{chrom}'
                     AND sr.quality = 1
                     AND sr.base IN (ms.ancestral, ms.derived)
                GROUP BY sr.chrom, sr.site, sr.sample_id
                
               ) AS rand ON rand.id = sample_reads.id
           SET called = 1
           """.format(population=population, chrom=chrom))


def apply_genotype_filters(chrom):
    """
    Apply MIN_GENO_QUAL quality filter.
    """
    dbc = db_conn()

    dbc.execute_sql("""
        UPDATE sample_reads sr
           SET sr.quality = 1,
               sr.called = 1 
         WHERE sr.chrom = '{chrom}'
           AND sr.genoq >= {genoq}
           """.format(chrom=chrom, genoq=MIN_GENO_QUAL))


def discover_snps(population):
    """
    Run queries to mark callable SNPs in the ancient populations.

    Uses mysql table partitioning (w/ MyISAM engine) at the chromosome level to prevent disk swapping caused by massive
    RAM usage for tables with billions of records.
    """

    start = began = time()

    # chunk all the queries by chrom (otherwise we get massive temp tables as the results can't be held in memory)
    chroms = CHROM_SIZE[SPECIES].keys()

    print("INFO: Starting SNP discovery")

    print("INFO: Resetting analysis flags... ", end='')

    for chrom in chroms:
        reset_flags(chrom)

    print("({}).".format(timedelta(seconds=time() - start)))
    start = time()

    print("INFO: Applying quality filters... ", end='')

    for chrom in chroms:
        apply_quality_filters(chrom)

    print("({}).".format(timedelta(seconds=time() - start)))
    start = time()

    print("INFO: Choosing a random read from those that pass quality filters... ", end='')

    for chrom in chroms:
        choose_random_read(chrom, population)

    print("({}).".format(timedelta(seconds=time() - start)))
    start = time()

    print("INFO: Applying genotype quality filters to diploid calls... ", end='')

    for chrom in chroms:
        apply_genotype_filters(chrom)

    print("({}).".format(timedelta(seconds=time() - start)))

    print("SUCCESS: Finished the SNP discovery ({})".format(timedelta(seconds=time() - began)))
