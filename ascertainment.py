#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function

from db_conn import db_conn
from time import time
from datetime import timedelta

# the number of flanking SNPs (on either side) to include
NUM_FLANKING_SNPS = 2

def fetch_gwas_peaks(species):
    """
    Fetch all the GWAS peaks from the QTL database.
    """
    dbc = db_conn()

    start = time()

    print("INFO: Fetching all the GWAS peaks from the QTL database... ", end='')

    dbc.execute_sql("""
        INSERT INTO ascertainment (qtl_id, qtl_pos, rsnumber, chrom, site, ref, alt, chip_name, snp_name)
        SELECT q.id AS qtl_id, 'peak' AS qtl_pos,
               ev.rsnumber, ev.chrom, ev.start AS site, ev.ref, ev.alt, 
               ds.chip_name, GROUP_CONCAT(ds.snp_name) AS snp_name
          FROM (
                  # get a unique list of GWAS peaks     
                  SELECT min(id) AS id, peak, site
                    FROM qtls
                   WHERE species = '{species}'
                     AND associationType = 'Association'
                     AND valid = 1
                GROUP BY peak 

               ) AS q
          JOIN ensembl_variants ev
            ON ev.rsnumber = q.peak
     LEFT JOIN dbsnp_snpchip ds
            ON ds.rsnumber = ev.rsnumber
      GROUP BY ev.rsnumber""".format(species=species))

    print("({}).".format(timedelta(seconds=time() - start)))


def fetch_gwas_flanking_snps(species):
    """
    Fetch the best flanking SNPs for each GWAS peak.
    """
    dbc = db_conn()

    start = time()

    print("INFO: Fetching the best flanking SNPs for each GWAS peak.... ", end='')

    # get the best flanking SNPs for each QTL
    qtls = dbc.get_records_sql("""
        SELECT q.id, q.chrom, q.site,
                SUBSTRING_INDEX(
                    GROUP_CONCAT(
                        IF(ms.site < q.site, ms.variant_id, NULL) ORDER BY ms.snpchip_id IS NULL, qs.num_reads DESC, RAND()
                    ), ',',  {num_snps}) left_flank,
                SUBSTRING_INDEX(
                    GROUP_CONCAT(
                        IF(ms.site > q.site, ms.variant_id, NULL) ORDER BY ms.snpchip_id IS NULL, qs.num_reads DESC, RAND()
                    ), ',',  {num_snps}) right_flank
          FROM (
                    # get a unique list of GWAS peaks
                    SELECT min(id) AS id, chrom, peak, site
                      FROM qtls
                     WHERE associationType = 'Association'
                       AND valid = 1
                  GROUP BY peak
                  
               ) AS q
          JOIN qtl_snps qs
            ON qs.qtl_id = q.id
          JOIN modern_snps ms
            ON ms.id = qs.modsnp_id
      GROUP BY q.id""".format(num_snps=NUM_FLANKING_SNPS))

    print("({}).".format(timedelta(seconds=time() - start)))

    start = time()

    print("INFO: Loading the flanking SNPs into the database.... ", end='')

    # we have to do this iteratively, as FIND_IN_SET() performs terribly
    for qtl_id, qtl in qtls.iteritems():

        # merge the flanking SNP modsnp_ids
        modsnps = qtl['left_flank'] + "," + qtl['right_flank']

        dbc.execute_sql("""
            INSERT INTO ascertainment (qtl_id, qtl_pos, rsnumber, chrom, site, ref, alt, chip_name, snp_name)
            SELECT {qtl_id},
                   IF(ev.start > {site}, 'right', 'left') AS qtl_pos,
                   ev.rsnumber, ev.chrom, ev.start AS site, ev.ref, ev.alt,
                   ds.chip_name, GROUP_CONCAT(ds.snp_name) AS snp_name
              FROM ensembl_variants ev
         LEFT JOIN dbsnp_snpchip ds
                ON ds.rsnumber = ev.rsnumber
             WHERE ev.type = 'SNV'
               AND ev.id IN ({modsnps})
          GROUP BY ev.rsnumber
               """.format(qtl_id=qtl['id'], site=qtl['site'], modsnps=modsnps))

    print("({}).".format(timedelta(seconds=time() - start)))


def fetch_selective_sweep_snps(species):
    """
    Fetch the best SNPs from the selective sweep loci.
    """
    dbc = db_conn()

    dbc.execute_sql("""
    """)


def perform_ascertainment(species):
    """
    Ascertain the best SNPs for a catpure array.
    """

    print("INFO: Ascertainment process for {}".format(species))

    # fetch all the GWAS peaks from the QTL database
    # fetch_gwas_peaks(species)

    # fetch the best flanking SNPs for each GWAS peak
    fetch_gwas_flanking_snps(species)

    # fetch the best SNPs from the selective sweep loci
    fetch_selective_sweep_snps(species)