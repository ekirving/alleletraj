#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function

from pipeline_utils import *


def fetch_gwas_peaks():
    """
    Fetch all the GWAS peaks from the QTL database.
    """
    dbc = db_conn()

    start = time()

    print("INFO: Fetching all the GWAS peaks from the QTL database... ", end='')

    dbc.execute_sql("""
        INSERT 
          INTO ascertainment (qtl_id, type, rsnumber, chrom, site, ref, alt, chip_name, snp_name)
        SELECT q.id AS qtl_id, 'peak' AS type,
               ev.rsnumber, ev.chrom, ev.start AS site, ev.ref, ev.alt, 
               sc.chip_name, GROUP_CONCAT(sc.snp_name) AS snp_name
          FROM (
                  # get a unique list of GWAS peaks     
                  SELECT min(id) AS id, peak, site
                    FROM qtls
                   WHERE associationType = 'Association'
                     AND valid = 1
                GROUP BY peak 

               ) AS q
          JOIN ensembl_variants ev
            ON ev.rsnumber = q.peak
     LEFT JOIN snpchip sc
            ON sc.rsnumber = ev.rsnumber
      GROUP BY ev.rsnumber""")

    print("({}).".format(timedelta(seconds=time() - start)))


def fetch_gwas_flanking_snps():
    """
    Fetch the best flanking SNPs for each GWAS peak.
    """
    dbc = db_conn()

    start = time()

    print("INFO: Fetching the best flanking SNPs for each GWAS peak... ", end='')

    # get the best flanking SNPs for each QTL
    qtls = dbc.get_records_sql("""
        SELECT q.id, q.chrom, q.site,
                
                SUBSTRING_INDEX(
                    GROUP_CONCAT(
                        IF(ms.site < q.site, ms.variant_id, NULL) 
                          ORDER BY ms.snpchip_id IS NULL, qs.num_reads DESC, RAND()
                    ), ',',  {num_snps}) left_flank,
                
                SUBSTRING_INDEX(
                    GROUP_CONCAT(
                        IF(ms.site > q.site, ms.variant_id, NULL) 
                          ORDER BY ms.snpchip_id IS NULL, qs.num_reads DESC, RAND()
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
      GROUP BY q.id""".format(num_snps=QTL_FLANK_NUM_SNPS))

    print("({}).".format(timedelta(seconds=time() - start)))

    start = time()

    print("INFO: Loading the flanking SNPs into the database... ", end='')

    # we have to do this iteratively, as FIND_IN_SET() performs terribly
    for qtl_id, qtl in qtls.iteritems():

        # merge the flanking SNP modsnp_ids
        modsnps = ','.join(flank for flank in [qtl['left_flank'], qtl['right_flank']] if flank)

        dbc.execute_sql("""
            INSERT 
              INTO ascertainment (qtl_id, type, rsnumber, chrom, site, ref, alt, chip_name, snp_name)
            SELECT {qtl_id},
                   IF(ev.start > {site}, 'right', 'left') AS type,
                   ev.rsnumber, ev.chrom, ev.start AS site, ev.ref, ev.alt,
                   sc.chip_name, GROUP_CONCAT(sc.snp_name) AS snp_name
              FROM ensembl_variants ev
         LEFT JOIN snpchip sc
                ON sc.rsnumber = ev.rsnumber
             WHERE ev.type = 'SNV'
               AND ev.id IN ({modsnps})
          GROUP BY ev.rsnumber
               """.format(qtl_id=qtl['id'], site=qtl['site'], modsnps=modsnps))

    print("({}).".format(timedelta(seconds=time() - start)))


def fetch_selective_sweep_snps():
    """
    Fetch the best SNPs from the selective sweep loci.
    """
    dbc = db_conn()

    start = time()

    print("INFO: Fetching the best SNPs for each selective sweep locus... ", end='')

    # get the best SNPs for each sweep
    qtls = dbc.get_records_sql("""
        SELECT ss.qtl_id AS id, 
               SUBSTRING_INDEX(
                  GROUP_CONCAT(ms.id ORDER BY ss.p, ABS(ss.site-ms.site)), 
                  ',',  {num_snps}) AS snps
          FROM sweep_snps ss
          JOIN modern_snps ms
            ON ss.chrom = ms.chrom
           AND ms.site BETWEEN ss.site - {offset} AND ss.site + {offset}
          JOIN qtl_snps qs
            ON qs.modsnp_id = ms.id
           AND qs.qtl_id = ss.qtl_id 
      GROUP BY ss.qtl_id
           """.format(num_snps=SWEEP_NUM_SNPS, offset=SWEEP_PEAK_WIDTH/2))

    print("({}).".format(timedelta(seconds=time() - start)))

    start = time()

    print("INFO: Loading the sweep SNPs into the database... ", end='')

    # we have to do this iteratively, as FIND_IN_SET() performs terribly
    for qtl_id, qtl in qtls.iteritems():

        dbc.execute_sql("""
            INSERT 
              INTO ascertainment (qtl_id, type, rsnumber, chrom, site, ref, alt, chip_name, snp_name)
            SELECT {qtl_id}, 'sweep' AS type, ev.rsnumber, ms.chrom, ms.site, 
                   COALESCE(ev.ref, ms.ancestral) ref, COALESCE(ev.alt, ms.derived) alt, 
                   sc.chip_name, GROUP_CONCAT(sc.snp_name) AS snp_name
              FROM modern_snps ms
         LEFT JOIN ensembl_variants ev
                ON ev.id = ms.variant_id
         LEFT JOIN snpchip sc
                ON sc.rsnumber = ev.rsnumber
             WHERE ms.id IN ({modsnps})
          GROUP BY ms.id
               """.format(qtl_id=qtl['id'], modsnps=qtl['snps']))


    print("({}).".format(timedelta(seconds=time() - start)))


def fetch_mc1r_snps():
    """
    Get all the dnsnp SNPs which fall within the MC1R gene.

    See https://www.ensembl.org/sus_scrofa/Gene/Summary?g=ENSSSCG00000020924&db=core
    """

    dbc = db_conn()

    start = time()

    print("INFO: Fetching all the dnsnp SNPs within the MC1R gene... ", end='')

    # get the MC1R qtl
    qtl = dbc.get_record('qtls', {'associationType': 'MC1R'})

    dbc.execute_sql("""
        INSERT 
          INTO ascertainment (qtl_id, type, rsnumber, chrom, site, ref, alt, chip_name, snp_name)
        SELECT {qtl_id},
               'MC1R',
               ev.rsnumber, ev.chrom, ev.start AS site, ev.ref, ev.alt,
               sc.chip_name, GROUP_CONCAT(sc.snp_name) AS snp_name
          FROM ensembl_genes eg
          JOIN ensembl_variants ev
            ON ev.chrom = eg.chrom
           AND ev.start BETWEEN eg.start AND eg.end
           AND ev.type = 'SNV'
           AND LENGTH(ev.alt) = 1
     LEFT JOIN snpchip sc
            ON sc.rsnumber = ev.rsnumber
         WHERE eg.gene_id = '{gene_id}'
      GROUP BY ev.rsnumber
           """.format(qtl_id=qtl['id'], gene_id=MC1R_GENE_ID))

    print("({}).".format(timedelta(seconds=time() - start)))


def fetch_neutral_snps():
    """
    Get neutral SNPs (excluding all QTLs and gene regions, w/ buffer)
    """

    dbc = db_conn()

    start = time()

    print("INFO: Fetching neutral SNPs... ", end='')


    dbc.execute_sql("""
        INSERT
          INTO ascertainment (qtl_id, type, rsnumber, chrom, site, ref, alt, chip_name, snp_name)
        SELECT q.id, q.associationType,
               ev.rsnumber, ev.chrom, ev.start AS site, ev.ref, ev.alt,
               sc.chip_name, GROUP_CONCAT(sc.snp_name) AS snp_name
          FROM qtls q
          JOIN qtl_snps qs
            ON q.id = qs.qtl_id
          JOIN modern_snps ms
            ON ms.id = qs.modsnp_id
          JOIN ensembl_variants ev
            ON ev.id = ms.variant_id
     LEFT JOIN snpchip sc
            ON sc.rsnumber = ev.rsnumber
         WHERE q.associationType = 'Neutral'
           AND q.valid = 1
      GROUP BY ms.id
      ORDER BY ms.snpchip_id IS NULL, qs.num_reads DESC, rand()
         LIMIT {num_neutral}
           """.format(num_neutral=NUM_NEUTRAL_SNPS))

    print("({}).".format(timedelta(seconds=time() - start)))


def perform_ascertainment():
    """
    Ascertain the best SNPs for a catpure array.
    """

    print("INFO: Starting ascertainment process")

    # clear any existing SNPs
    dbc = db_conn()
    dbc.execute_sql("TRUNCATE TABLE ascertainment")

    # fetch all the GWAS peaks from the QTL database
    fetch_gwas_peaks()

    # fetch the best flanking SNPs for each GWAS peak
    fetch_gwas_flanking_snps()

    # fetch the best SNPs from the selective sweep loci
    fetch_selective_sweep_snps()

    # get all MC1R snps
    fetch_mc1r_snps()

    # get neutral SNPs (excluding all QTLs and gene regions, w/ buffer)
    fetch_neutral_snps()

    # TODO make sure none of our SNPs are close to indels

    # TODO get the ancestral SNPs