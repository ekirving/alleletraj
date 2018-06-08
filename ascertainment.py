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
        SELECT q.id, 'peak',
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
            SELECT {qtl_id}, IF(ev.start > {site}, 'right', 'left'),
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

    print("INFO: Fetching the {} best SNPs for each selective sweep locus... ".format(SWEEP_NUM_SNPS), end='')

    qtls = dbc.get_records_sql("""
        SELECT qtl_id AS id,
               SUBSTRING_INDEX(
                   GROUP_CONCAT(modsnp_id ORDER BY ss.p, ABS(ss.site-near.site)), 
                   ',',  3) AS snps
          FROM (
                    # find the nearest sweep SNP to each QTL SNP in dbsnp 
                    SELECT ms.id AS modsnp_id,
                           ms.site,
                           SUBSTRING_INDEX(
                               GROUP_CONCAT(ss.id ORDER BY ABS(ss.site-ms.site)), 
                               ',',  1) AS ss_id
                      FROM sweep_snps ss
                      JOIN qtl_snps qs
                        ON qs.qtl_id = ss.qtl_id 
                      JOIN modern_snps ms
                        ON ms.id = qs.modsnp_id
                     WHERE ms.variant_id IS NOT NULL
                  GROUP BY ms.id
                  
               ) AS near
          JOIN sweep_snps ss
            ON ss.id = near.ss_id
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
            SELECT {qtl_id}, 'sweep', 
                   ev.rsnumber, ev.chrom, ev.start AS site, ev.ref, ev.alt,
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
        SELECT {qtl_id}, 'mc1r', 
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
           """.format(qtl_id=qtl['id'], gene_id=MC1R_GENE_ID[SPECIES]))

    print("({}).".format(timedelta(seconds=time() - start)))


def fetch_neutral_snps():
    """
    Get neutral SNPs (excluding all QTLs and gene regions, w/ buffer)
    """

    dbc = db_conn()

    start = time()

    print("INFO: Fetching {:,} neutral SNPs... ".format(NUM_NEUTRAL_SNPS), end='')


    # get the total size of the autosomes
    total = sum(size for chrom, size in CHROM_SIZE[SPECIES].iteritems() if chrom not in ['X', 'Y'])

    # calculate the proportional size of each autosomes
    autosomes = OrderedDict((chrom, float(size) / total)
                         for chrom, size in CHROM_SIZE[SPECIES].iteritems() if chrom not in ['X', 'Y'])

    for chrom, perct in autosomes.iteritems():

        # get the weighted number of SNPs for this chrom
        num_snps = int(round(NUM_NEUTRAL_SNPS/perct))

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
             WHERE q.chrom = '{chrom}'
               AND q.associationType = 'Neutral'
               AND q.valid = 1
          GROUP BY ms.id
          ORDER BY ms.snpchip_id IS NULL, qs.num_reads DESC, rand()
             LIMIT {num_snps}
               """.format(chrom=chrom, num_snps=num_snps))

    print("({}).".format(timedelta(seconds=time() - start)))


def fetch_ancestral_snps():
    """
    Get ancestral SNPs which are in variable in ASD and SUM (Sumatran Sus scrofa)
    """

    dbc = db_conn()

    start = time()

    print("INFO: Fetching {:,} ancestral SNPs... ".format(NUM_ANCESTRAL_SNPS), end='')

    # get the total size of the autosomes
    total = sum(size for chrom, size in CHROM_SIZE[SPECIES].iteritems() if chrom not in ['X', 'Y'])

    # calculate the proportional size of each autosomes
    autosomes = OrderedDict((chrom, float(size) / total)
                         for chrom, size in CHROM_SIZE[SPECIES].iteritems() if chrom not in ['X', 'Y'])

    for chrom, perct in autosomes.iteritems():

        # get the weighted number of SNPs for this chrom
        num_snps = int(round(NUM_ANCESTRAL_SNPS / perct))

        dbc.execute_sql("""
            INSERT
              INTO ascertainment (qtl_id, type, rsnumber, chrom, site, ref, alt, chip_name, snp_name)
            SELECT NULL, 'Ancestral',
                   ev.rsnumber, ev.chrom, ev.start AS site, ev.ref, ev.alt,
                   sc.chip_name, GROUP_CONCAT(sc.snp_name) AS snp_name
              FROM modern_snps asd
              JOIN modern_snps sum
                ON sum.population = 'SUM'
               AND sum.chrom = asd.chrom
               AND sum.variant_id = asd.variant_id
              JOIN ensembl_variants ev
                ON ev.id = asd.variant_id
         LEFT JOIN snpchip sc
                ON sc.rsnumber = ev.rsnumber 
             WHERE asd.population = 'ASD'
               AND asd.chrom = '{chrom}'
               AND asd.derived_count > 1
          GROUP BY asd.id
          ORDER BY asd.snpchip_id IS NULL, rand()
             LIMIT {num_snps}
               """.format(chrom=chrom, num_snps=num_snps))

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

    # get ancestral SNPs which are in variable in ASD and Sumatran scrofa
    fetch_ancestral_snps()

    # TODO make sure none of our SNPs are close to indels