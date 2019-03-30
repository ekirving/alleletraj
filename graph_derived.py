#!/usr/bin/env python
# -*- coding: utf-8 -*-

import unicodecsv as csv

from pipeline_utils import *


def graph_derived():
    """
    Produce a scatter plot of the oldest observed derived allele vs. the oldest sample at that site.
    """

    # open a db connection
    dbc = db_conn()

    # sql fragment to calculate the median age of each sample
    median = "COALESCE(sd.median, (COALESCE(c14.lower, sd.lower)+COALESCE(c14.upper, sd.upper))/2)"

    # get the age of every covered snp
    reads = dbc.get_records_sql("""
        SELECT sr.id AS read_id,
               {median} AS median
           FROM (
                    # get a unique list of GWAS peaks
                    SELECT min(id) AS id
                      FROM qtls
                     WHERE associationType = 'Association'
                       AND valid = 1
                  GROUP BY peak
                      
                ) AS q
           JOIN qtl_snps qs
             ON q.id = qs.qtl_id
           JOIN modern_snps ms
             ON ms.id = qs.modsnp_id
           JOIN sample_reads sr
             ON sr.chrom = ms.chrom
            AND sr.site = ms.site
            AND sr.called = 1
           JOIN samples s
             ON s.id = sr.sample_id
      LEFT JOIN sample_dates sd
             ON s.age = sd.age
      LEFT JOIN sample_dates_c14 c14
             ON c14.accession = s.accession
        HAVING median IS NOT NULL""".format(median=median), key=None)

    with open("tsv/{}-snps-counts.tsv".format(SPECIES), "wb") as tsv_file:

        fields = ['read_id', 'median']
        writer = csv.DictWriter(tsv_file, fieldnames=fields, delimiter='\t')
        writer.writeheader()

        # write the data to disk
        for read in reads:
            writer.writerow(read)

    # get the dates of the oldest sample and oldest observation of the derived allele for every QTL SNP
    snps = dbc.get_records_sql("""
         SELECT qs.modsnp_id,
                t.class AS trait,
                ms.ancestral,
                ms.derived,
                ms.daf,
                MAX({median}) AS oldest_sample,
                MAX(IF(sr.base = ms.derived, {median}, NULL)) AS oldest_derived
           FROM (
                    # get a unique list of GWAS peaks
                    SELECT min(id) AS id, trait_id
                      FROM qtls
                     WHERE associationType = 'Association'
                       AND valid = 1
                  GROUP BY peak
                      
                ) AS q
           JOIN qtl_snps qs
             ON q.id = qs.qtl_id
           JOIN traits t
             ON t.id = q.trait_id
           JOIN modern_snps ms
             ON ms.id = qs.modsnp_id
           JOIN sample_reads sr
             ON sr.chrom = ms.chrom
            AND sr.site = ms.site
            AND sr.called = 1
           JOIN samples s
             ON s.id = sr.sample_id
      LEFT JOIN sample_dates sd
             ON s.age = sd.age
      LEFT JOIN sample_dates_c14 c14
             ON c14.accession = s.accession
          WHERE qs.num_reads IS NOT NULL
       GROUP BY qs.modsnp_id""".format(median=median), key=None)

    with open("tsv/{}-snps-ages.tsv".format(SPECIES), "wb") as tsv_file:

        fields = ['modsnp_id', 'trait', 'ancestral', 'derived', 'daf', 'oldest_sample', 'oldest_derived']
        writer = csv.DictWriter(tsv_file, fieldnames=fields, delimiter='\t')
        writer.writeheader()

        # write the data to disk
        for snp in snps:
            writer.writerow(snp)

    # now generate the plot
    run_cmd(['Rscript', 'rscript/plot-age-derived.R', '{}-snps-ages'.format(SPECIES), 'Derived Allele vs. Sample Age'])

    # TODO sips is OS X only, replace with "convert"
    # convert to PNG
    run_cmd(["sips -s format png pdf/{species}-snps-ages.pdf --out pdf/{species}-snps-ages.png".format(species=SPECIES)], shell=True)
