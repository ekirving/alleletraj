#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function

import unicodecsv as csv

from pipeline_utils import *

def generate_mc1r_snp_input(population):
    """
    Custom code for handling the PCR data for modsnp #71891
    """

    dbc = Database()

    print("INFO: Generating sample input file for SNP #71891")

    samples = dbc.get_records_sql("""
        # get the ancient frequencies in each bin
        SELECT SUM(CASE s.mc1r_snp 
                       WHEN 'A/A' THEN 2 
                       WHEN 'G/A' THEN 1 
                       ELSE 0 
                   END) AS derived_count,
               count(s.id)*2 AS sample_size,
               -CAST(SUBSTRING_INDEX(sb.bin, ' - ',  1) AS SIGNED INTEGER) AS bin_high,
               -CAST(SUBSTRING_INDEX(sb.bin, ' - ', -1) AS SIGNED INTEGER) - 1 AS bin_low
          FROM samples s
          JOIN sample_bins sb
            ON sb.sample_id = s.id
          WHERE s.valid = 1 
            AND s.status IN ('Domestic')
            AND s.mc1r_snp IS NOT NULL 
       GROUP BY sb.bin

          UNION

          # add the modern frequency 
          SELECT derived_count, ancestral_count + derived_count, 0, 0
            FROM modern_snps ms
           WHERE ms.id = 71891

        ORDER BY bin_high""", key=None)

    # write the sample input file
    with open("selection/{}-{}-modsnp_71891.input".format(SPECIES, population), "wb") as tsv_file:

        fields = ['derived_count', 'sample_size', 'bin_high', 'bin_low']
        writer = csv.DictWriter(tsv_file, fieldnames=fields, delimiter='\t')

        # write the data to disk
        for sample in samples:
            writer.writerow(sample)