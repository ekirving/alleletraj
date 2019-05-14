#!/usr/bin/env python
# -*- coding: utf-8 -*-

# third party modules
import unicodecsv as csv

# local modules
from alleletraj.db.conn import Database


def generate_mc1r_snp_input(species, population):
    """
    Custom code for handling the PCR data for modsnp #71891
    """

    dbc = Database(species)

    print("INFO: Generating sample input file for SNP #71891")

    samples = dbc.get_records_sql("""
        # get the ancient frequencies in each bin
        SELECT SUM(CASE a.mc1r_snp 
                       WHEN 'A/A' THEN 2 
                       WHEN 'G/A' THEN 1 
                       ELSE 0 
                   END) AS derived_count,
               count(a.id)*2 AS sample_size,
               -CAST(SUBSTRING_INDEX(ab.name, ' - ',  1) AS SIGNED INTEGER) AS bin_high,
               -CAST(SUBSTRING_INDEX(ab.name, ' - ', -1) AS SIGNED INTEGER) - 1 AS bin_low
          FROM ancient a
          JOIN ancient_bins ab
            ON a.ancient_id = ab.id
          WHERE a.valid = 1 
            AND a.population IN ('Domestic')  # TODO fix me
            AND a.mc1r_snp IS NOT NULL 
       GROUP BY ab.bin

          UNION

          # add the modern frequency 
          SELECT derived_count, ancestral_count + derived_count, 0, 0
            FROM modern_snps ms
           WHERE ms.id = 71891

        ORDER BY bin_high""", key=None)

    # write the sample input file
    with open("data/selection/{}-{}-modsnp_71891.input".format(species, population), "wb") as tsv_file:
        fields = ['derived_count', 'sample_size', 'bin_high', 'bin_low']
        writer = csv.DictWriter(tsv_file, fieldnames=fields, delimiter='\t')

        # write the data to disk
        for sample in samples:
            writer.writerow(sample)
