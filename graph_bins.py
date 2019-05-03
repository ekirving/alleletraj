#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
from alleletraj.database.conn import Database

import unicodecsv as csv

# the arbitrary +/- age uncertainty for median age dates
MEDIAN_AGE_UNCERT = 100

BIN_START = 12000
BIN_END = 0

BIN_WIDTH = 500
BIN_PERCENT = 0.5  # samples must overlap a bin by >= 50%

# open a db connection
dbc = Database('horse')

with open("tsv/all-bins.tsv", "wb") as tsv_file:
    fields = ['bin', 'accession', 'map_prcnt', 'population', 'age', 'confident', 'lower', 'upper', 'width', 'overlap',
              'perct_overlap']

    writer = csv.DictWriter(tsv_file, fieldnames=fields, delimiter='\t')
    writer.writeheader()

    for bin_lower in range(BIN_START, BIN_END, - BIN_WIDTH):
        bin_upper = bin_lower - BIN_WIDTH + 1

        # get all the samples which overlap this bin by >= BIN_OVERLAP
        samples = dbc.get_records_sql("""
             SELECT *,
                    @width := lower - upper AS width,
                    @overlap := LEAST(lower, {binlower}) - GREATEST(upper, {binupper}) AS overlap,
                    @overlap / @width AS perct_overlap
               FROM (
                      # get the most precise dates for each sample, and the normalised wild/dom status
                      SELECT s.accession,
                             s.map_prcnt,
                             CASE
                                 WHEN COALESCE(s.gmm_status, s.status) LIKE '%wild%'     THEN 'Wild'
                                 WHEN COALESCE(s.gmm_status, s.status) LIKE '%domestic%' THEN 'Domestic'
                                 WHEN SUBSTR(`group`, 3, 1) = 'W' THEN 'Wild'
                                 WHEN SUBSTR(`group`, 3, 1) = 'D' THEN 'Domestic'
                                 WHEN haplogroup = 'Y1' THEN 'Domestic'
                                 ELSE 'NA'
                             END AS population,
                             s.age,
                             COALESCE(c14.confident, sd.confident) confident,
                             COALESCE(c14.lower, sd.lower, sd.median + {uncert}) lower,
                             COALESCE(c14.upper, sd.upper, sd.median - {uncert}) upper
                        FROM samples s
                   LEFT JOIN sample_dates sd
                          ON s.age = sd.age
                   LEFT JOIN sample_dates_c14 c14
                          ON c14.accession = s.accession
                       WHERE s.valid = 1
    
                    ) as age
    
              WHERE lower >= {binupper}
                AND upper <= {binlower}
             HAVING overlap/width >= {binpercent}""".format(uncert=MEDIAN_AGE_UNCERT,
                                                            binpercent=BIN_PERCENT,
                                                            binlower=bin_lower,
                                                            binupper=bin_upper), key=None)

        stub = "{lower:05d}-{upper:05d}".format(lower=bin_lower, upper=bin_upper)
        label = "{lower:,} - {upper:,} BP".format(lower=bin_lower, upper=bin_upper)

        if not samples:
            print("WARNING: No samples for bin {label}".format(label=label))
            continue

        # write the data to disk
        for sample in samples:
            sample['bin'] = stub
            writer.writerow(sample)

        # now generate the plot
        # utils.run_cmd(['Rscript', 'rscript/plot-bin-mapping.R', stub, label])
