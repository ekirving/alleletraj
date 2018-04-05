#!/usr/bin/env python
# -*- coding: utf-8 -*-

from db_conn import db_conn

import unicodecsv as csv

species = 'pig'

# the arbitrary +/- age uncertainty for median age dates
MEDIAN_AGE_UNCERT = 100

BIN_START = 12000
BIN_END = 0

BIN_WIDTH = 500
BIN_PERCENT = 0.5  # samples must overlap a bin by >= 50%

# open a db connection
dbc = db_conn()

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
                             WHEN s.status LIKE '%domestic%' THEN 'Domestic'
                             WHEN s.status LIKE '%wild%'     THEN 'Wild'
                             ELSE NULL
                         END AS status,
                         s.age,
                         COALESCE(c14.confident, sd.confident) confident,
                         COALESCE(c14.lower, sd.lower, sd.median + {uncert}) lower,
                         COALESCE(c14.upper, sd.upper, sd.median - {uncert}) upper
                    FROM samples s
                    JOIN sample_dates sd
                      ON s.age <=> sd.age
               LEFT JOIN sample_dates_c14 c14
                      ON c14.accession = s.accession
                   WHERE s.species = '{species}'
                     AND s.valid = 1

                ) as age

          WHERE lower >= {binupper}
            AND upper <= {binlower}
         HAVING overlap/width >= {binpercent}""".format(species=species,
                                                        uncert=MEDIAN_AGE_UNCERT,
                                                        binpercent=BIN_PERCENT,
                                                        binlower=bin_lower,
                                                        binupper=bin_upper), key=None)

    tsv_path = "tsv/sample_bin-{lower:05d}-{upper:05d}.tsv".format(lower=bin_lower, upper=bin_upper)
    fields = ['accession', 'map_prcnt', 'status', 'age', 'confident', 'lower', 'upper', 'width', 'overlap', 'perct_overlap']

    # write the data to disk
    with open(tsv_path, "wb") as tsv_file:
        writer = csv.DictWriter(tsv_file, fieldnames=fields, delimiter='\t')
        writer.writeheader()

        for sample in samples:
            writer.writerow(sample)

