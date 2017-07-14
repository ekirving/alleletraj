#!/usr/bin/env python
# -*- coding: utf-8 -*-

from db_conn import db_conn
from fetch_metadata import fetch_metadata

import pysam as ps

# observations about the QTL db
# - there are many overlapping QTL windows (some very large windows span multiple smaller ones)
# - there are many duplicate windows, ~8k (~44%)


def check_coverage():
    # open a db connection
    dbc = db_conn()

    # get a sorted list of unique QTLs
    qtls = dbc.get_records_sql(
        """SELECT GROUP_CONCAT(id) AS id, chromosome as chr, genomeLoc_start as start, genomeLoc_end as end
             FROM qtls
            WHERE genomeLoc_start IS NOT NULL
              AND genomeLoc_end IS NOT NULL
            GROUP BY chromosome, genomeLoc_start, genomeLoc_end
            ORDER BY chromosome, genomeLoc_start
            LIMIT 10"""  # TODO remove this
    )

    print "INFO: Found %s unique QTLs" % len(qtls)

    # fetch the metadata
    df = fetch_metadata()

    print "INFO: Found %s samples to compute coverage for" % df.shape[0]

    # process each QTL
    for ids, qtl in qtls.iteritems():

        # check all the samples for coverage in this QTL
        for accession, sample in df.iterrows():

            ps.AlignmentFile.fetch()

            # open the BAM file for reading
            with ps.AlignmentFile(sample['path'], 'rb') as bamfile:
                # extract the qtl window
                for read in bamfile.fetch(qtl['chr'], qtl['start'], qtl['end']):
                    print read

            # pool the data
        # for each position in the qtl window
            # if the site is variable across the pool
                # store result
                # print info

    # pickle result

check_coverage()