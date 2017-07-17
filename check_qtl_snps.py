#!/usr/bin/env python
# -*- coding: utf-8 -*-

from db_conn import db_conn
from fetch_metadata import fetch_metadata

from collections import defaultdict
from pprint import pprint

import pysam as ps
import random

# observations about the QTL db
# - there are many overlapping QTL windows (some very large windows span multiple smaller ones)
# - there are many duplicate windows, ~8k (~44%)

VERBOSE = False

# maximum length of the QTL to process (100 Kb)
MAX_QTL_LENGTH = 100000


def merge_intervals(ranges):
    """
    Merge overlapping intervals, so we only check each site once
    """
    saved = list(ranges[0])
    for start, end in sorted([sorted(t) for t in ranges]):
        if start <= saved[1]:
            saved[1] = max(saved[1], end)
        else:
            yield tuple(saved)
            saved[0] = start
            saved[1] = end

    yield tuple(saved)


def populate_sample_coverage():
    """
    Extract the list of unique QTL intervals, scan all the samples for coverage and save the results to the databse.
    """

    # open a db connection
    dbc = db_conn()

    # get all the QTL intervals
    results = dbc.get_records_sql(
        """
        SELECT DISTINCT chromosome as chrom, genomeLoc_start as start, genomeLoc_end as end
          FROM qtls
         WHERE (genomeLoc_end - genomeLoc_start) <= %s
      ORDER BY chrom, start, end
        """ % MAX_QTL_LENGTH, key=None
    )

    intvals = defaultdict(list)

    # group the intervals by chromosome
    for result in results:
        intvals[result['chrom']].append((result['start'], result['end']))

    num_sites = 0
    num_intvals = 0

    for chrom in intvals:
        # merge the intervals
        intvals[chrom] = list(merge_intervals(intvals[chrom]))

        # count the intervals
        num_intvals += len(intvals[chrom])

        # count the sites
        num_sites += sum([intval[1]-intval[0] for intval in intvals[chrom]])

    print "INFO: Found {:,} intervals to scan, totalling {:,} bp".format(num_intvals, num_sites)

    # TODO fetch the metadata
    df = fetch_metadata()

    print "INFO: Found %s samples to extract coverage for" % df.shape[0]

    # process each interval
    for chrom in intvals:

        print "INFO: Processing chromosome %s" % chrom

        for start, end in intvals[chrom]:

            if VERBOSE:
                print "INFO: Checking interval chr%s:%s-%s" % (chrom, start, end)

            # check all the samples for coverage in this interval
            for accession, sample in df.iterrows():

                if VERBOSE:
                    print "INFO: Checking sample %s" % accession

                # open the BAM file for reading
                with ps.AlignmentFile(sample['path'], 'rb') as bamfile:

                    # extract the qtl window
                    for pileupcolumn in bamfile.pileup(str(chrom), start, end+1):

                        # what position are we at in the BAM file
                        pos = pileupcolumn.reference_pos

                        reads = []

                        # iterate over all the reads for this site
                        for pileupread in pileupcolumn.pileups:

                            # skip alignments that don't have a base at this site (i.e. indels)
                            if pileupread.is_del or pileupread.is_refskip:
                                continue

                            # setup the record to insert
                            read = dict()
                            read['sampleID'] = accession
                            read['chrom'] = chrom
                            read['pos'] = pos

                            # get the read position
                            read_pos = pileupread.query_position

                            # get the aligned base for this read
                            read['base'] = pileupread.alignment.query_sequence[read_pos]

                            # get the map quality
                            read['mapq'] = pileupread.alignment.mapping_quality

                            # get the base quality
                            read['baseq'] = pileupread.alignment.query_qualities[read_pos]

                            # get the overall length of the read
                            read_length = len(pileupread.alignment.query_sequence)

                            # how close is the base to the edge of the read
                            read['dist'] = min(read_pos, read_length-read_pos)

                            reads.append(read)

                        if reads:
                            # choose one read at random
                            read = random.choice(reads)
                            read['random'] = 1

                        # save all the reads to the db
                        for read in reads:
                            dbc.save_record('sample_reads', read)

    print "FINISHED: Fully populated the database"

populate_sample_coverage()
