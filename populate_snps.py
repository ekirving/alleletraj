#!/usr/bin/env python
# -*- coding: utf-8 -*-

from db_conn import db_conn

from collections import defaultdict

import pysam as ps
import multiprocessing as mp
import itertools

from natsort import natsorted

# observations about the QTL db
# - there are many overlapping QTL windows (some very large windows span multiple smaller ones)
# - there are many duplicate windows, ~8k (~44%)

# show lots of debugging output
VERBOSE = False

# maximum length of the QTL to process (100 Kb)
MAX_QTL_LENGTH = 100000

# should we use multi-threading to speed up record insertion
MULTI_THREADED_SEARCH = True

# no single worker should use more than 30% of the available cores
MAX_CPU_CORES = int(mp.cpu_count() * 0.3)


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


def populate_snps():
    """
    Extract the list of unique QTL intervals, scan all the samples for coverage and save the results to the database.
    """

    # open a db connection
    dbc = db_conn()

    # get all the QTL intervals
    results = dbc.get_records_sql("""
        SELECT DISTINCT chromosome as chrom, genomeLoc_start as start, genomeLoc_end as end
          FROM qtls
         WHERE (genomeLoc_end - genomeLoc_start) <= %s
      ORDER BY chrom, start, end""" % MAX_QTL_LENGTH, key=None
    )

    intervals = defaultdict(list)

    # group the intervals by chromosome
    for result in results:
        intervals[result['chrom']].append((result['start'], result['end']))

    num_sites = 0
    num_intvals = 0

    for chrom in intervals:
        # merge the intervals
        intervals[chrom] = list(merge_intervals(intervals[chrom]))

        # count the intervals
        num_intvals += len(intervals[chrom])

        # count the sites across all intervals
        num_sites += sum([intval[1]-intval[0] for intval in intervals[chrom]])

    print "INFO: Found {:,} intervals to scan, totalling {:,} bp".format(num_intvals, num_sites)

    # get all the samples w/ BAM files
    samples = dbc.get_records_sql("SELECT * FROM samples "
                                  "WHERE path IS NOT NULL")

    print "INFO: Found %s samples to extract coverage for" % len(samples)

    # process the chromosomes in order
    chroms = natsorted(intervals.keys())

    if MULTI_THREADED_SEARCH:
        # process the chromosomes with multi-threading to make this faster
        pool = mp.Pool(MAX_CPU_CORES)
        pool.map(process_chrom, itertools.izip(chroms, itertools.repeat((intervals, samples))))
    else:
        # process the chromosomes without multi-threading
        for chrom in chroms:
            process_chrom((chrom, (intervals, samples)))

    print "FINISHED: Fully populated all the QTL coverage"


def process_chrom(args):
    """
    Scan all the intervals across this chromosome and add the covered bases to the DB.
    """

    # extract the nested tuple of arguments (an artifact of using izip to pass args to mp.Pool)
    chrom, (intervals, samples) = args

    # open a db connection
    dbc = db_conn()

    # we're going to do a lot of inserts, so let's speed things up a little
    dbc.cursor.execute("SET unique_checks=0")
    dbc.cursor.execute("SET foreign_key_checks=0")

    print "INFO: Processing chromosome %s (%s intervals)" % (chrom, len(intervals[chrom]))

    # process each interval
    for start, end in intervals[chrom]:

        print "INFO: Checking interval chr%s:%s-%s" % (chrom, start, end)

        # check all the samples for coverage in this interval
        for sample_id, sample in samples.iteritems():

            if VERBOSE:
                print "INFO: Checking sample %s" % sample['accession']

            # open the BAM file for reading
            with ps.AlignmentFile(sample['path'], 'rb') as bamfile:

                reads = []

                # extract the qtl window
                for pileupcolumn in bamfile.pileup(str(chrom), start, end + 1):

                    # what position are we at in the BAM file
                    pos = pileupcolumn.reference_pos

                    # iterate over all the reads for this site
                    for pileupread in pileupcolumn.pileups:

                        # skip alignments that don't have a base at this site (i.e. indels)
                        if pileupread.is_del or pileupread.is_refskip:
                            continue

                        # setup the record to insert
                        read = dict()
                        read['sampleID'] = sample_id
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
                        read['dist'] = min(read_pos, read_length - read_pos)

                        # store the read so we can batch insert later
                        reads.append(read)

                if reads:
                    dbc.save_records('sample_reads', reads)

    # now lets put it all back the way it was before
    dbc.cursor.execute("SET unique_checks=1")
    dbc.cursor.execute("SET foreign_key_checks=1")
