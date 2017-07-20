#!/usr/bin/env python
# -*- coding: utf-8 -*-

from db_conn import db_conn

from collections import defaultdict

import pysam as ps
import multiprocessing as mp
import itertools
import mysql.connector

from natsort import natsorted

# observations about the QTL db
# - there are many overlapping QTL windows (some very large windows span multiple smaller ones)
# - there are many duplicate windows, ~8k (~44%)

# show lots of debugging output
VERBOSE = False

# maximum length of the QTL to process (100 Kb)
MAX_QTL_LENGTH = 100000

# should we use multi-threading to speed up record insertion
MULTI_THREADED = False

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


def populate_intervals():

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

    for chrom in natsorted(intervals.keys()):
        # merge the intervals
        intervals[chrom] = list(merge_intervals(intervals[chrom]))

        # count the intervals
        num_intvals += len(intervals[chrom])

        for start, end in intervals[chrom]:
            # count the sites across all intervals
            num_sites += end - start

            # add the intervals to the db
            record = {'chrom': chrom, 'start': start, 'end': end }
            dbc.save_record('intervals', record)

    print "INFO: Found {:,} intervals, totalling {:,} bp".format(num_intvals, num_sites)


def populate_coverage():
    """
    Extract the list of unique QTL intervals, scan all the samples for coverage and save the results to the database.
    """

    # open a db connection
    dbc = db_conn()

    # get all the intervals we've not finished processing yet
    intervals = dbc.get_records('intervals', {'finished': 0})

    print "INFO: Found {:,} intervals to scan".format(len(intervals))

    # get all the samples w/ BAM files
    samples = dbc.get_records_sql(
        "SELECT * FROM samples "
        "WHERE path IS NOT NULL"
    )

    print "INFO: Found {:,} samples to extract coverage for".format(len(samples))

    # before we start, tidy up any records from intervals that were not finished
    dbc.cursor.execute("""
        DELETE sample_reads
          FROM sample_reads
          JOIN intervals 
            ON intervals.id = sample_reads.intervalID
         WHERE intervals.finished != 1"""
    )
    dbc.cnx.commit()

    if MULTI_THREADED:
        # process the chromosomes with multi-threading to make this faster
        pool = mp.Pool(MAX_CPU_CORES)
        pool.map(process_interval, itertools.izip(intervals.values(), itertools.repeat(samples)))
    else:
        # process the chromosomes without multi-threading
        for interval in intervals.values():
            process_interval((interval, samples))

    # add the indexs back in now that we're done
    dbc.cursor.execute("ALTER TABLE sample_reads ADD INDEX (chrom, pos)")
    dbc.cnx.commit()

    print "FINISHED: Fully populated all the QTL coverage"


def process_interval(args):
    """
    Scan all the intervals across this chromosome and add the covered bases to the DB.
    """

    # extract the nested tuple of arguments (an artifact of using izip to pass args to mp.Pool)
    (interval, samples) = args

    # open a db connection
    dbc = db_conn()

    interval_id, chrom, start, end = interval['id'], interval['chrom'], interval['start'], interval['end']

    print "INFO: Checking interval chr%s:%s-%s" % (chrom, start, end)

    pos_alleles = defaultdict(set)
    pos_samples = defaultdict(set)
    reads = defaultdict(list)

    # check all the samples for coverage in this interval
    for sample_id, sample in samples.iteritems():

        if VERBOSE:
            print "INFO: Checking sample %s" % sample['accession']

        # open the BAM file for reading
        with ps.AlignmentFile(sample['path'], 'rb') as bamfile:

            # extract the qtl window
            for pileupcolumn in bamfile.pileup(str(chrom), start, end + 1):

                # what position are we at in the BAM file
                pos = pileupcolumn.reference_pos

                if pos < start or pos > end:
                    # skip columns not in range (happens because of overlapping reads
                    continue

                # iterate over all the reads for this site
                for pileupread in pileupcolumn.pileups:

                    # skip alignments that don't have a base at this site (i.e. indels)
                    if pileupread.is_del or pileupread.is_refskip:
                        continue

                    # get the read position
                    read_pos = pileupread.query_position

                    # get the aligned base for this read
                    base = pileupread.alignment.query_sequence[read_pos]

                    # get the map quality
                    mapq = pileupread.alignment.mapping_quality

                    # get the base quality
                    baseq = pileupread.alignment.query_qualities[read_pos]

                    # get the overall length of the read
                    read_length = len(pileupread.alignment.query_sequence)

                    # how close is the base to the edge of the read
                    dist = min(read_pos, read_length - read_pos)

                    # add this allele and sample to the sets
                    pos_alleles[pos].add(base)
                    pos_samples[pos].add(sample_id)

                    # setup the record to insert, in this order
                    read = (interval_id, sample_id, chrom, pos, base, mapq, baseq, dist)

                    # store the read so we can batch insert later
                    reads[pos].append(read)

    # check which sites have more than one allele and sample
    insert = [reads[pos] for pos in pos_alleles if len(pos_alleles[pos]) > 1 and len(pos_samples[pos]) > 1]

    if insert:
        # the column headers for the list of tuples
        fields = ('intervalID', 'sampleID', 'chrom', 'pos', 'base', 'mapq', 'baseq', 'dist')

        # collapse the list of lists
        insert = itertools.chain.from_iterable(insert)

        # bulk insert the records
        dbc.save_records('sample_reads', fields, insert)

    # update the interval to show we're done
    interval['finished'] = 1
    dbc.save_record('intervals', interval)


