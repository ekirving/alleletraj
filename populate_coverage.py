#!/usr/bin/env python
# -*- coding: utf-8 -*-

from db_conn import db_conn

from collections import defaultdict

import pysam as ps
import multiprocessing as mp
import itertools
import mysql.connector
import socket

from natsort import natsorted

# observations about the QTL db
# - there are many overlapping QTL windows (some very large windows span multiple smaller ones)
# - there are many duplicate windows, ~8k (~44%)

# show lots of debugging output
VERBOSE = False

# minimum depth of coverage to call diploid genotypes
MIN_GENO_DEPTH = 10

# minumum mapping quality (hard filtered)
HARD_MAPQ_CUTOFF = 20

# minimum base quality (hard filtered)
HARD_BASEQ_CUTOFF = 20

# should we use multi-threading to speed up record insertion
MULTI_THREADED = True if socket.gethostname() != 'macbookpro.local' else False

# the minimum minor allele frequency of modern SNPs to include
MIN_MAF = 0.05

# no single worker should use more than 30% of the available cores
MAX_CPU_CORES = int(mp.cpu_count() * 0.3)

# TODO ...
EUROPE = ['Armenia', 'Belgium', 'Bulgaria', 'Croatia', 'Czech Rep.', 'Denmark', 'England', 'Faroes', 'France', 'Georgia',
          'Germany', 'Greece', 'Hungary', 'Iceland', 'Italy', 'Moldova', 'Netherlands', 'Poland', 'Portugal', 'Romania',
          'Serbia', 'Slovakia', 'Spain', 'Sweden', 'Switzerland', 'Ukraine']


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


def populate_intervals(species):

    # open a db connection
    dbc = db_conn()

    # tidy up any existing intervals
    dbc.delete_records('intervals', {'species': species})

    # get all the unique QTL windows
    results = dbc.get_records_sql("""
        SELECT q.chromosome AS chrom, q.start, q.end
          FROM qtls q
         WHERE q.species = '{}'
           AND q.valid = 1
      GROUP BY q.chromosome, q.start, q.end
    """.format(species), key=None
    )

    intervals = defaultdict(list)

    for result in results:
        # group the intervals by chromosome
        intervals[result['chrom']].append((result['start'], result['end']))

    num_sites = 0
    num_intvals = 0

    for chrom in natsorted(intervals.keys()):
        # merge overlapping intervals
        intervals[chrom] = list(merge_intervals(intervals[chrom]))

        # count the resulting merged intervals
        num_intvals += len(intervals[chrom])

        for start, end in intervals[chrom]:
            # sum the number of sites across all intervals
            num_sites += end - start

            # add the merged intervals to the db
            record = {'species': species, 'chrom': chrom, 'start': start, 'end': end }
            dbc.save_record('intervals', record)

    print "INFO: Extracted {:,} {} intervals, totalling {:,} bp".format(num_intvals, species, num_sites)


def populate_interval_snps(species):
    """
    Now we have ascertained all the modern SNPs, let's find those that intersect with the unique intervals.
    """

    dbc = db_conn()

    print "INFO: Populating all the {} interval SNPs".format(species)

    # insert linking records to make future queries much quicker
    dbc.execute_sql("""
        INSERT INTO intervals_snps (interval_id, modsnp_id)
             SELECT i.id, ms.id
               FROM intervals i
               JOIN modern_snps ms 
                 ON ms.species = i.species
                AND ms.chrom = i.chrom
                AND ms.site BETWEEN i.start AND i.end
              WHERE i.species = '%s'
                AND ms.maf >= %s""" % (species, MIN_MAF)
    )


def populate_coverage(species):
    """
    Scan all the samples for coverage of the QTLs and save the results to the database.
    """

    # open a db connection
    dbc = db_conn()

    # count all the intervals we've not finished processing yet
    intervals = dbc.get_records('intervals', {'species': species, 'finished': 0})

    # make a list of permissible countries
    countries = "','".join(EUROPE)

    # get all the samples w/ BAM files
    samples = dbc.get_records_sql(
        """SELECT *
             FROM samples
            WHERE species = '%s'
              AND country IN ('%s')
              AND path IS NOT NULL""" % (species, countries)
    )

    print "INFO: Processing {:,} intervals in {:,} {} samples".format(len(intervals), len(samples), species)

    # before we start, tidy up any records from intervals that were not finished
    dbc.execute_sql("""
        DELETE sample_reads
          FROM sample_reads
          JOIN intervals 
            ON intervals.id = sample_reads.intervalID
         WHERE intervals.species = '%s'
           AND intervals.finished = 0""" % species
    )

    if MULTI_THREADED:
        # process the chromosomes with multi-threading to make this faster
        pool = mp.Pool(MAX_CPU_CORES)
        pool.map(process_interval, itertools.izip(intervals.values(), itertools.repeat(samples)))
    else:
        # process the chromosomes without multi-threading
        for interval in intervals.values():
            process_interval((interval, samples))

    print "FINISHED: Fully populated all the %s samples for %s intervals" % (species, len(intervals))


def process_interval(args):
    """
    Scan all the intervals across this chromosome and add the covered bases to the DB, as long as they are variable in
    the modern data.
    """

    # extract the nested tuple of arguments (an artifact of using izip to pass args to mp.Pool)
    (interval, samples) = args

    # open a db connection
    dbc = db_conn()

    # unpack the interval
    interval_id, chrom, start, end = interval['id'], interval['chrom'], interval['start'], interval['end']

    # get all the modern SNPs in this interval
    snps = dbc.get_records_sql("""
        SELECT ms.site
          FROM intervals i
          JOIN intervals_snps s
            ON s.interval_id = i.id
          JOIN modern_snps ms
            ON ms.id = s.modsnp_id
         WHERE i.id = %s""" % interval_id, key='site'
    ).keys()

    print "INFO: Scanning interval chr{}:{}-{} for {:,} SNPs".format(chrom, start, end, len(snps))

    # the column headers for batch inserting into the db
    fields = ('intervalID', 'sampleID', 'chrom', 'pos', 'base', 'mapq', 'baseq', 'dist')

    num_reads = 0

    # check all the samples for coverage in this interval
    for sample_id, sample in samples.iteritems():

        if VERBOSE:
            print "INFO: Scanning interval chr{}:{}-{} in sample {}".format(chrom, start, end, sample['accession'])

        reads = list()

        # open the BAM file for reading
        with ps.AlignmentFile(sample['path'], 'rb') as bamfile:

            # get the full interval
            for pileupcolumn in bamfile.pileup(chrom, start, end + 1):

                # what position are we at in the BAM file
                pos = pileupcolumn.reference_pos

                # skip all non-SNP sites
                if pos not in snps:
                    continue

                if len(pileupcolumn.pileups) >= MIN_GENO_DEPTH:
                    # TODO genotype this site properly
                    pass

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

                    if mapq < HARD_MAPQ_CUTOFF:
                        continue

                    # get the base quality
                    baseq = pileupread.alignment.query_qualities[read_pos]

                    if baseq < HARD_BASEQ_CUTOFF:
                        continue

                    # get the overall length of the read
                    read_length = len(pileupread.alignment.query_sequence)

                    # how close is the base to the edge of the read
                    dist = min(read_pos, read_length - read_pos)

                    # setup the record to insert, in this order
                    read = (interval_id, sample_id, chrom, pos, base, mapq, baseq, dist)

                    # store the read so we can batch insert later
                    reads.append(read)

        if reads:
            # count the total number of reads
            num_reads += len(reads)

            # bulk insert all the reads for this sample
            dbc.save_records('sample_reads', fields, reads)

    # update the interval to show we're done
    interval['finished'] = 1
    dbc.save_record('intervals', interval)

    print "INFO: Found {:,} reads for interval chr{}:{}-{}".format(num_reads, chrom, start, end)
