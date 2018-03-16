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
VERBOSE = True

# maximum length of the QTL to process (100 Kb)
MAX_QTL_LENGTH = 100000

# minimum depth of coverage to call diploid genotypes
MIN_GENO_DEPTH = 10

# minumum mapping quality (hard filtered)
HARD_MAPQ_CUTOFF = 20

# minimum base quality (hard filtered)
HARD_BASEQ_CUTOFF = 20

# should we use multi-threading to speed up record insertion
MULTI_THREADED = True if socket.gethostname() != 'macbookpro.local' else False

# no single worker should use more than 30% of the available cores
MAX_CPU_CORES = int(mp.cpu_count() * 0.3)

# sizes of each chrom in the given assemblies
CHROM_SIZE = {

    # UMD_3.1.1
    'cattle': {'1': 158337067, '2': 137060424, '3': 121430405, '4': 120829699, '5': 121191424, '6': 119458736,
               '7': 112638659, '8': 113384836, '9': 105708250, '10': 104305016, '11': 107310763, '12': 91163125,
               '13': 84240350, '14': 84648390, '15': 85296676, '16': 81724687, '17': 75158596, '18': 66004023,
               '19': 64057457, '20': 72042655, '21': 71599096, '22': 61435874, '23': 52530062, '24': 62714930,
               '25': 42904170, '26': 51681464, '27': 45407902, '28': 46312546, '29': 51505224, 'X': 148823899},

    # Sscrofa11.1
    'pig':    {'1': 274330532, '2': 151935994, '3': 132848913, '4': 130910915, '5': 104526007, '6': 170843587,
               '7': 121844099, '8': 138966237, '9': 139512083, '10': 69359453, '11': 79169978, '12': 61602749,
               '13': 208334590, '14': 141755446, '15': 140412725, '16': 79944280, '17': 63494081, '18': 55982971,
               'X': 125939595, 'Y': 43547828},

    # EquCab2.0
    'horse':  {'1': 185838109, '2': 120857687, '3': 119479920, '4': 108569075, '5': 99680356, '6': 84719076,
               '7': 98542428, '8': 94057673, '9': 83561422, '10': 83980604, '11': 61308211, '12': 33091231,
               '13': 42578167, '14': 93904894, '15': 91571448, '16': 87365405, '17': 80757907, '18': 82527541,
               '19': 59975221, '20': 64166202, '21': 57723302, '22': 49946797, '23': 55726280, '24': 46749900,
               '25': 39536964, '26': 41866177, '27': 39960074, '28': 46177339, '29': 33672925, '30': 30062385,
               '31': 24984650, 'X': 124114077},

    # CHIR_1.0
    'goat':   {'1': 155011307, '2': 135415751, '3': 116796116, '4': 115961478, '5': 111055201, '6': 114334461,
               '7': 106547263, '8': 111020524, '9': 90293942, '10': 99198151, '11': 105305070, '12': 82535142,
               '13': 80625018, '14': 92306894, '15': 78986926, '16': 77678508, '17': 71877645, '18': 61067880,
               '19': 62130014, '20': 71279863, '21': 66773250, '22': 57956300, '23': 49403180, '24': 61756751,
               '25': 41496684, '26': 50169583, '27': 44118849, '28': 43231948, '29': 48376377, 'X': 121952644},
}


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
    dbc.delete_records('intervals', {'species':species})

    # get all the unique QTL peaks
    results = dbc.get_records_sql("""
        SELECT d.*
          FROM qtls q
          JOIN dbsnp d on d.rsnumber = q.peak
         WHERE q.species = '{}' 
           AND q.peak like 'rs%'                  # only QTLs with a dbsnp peak
           AND q.associationType = 'Association'  # only GWAS studies
           AND q.significance = 'Significant'     # only significant hits
           AND IFNULL(q.genomeLoc_end - q.genomeLoc_start, 0) <= {}
      GROUP BY d.id
    """.format(species, MAX_QTL_LENGTH), key=None
    )

    intervals = defaultdict(list)

    # set offset to be half the max QTL length
    offset = MAX_QTL_LENGTH / 2

    for result in results:
        # get the size of the current chrom
        chom_size = CHROM_SIZE[species][result['chrom']]

        # calculate the bounded window size
        start = result['site'] - offset if result['site'] > offset else 0
        end   = result['site'] + offset if result['site'] + offset < chom_size else chom_size

        # group the intervals by chromosome
        intervals[result['chrom']].append((start, end))

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


def populate_coverage(species):
    """
    Scan all the samples for coverage of the QTLs and save the results to the database.
    """

    # open a db connection
    dbc = db_conn()

    # count all the intervals we've not finished processing yet
    intervals = dbc.get_records('intervals', {'species': species, 'finished': 0})

    # TODO this needs better filtering (location, domestic status, etc)
    # get all the samples w/ BAM files
    samples = dbc.get_records_sql(
        """SELECT *
             FROM samples
            WHERE species = '%s'
              AND path IS NOT NULL""" % species
    )

    print "INFO: Processing {:,} intervals in {:,} {} samples".format(len(intervals), len(samples), species)

    # before we start, tidy up any records from intervals that were not finished
    dbc.cursor.execute("""
        DELETE sample_reads
          FROM sample_reads
          JOIN intervals 
            ON intervals.id = sample_reads.intervalID
         WHERE intervals.species = '%s'
           AND intervals.finished = 0""" % species
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
          JOIN modern_snps ms 
            ON ms.species = i.species
           AND ms.chrom = i.chrom
           AND ms.site between i.start and i.end
         WHERE i.id = %s""" % interval_id, key='site'
    ).keys()

    print "INFO: Scanning interval chr{}:{}-{} for {:,} SNPs".format(chrom, start, end, len(snps))

    # the column headers for the list of tuples
    fields = ('intervalID', 'sampleID', 'chrom', 'pos', 'base', 'mapq', 'baseq', 'dist')

    num_reads = 0

    # check all the samples for coverage in this interval
    for sample_id, sample in samples.iteritems():

        if VERBOSE:
            print "INFO: Scanning interval chr{}:{}-{} for sample {}".format(chrom, start, end, sample['accession'])

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
