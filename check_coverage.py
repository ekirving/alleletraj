#!/usr/bin/env python
# -*- coding: utf-8 -*-

from db_conn import db_conn
from fetch_metadata import fetch_metadata

from collections import defaultdict
from pprint import pprint

import itertools
import pysam as ps

import pickle

# observations about the QTL db
# - there are many overlapping QTL windows (some very large windows span multiple smaller ones)
# - there are many duplicate windows, ~8k (~44%)

# the minimum phred scaled genotype quality (30 = 99.9%)
MIN_BASE_QUAL = 30
MIN_MAPPING_QUAL = 30

# number of bases to soft clip
SOFT_CLIP_DIST = 5

VERBOSE = False

COVERAGE_FILE = 'coverage-pigs.pickle'

# maximum length of the QTL to process (100 kb)
MAX_QTL_LENGTH = 100000

def fetch_coverage():

    try:
        # load the coverage data
        data = pickle.load(open(COVERAGE_FILE, 'r'))

        print "INFO: Loaded data from '%s'" % COVERAGE_FILE

    except IOError:
        # file doesn't exist, so compute the coverage
        data = check_coverage()

        pickle.dump(data, open(COVERAGE_FILE, 'w'))

    return data


def check_coverage():

    # open a db connection
    dbc = db_conn()

    # get a sorted list of unique QTLs
    qtls = dbc.get_records_sql(
        """SELECT GROUP_CONCAT(id) AS id, chromosome as chr, genomeLoc_start as start, genomeLoc_end as end
             FROM qtls
            WHERE genomeLoc_start IS NOT NULL
              AND genomeLoc_end IS NOT NULL
              AND (genomeLoc_end - genomeLoc_start) = %s
            GROUP BY chromosome, genomeLoc_start, genomeLoc_end
            ORDER BY chromosome, genomeLoc_start""" % MAX_QTL_LENGTH
    )

    print "INFO: Found %s unique QTLs" % len(qtls)

    # fetch the metadata
    df = fetch_metadata()

    print "INFO: Found %s samples to compute coverage for" % df.shape[0]

    # process each QTL
    for ids, qtl in qtls.iteritems():
        # get the window coordinates
        chrom, start, end = str(qtl['chr']), qtl['start'], qtl['end'] + 1

        print "INFO: Checking QTLs with window chr%s:%s-%s (%s)" % (chrom, start, end-1, ids)

        data = defaultdict(dict)

        # check all the samples for coverage in this QTL
        for accession, sample in df.iterrows():

            if VERBOSE:
                print "INFO: Checking sample %s" % accession

            # open the BAM file for reading
            with ps.AlignmentFile(sample['path'], 'rb') as bamfile:

                # extract the qtl window
                for pileupcolumn in bamfile.pileup(chrom, start, end):

                    # what position are we at in the BAM file
                    pos = pileupcolumn.reference_pos

                    # iterate over all the reads for this site
                    for pileupread in pileupcolumn.pileups:

                        # skip alignments that don't have a base at this site (i.e. indels)
                        if pileupread.is_del or pileupread.is_refskip:
                            if VERBOSE:
                                print "WARNING: chr%s:%s - skipping indel site" % (chrom, pos)
                            continue

                        # skip low quality alignments
                        map_qual = pileupread.alignment.mapping_quality
                        if map_qual < MIN_MAPPING_QUAL:
                            if VERBOSE:
                                print "WARNING: chr%s:%s - skipping low mapq (%s)" % (chrom, pos, map_qual)
                            continue

                        # get the read postion
                        read_pos = pileupread.query_position

                        # skip low quality base calls
                        base_qual = pileupread.alignment.query_qualities[read_pos]
                        if base_qual < MIN_BASE_QUAL:
                            if VERBOSE:
                                print "WARNING: chr%s:%s - skipping low baseq (%s)" % (chrom, pos, base_qual)
                            continue

                        # get the overall length of the read
                        read_length = len(pileupread.alignment.query_sequence)

                        # soft clip bases near the edges of the read
                        if read_pos <= SOFT_CLIP_DIST or read_pos >= (read_length - SOFT_CLIP_DIST):
                            if VERBOSE:
                                print "WARNING: chr%s:%s - skipping soft clipped base (%s)" % (chrom, pos, read_pos)
                            continue

                        # get the aligned base for this read
                        base = pileupread.alignment.query_sequence[read_pos]

                        # initialise the dictionary
                        if accession not in data[(chrom, pos)]:
                            data[(chrom, pos)] = {accession: []}

                        # add the base to the list
                        data[(chrom, pos)][accession].append(base)

        # now check what overlap we have between the samples
        for chrom, pos in data:

            # get the unique list of alleles at this locus
            alleles = set(itertools.chain.from_iterable(data[(chrom, pos)].values()))

            # TODO single sample variants?
            if len(alleles) > 1 and len(data[(chrom, pos)]) > 1:
                print "SUCCESS: Found a SNP at chr%s:%s with %s samples and alleles %s" % \
                      (chrom, pos, len(data[(chrom, pos)]), list(alleles))

fetch_coverage()
