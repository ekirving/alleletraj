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

def fetch_qtl_snps():

    try:
        # load the QTL SNPs data
        snps = pickle.load(open(COVERAGE_FILE, 'r'))

        print "INFO: Loaded data from '%s'" % COVERAGE_FILE

    except IOError:
        # file doesn't exist, so find the QTL SNPs
        snps = find_qtl_snps()

        # pickle.dump(snps, open(COVERAGE_FILE, 'w'))

    return snps


def find_qtl_snps():

    # open a db connection
    dbc = db_conn()

    # get a sorted list of unique QTLs
    qtls = dbc.get_records_sql(
        """SELECT GROUP_CONCAT(id) AS id, chromosome as chr, genomeLoc_start as start, genomeLoc_end as end
             FROM qtls
            WHERE genomeLoc_start IS NOT NULL
              AND genomeLoc_end IS NOT NULL
              AND (genomeLoc_end - genomeLoc_start) <= %s
              # AND id = 453
            GROUP BY chromosome, genomeLoc_start, genomeLoc_end
            ORDER BY chromosome, genomeLoc_start""" % MAX_QTL_LENGTH
    )

    print "INFO: Found %s unique QTLs" % len(qtls)

    # fetch the metadata
    df = fetch_metadata()

    print "INFO: Found %s samples to find SNPs in" % df.shape[0]

    # keep track of the SNPs we find
    snps = defaultdict(dict)

    # process each QTL
    for ids, qtl in qtls.iteritems():
        # get the window coordinates
        chrom, start, end = str(qtl['chr']), qtl['start'], qtl['end'] + 1

        # TODO if this is a very large window (i.e. > 100 Kb or 1 Mb) then chunk it

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
                                print "WARNING: chr%s:%s - skipping indel site" % (accession, chrom, pos)
                            continue

                        # skip low quality alignments
                        map_qual = pileupread.alignment.mapping_quality
                        if map_qual < MIN_MAPPING_QUAL:
                            if VERBOSE:
                                print "WARNING: chr%s:%s - skipping low mapq (%s)" % (accession, chrom, pos, map_qual)
                            continue

                        # get the read postion
                        read_pos = pileupread.query_position

                        # skip low quality base calls
                        base_qual = pileupread.alignment.query_qualities[read_pos]
                        if base_qual < MIN_BASE_QUAL:
                            if VERBOSE:
                                print "WARNING: chr%s:%s - skipping low baseq (%s)" % (accession, chrom, pos, base_qual)
                            continue

                        # get the overall length of the read
                        read_length = len(pileupread.alignment.query_sequence)

                        # soft clip bases near the edges of the read
                        if read_pos <= SOFT_CLIP_DIST or read_pos >= (read_length - SOFT_CLIP_DIST):
                            if VERBOSE:
                                print "WARNING: chr%s:%s - skipping soft clipped base (%s)" % (accession, chrom, pos, read_pos)
                            continue

                        # get the aligned base for this read
                        base = pileupread.alignment.query_sequence[read_pos]

                        # initialise the dictionary
                        if accession not in data[(chrom, pos)]:
                            data[(chrom, pos)][accession] = []

                        # add the base to the list
                        data[(chrom, pos)][accession].append(base)

        # now check the coverage for multi-sample SNPs
        for chrom, pos in data:

            # get the unique list of alleles at this locus
            alleles = set(itertools.chain.from_iterable(data[(chrom, pos)].values()))

            # TODO how should we handle single sample variants / are these informative?

            # is this a multi-sample polymorphic SNP
            if len(data[(chrom, pos)]) > 1 and len(alleles) > 1:

                # add the locus to the SNPs dictionary
                snps[ids][(chrom, pos)] = data[(chrom, pos)]

                print "SUCCESS: Found a SNP at chr%s:%s with %s samples and alleles %s" % \
                      (chrom, pos, len(data[(chrom, pos)]), list(alleles))

    return snps

fetch_qtl_snps()
