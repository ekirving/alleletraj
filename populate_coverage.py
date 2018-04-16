#!/usr/bin/env python
# -*- coding: utf-8 -*-

from db_conn import db_conn
from discover_snps import *

from collections import defaultdict

import pysam as ps
import multiprocessing as mp
import traceback
import sys
import itertools
import socket
import subprocess

from natsort import natsorted

# observations about the QTL db
# - there are many overlapping QTL windows (some very large windows span multiple smaller ones)
# - there are many duplicate windows, ~8k (~44%)

# show lots of debugging output
VERBOSE = False

# minimum depth of coverage to call diploid genotypes
MIN_GENO_DEPTH = 10

# TODO what about the others species
# location of reference genome
REF_FILE = "fasta/Sus_scrofa.Sscrofa10.2.dna.toplevel.fa"

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


def run_cmd(cmd, shell=False):
    """
    Executes the given command in a system subprocess

    :param cmd: The system command to run (list|string)
    :param shell: Use the native shell
    :return: The stdout stream
    """
    # subprocess only accepts strings
    cmd = [str(args) for args in cmd]

    # run the command
    proc = subprocess.Popen(cmd,
                            shell=shell,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)

    # fetch the output and error
    (stdout, stderr) = proc.communicate()

    # bail if something went wrong
    if proc.returncode:
        raise Exception(stderr)

    return stdout


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
        SELECT q.chrom AS chrom, q.start, q.end
          FROM qtls q
         WHERE q.species = '{species}'
           AND q.valid = 1
      GROUP BY q.chrom, q.start, q.end""".format(species=species), key=None)

    intervals = defaultdict(list)

    for result in results:
        # group the intervals by chrom
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
            record = {'species': species, 'chrom': chrom, 'start': start, 'end': end}
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
              WHERE i.species = '{species}'
                AND ms.maf >= {minmaf}""".format(species=species, minmaf=MIN_MAF))


def populate_coverage(species):
    """
    Scan all the samples for coverage of the QTLs and save the results to the database.
    """

    # open a db connection
    dbc = db_conn()

    # count all the intervals we've not finished processing yet
    intervals = dbc.get_records('intervals', {'species': species, 'finished': 0})

    # get all the valid samples
    samples = dbc.get_records_sql(
        """SELECT s.*, GROUP_CONCAT(sf.path) paths
             FROM samples s
             JOIN sample_files sf
               ON sf.sample_id = s.id
            WHERE s.species = '{species}'
              AND s.valid = 1
         GROUP BY s.id;""".format(species=species))

    print "INFO: Processing {:,} intervals in {:,} {} samples".format(len(intervals), len(samples), species)

    # before we start, tidy up any records from intervals that were not finished
    dbc.execute_sql("""
        DELETE sample_reads
          FROM sample_reads
          JOIN intervals 
            ON intervals.id = sample_reads.interval_id
         WHERE intervals.species = '{species}'
           AND intervals.finished = 0""".format(species=species))

    if MULTI_THREADED:
        # process the chromosomes with multi-threading to make this faster
        pool = mp.Pool(MAX_CPU_CORES)
        pool.map(process_interval, itertools.izip(intervals.values(), itertools.repeat(samples)))
    else:
        # process the chromosomes without multi-threading
        for interval in intervals.values():
            process_interval((interval, samples))

    print "FINISHED: Fully populated all the {} samples for {:,} intervals".format(species, len(intervals))


def process_interval(args):
    """
    Scan all the intervals across this chromosome and add the covered bases to the DB, as long as they are variable in
    the modern data.
    """

    try:

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
             WHERE i.id = {id}""".format(id=interval_id), key='site').keys()

        print "INFO: Scanning interval chr{}:{}-{} for {:,} SNPs".format(chrom, start, end, len(snps))

        # the column headers for batch inserting into the db
        fields = ('interval_id', 'sample_id', 'chrom', 'site', 'base', 'mapq', 'baseq', 'dist')

        num_reads = 0

        # check all the samples for coverage in this interval
        for sample_id, sample in samples.iteritems():

            if VERBOSE:
                print "INFO: Scanning interval chr{}:{}-{} in sample {}".format(chrom, start, end, sample['accession'])

            # buffer the reads so we can bulk insert them into the db
            reads = defaultdict(list)

            # there may be multiple BAM files for each sample
            for path in sample['paths'].splt(','):

                # open the BAM file for reading
                with ps.AlignmentFile(path, 'rb') as bamfile:

                    # get the full interval
                    for pileupcolumn in bamfile.pileup(chrom, start, end + 1):

                        # IMPORTANT `reference_pos` is 0 based !!!!
                        # see http://pysam.readthedocs.io/en/latest/api.html#pysam.PileupColumn.reference_pos
                        site = pileupcolumn.reference_pos + 1

                        # skip all non-SNP sites
                        if site not in snps:
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
                            read = (interval_id, sample_id, chrom, site, base, mapq, baseq, dist)

                            # store the read so we can batch insert later
                            reads[(chrom, site)].append(read)

            # now we're buffered all the reads, lets call diploid genotypes on those that pass our depth threshold
            diploid = [idx for idx in reads if len(reads[idx]) >= MIN_GENO_DEPTH]

            if diploid:

                print "INFO: Calling diploid bases in {:,} sites for sample {}".format(len(diploid), sample_id)

                pos_file = 'vcf/diploid-int{}-sample{}.tsv'.format(interval_id, sample_id)
                vcf_file = 'vcf/diploid-int{}-sample{}.vcf'.format(interval_id, sample_id)

                # save all the callable positions to a file
                with open(pos_file, 'w') as fout:
                    fout.write("\n".join("{}\t{}".format(chrom, site) for (chrom, site) in diploid))

                # restrict the callable region using the interval start and end
                region = "{}:{}-{}".format(chrom, start, end)

                # use all the BAM files
                bam_files = " ".join(sample['path'].split(','))

                params = {'region': region, 'targets': pos_file, 'ref': REF_FILE, 'bams': bam_files, 'vcf': vcf_file}

                # call bases with bcftools (and drop indels and other junk)
                # uses both --region (random access) and --targets (streaming) for optimal speed
                # see https://samtools.github.io/bcftools/bcftools.html#mpileup
                cmd = "bcftools mpileup --region {region} --targets-file {targets} --fasta-ref {ref} {bams} " \
                      " | bcftools call --multiallelic-caller --output-type v " \
                      " | bcftools view --exclude-types indels,bnd,other --exclude INFO/INDEL=1 --output-file {vcf}".format(**params)

                # run the base calling
                run_cmd([cmd], shell=True)

                # parse the results with pysam
                for rec in ps.VariantFile(vcf_file).fetch():

                    # get the genotype call for this site
                    geno = rec.samples[sample['accession']]['GT']

                    if geno >= MIN_GENO_QUAL:
                        # the genotype is good, so drop the raw reads
                        reads.pop((rec.chrom, rec.pos))

                    # decode the GT notation into allele calls (e.g. 0/0, 0/1, 1/1)
                    alleles = [rec.alleles[idx] for idx in geno if idx is not None]

                    for allele in alleles:
                        # compose the read records
                        read = {
                            'interval_id': interval_id,
                            'sample_id': sample_id,
                            'chrom':    rec.chrom,
                            'site':     rec.pos,
                            'genoq':    int(rec.qual),
                            'base':     allele
                        }

                        dbc.save_record('sample_reads', read)

        # count the total number of reads
        num_reads += len(reads)

        # bulk insert all the reads for this sample
        if reads:
            dbc.save_records('sample_reads', fields, reads)

        # update the interval to show we're done
        interval['finished'] = 1
        dbc.save_record('intervals', interval)

        print "INFO: Found {:,} reads for interval chr{}:{}-{}".format(num_reads, chrom, start, end)

    except Exception:
        # Put all exception text into an exception and raise that
        raise Exception("".join(traceback.format_exception(*sys.exc_info())))
