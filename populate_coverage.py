#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function

from discover_snps import *
from populate_qtls import *

import pysam as ps
import multiprocessing as mp
import traceback
import sys
import itertools
import socket
import subprocess

from natsort import natsorted

from utilities import *

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

# enforce max interval size of 1 Gb
MAX_INTERVAL_SIZE = int(1e6)


def merge_intervals(ranges):
    """
    Merge overlapping intervals, so we only check each site once
    """
    saved = list(ranges[0])

    # sort the intervals by start, and iterate over the list
    for start, end in sorted([sorted(t) for t in ranges]):
        if start <= saved[1]:
            # if the current interval overlaps the saved one then take the largest end point
            saved[1] = max(saved[1], end)
        else:
            # enforce max interval size
            for i in range(saved[0], saved[1], MAX_INTERVAL_SIZE):
                yield (i, min(saved[1], i + MAX_INTERVAL_SIZE - 1))

            saved[0] = start
            saved[1] = end

    # return the final interval
    for i in range(saved[0], saved[1], MAX_INTERVAL_SIZE):
        yield (i, min(saved[1], i + MAX_INTERVAL_SIZE - 1))


def populate_intervals(species):

    # open a db connection
    dbc = db_conn()

    # get all the unique QTL windows
    results = dbc.get_records_sql("""
        SELECT q.chrom, q.start, q.end
          FROM qtls q
         WHERE q.species = '{species}'
           AND q.valid = 1
      GROUP BY q.chrom, q.start, q.end""".format(species=species), key=None)

    intervals = defaultdict(list)

    for result in results:
        # group the intervals by chrom
        intervals[result['chrom']].append((result['start'], result['end']))

    num_sites = 0
    add_intvals = 0
    del_intvals = 0

    for chrom in natsorted(intervals.keys()):
        # merge overlapping intervals
        intervals[chrom] = list(merge_intervals(intervals[chrom]))

        for start, end in intervals[chrom]:

            # compose the interval record
            record = {'species': species, 'chrom': chrom, 'start': start, 'end': end}

            # check if this interval already exists
            if dbc.get_record('intervals', record):
                continue

            # get any overlapping intervals
            overlap = dbc.get_records_sql("""
                SELECT *
                  FROM intervals
                 WHERE species = '{species}'
                   AND chrom = '{chrom}'
                   AND end > {start}
                   AND start < {end}
                   """.format(species=species, chrom=chrom, start=start, end=end))

            for interval_id in overlap:
                # delete the old intervals
                dbc.delete_records('sample_reads', {'interval_id': interval_id})
                dbc.delete_records('intervals_snps', {'interval_id': interval_id})
                dbc.delete_records('intervals', {'id': interval_id})

                del_intvals += len(overlap)

            # keep track of what we've done
            add_intvals += 1
            num_sites += end - start

            # save the new interval
            dbc.save_record('intervals', record)

    print("INFO: Added {:,} {} intervals ({:,} bp), deleted {:,} old intervals"
          .format(add_intvals, species, num_sites, del_intvals))


def populate_interval_snps(species):
    """
    Now we have ascertained all the modern SNPs, let's find those that intersect with the unique intervals.
    """
    dbc = db_conn()

    print("INFO: Populating all the {} interval SNPs".format(species))

    # tidy up an unfinished interval SNPs
    dbc.execute_sql("""
        DELETE s 
          FROM intervals_snps s
          JOIN intervals i
            ON i.id = s.interval_id
         WHERE i.finished = 0""")

    # process the query by chromosome to avoid buffering
    chroms = natsorted(CHROM_SIZE[species].keys())

    # insert linking records to make future queries much quicker
    for chrom in chroms:
        dbc.execute_sql("""
            INSERT INTO intervals_snps (interval_id, modsnp_id)
                 SELECT i.id, ms.id
                   FROM intervals i
                   JOIN modern_snps ms 
                     ON ms.species = i.species
                    AND ms.chrom = i.chrom
                    AND ms.site BETWEEN i.start AND i.end
                  WHERE i.finished = 0
                    AND i.species = '{species}'
                    AND i.chrom = '{chrom}'
                    AND ms.maf >= {minmaf}""".format(species=species, chrom=chrom, minmaf=MIN_MAF))


def populate_coverage(species):
    """
    Scan all the samples for coverage of the QTLs and save the results to the database.
    """

    # open a db connection
    dbc = db_conn()

    # get all the intervals we've not finished processing yet
    intervals = dbc.get_records('intervals', {'species': species, 'id': 1}, sort='end-start DESC')

    # get all the valid samples
    samples = dbc.get_records_sql(
        """SELECT s.*, GROUP_CONCAT(sf.path) paths
             FROM samples s
             JOIN sample_files sf
               ON sf.sample_id = s.id
            WHERE s.species = '{species}'
              AND s.valid = 1
              and s.id = 147
         GROUP BY s.id;""".format(species=species))

    print("INFO: Processing {:,} intervals in {:,} {} samples".format(len(intervals), len(samples), species))

    # before we start, tidy up any records from intervals that were not finished
    # dbc.execute_sql("""
    #     DELETE sample_reads
    #       FROM sample_reads
    #       JOIN intervals
    #         ON intervals.id = sample_reads.interval_id
    #      WHERE intervals.species = '{species}'
    #        AND intervals.finished = 0""".format(species=species))

    if MULTI_THREADED:
        # process the chromosomes with multi-threading to make this faster
        pool = mp.Pool(MAX_CPU_CORES)
        pool.map(process_interval, itertools.izip(intervals.values(), itertools.repeat(samples)))
    else:
        # process the chromosomes without multi-threading
        for interval in intervals.values():
            process_interval((interval, samples))

    print("FINISHED: Fully populated all the {} samples for {:,} intervals".format(species, len(intervals)))


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
            SELECT ms.site, ms.ancestral, ms.derived
              FROM intervals i
              JOIN intervals_snps s
                ON s.interval_id = i.id
              JOIN modern_snps ms
                ON ms.id = s.modsnp_id
             WHERE i.id = {id}""".format(id=interval_id), key='site')

        print("INFO: Scanning interval chr{}:{}-{} for {:,} SNPs".format(chrom, start, end, len(snps)))

        # the column headers for batch inserting into the db
        fields = ('interval_id', 'sample_id', 'chrom', 'site', 'base', 'mapq', 'baseq', 'dist')

        num_reads = 0

        # check all the samples for coverage in this interval
        for sample_id, sample in samples.iteritems():

            if VERBOSE:
                print("INFO: Scanning interval chr{}:{}-{} in sample {}".format(chrom, start, end, sample['accession']))

            # buffer the reads so we can bulk insert them into the db
            reads = defaultdict(list)

            # there may be multiple BAM files for each sample
            for path in sample['paths'].split(','):

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

                            # get the base quality
                            baseq = pileupread.alignment.query_qualities[read_pos]

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

                print("INFO: Calling diploid bases in {:,} sites for sample {}".format(len(diploid), sample_id))

                pos_file = 'vcf/diploid-int{}-sample{}.tsv'.format(interval_id, sample_id)
                vcf_file = 'vcf/diploid-int{}-sample{}.vcf'.format(interval_id, sample_id)

                # sort the diploid positions
                diploid.sort()

                # save all the callable positions to a file
                with open(pos_file, 'w') as fout:
                    # TODO ancestral may not be REF allele
                    fout.write("\n".join("{}\t{}\t{},{}".format(chrom, site, snps[site]['ancestral'], snps[site]['derived'])
                                         for (chrom, site) in diploid))

                targets = "{}.gz".format(pos_file)

                # bgzip and index the target file
                run_cmd(["bgzip -c > {}".format(targets)], shell=True)
                run_cmd(["tabix -s1 -b2 -e2 {}".format(targets)], shell=True)

                # restrict the callable region using the interval start and end
                region = "{}:{}-{}".format(chrom, start, end)

                # use all the BAM files
                bam_files = " ".join(sample['paths'].split(','))

                params = {'region': region, 'targets': targets, 'ref': REF_FILE, 'bams': bam_files, 'vcf': vcf_file}

                # call bases with bcftools (and drop indels and other junk)
                # uses both --region (random access) and --targets (streaming) for optimal speed
                # see https://samtools.github.io/bcftools/bcftools.html#mpileup
                cmd = "bcftools mpileup --region {region} --targets-file {targets} --fasta-ref {ref} {bams} " \
                      " | bcftools call --multiallelic-caller --targets-file {targets} --constrain alleles --output-type v " \
                      " | bcftools view --exclude-types indels,bnd,other --exclude INFO/INDEL=1 --output-file {vcf}" \
                      .format(**params)

                # run the base calling
                run_cmd([cmd], shell=True)

                # parse the results with pysam
                for rec in ps.VariantFile(vcf_file).fetch():

                    # get the genotype call for this site
                    geno = rec.samples[sample['accession']]['GT']

                    # get the genotype quality
                    genoq = int(rec.qual)

                    if genoq >= MIN_GENO_QUAL:
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
                            'genoq':    genoq,
                            'base':     allele
                        }

                        dbc.save_record('sample_reads', read)

            # apply hard filters before inserting (otherwise we swamp the DB with too many low quality reads)
            reads = [read for (chrom, site) in reads for read in reads[(chrom, site)]
                        if read[fields.index('mapq')] >= HARD_MAPQ_CUTOFF and
                           read[fields.index('baseq')] >= HARD_BASEQ_CUTOFF]

            # count the total number of reads
            num_reads += len(reads)

            # bulk insert all the reads for this sample
            if reads:
                dbc.save_records('sample_reads', fields, reads)

        # update the interval to show we're done
        interval['finished'] = 1
        dbc.save_record('intervals', interval)

        print("INFO: Found {:,} reads for interval chr{}:{}-{}".format(num_reads, chrom, start, end))

    except Exception:
        # Put all exception text into an exception and raise that
        raise Exception("".join(traceback.format_exception(*sys.exc_info())))
