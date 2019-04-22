#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function

import luigi
import pysam
import os
from random import shuffle

from collections import defaultdict
from natsort import natsorted

from pipeline_consts import MIN_GENO_QUAL, MIN_DAF

from pipeline_qtls import PopulateAllLoci
from pipeline_modern_snps import ModernSNPsPipeline
from pipeline_ensembl import EnsemblPipeline
from pipeline_snp_call import ExternalFASTA
from pipeline_samples import LoadSamples
from pipeline_utils import PipelineTask, merge_intervals, run_cmd

# minimum depth of coverage to call diploid genotypes
MIN_GENO_DEPTH = 10

# minimum mapping quality (hard filtered)
HARD_MAPQ_CUTOFF = 20

# minimum base quality (hard filtered)
HARD_BASEQ_CUTOFF = 20


class PopulateIntervals(PipelineTask):
    """
    Merge all the overlapping QTL and pseudo-QTL windows, to determine the unique list of intervals for which we need to
    load sample reads.

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):
        yield PopulateAllLoci(self.species)

    def output(self):
        return luigi.LocalTarget('db/{}-{}.log'.format(self.basename, self.classname))

    def run(self):
        # open a db connection
        dbc = self.db_conn()

        # get all the unique QTL windows
        results = dbc.get_records_sql("""
            SELECT DISTINCT q.chrom, q.start, q.end
              FROM qtls q
             WHERE q.valid = 1
          ORDER BY q.chrom, q.start, q.end""", key=None)

        intervals = defaultdict(list)

        for result in results:
            # group the intervals by chrom
            intervals[result['chrom']].append((result['start'], result['end']))

        num_sites = 0
        added = 0
        removed = 0

        for chrom in natsorted(intervals.keys()):
            # merge overlapping intervals (up to a maximum size)
            intervals[chrom] = list(merge_intervals(intervals[chrom]))

            for start, end in intervals[chrom]:

                # compose the interval record
                record = {'chrom': chrom, 'start': start, 'end': end}

                # check if this interval already exists
                if dbc.get_record('intervals', record):
                    continue

                # get any overlapping intervals
                overlap = dbc.get_records_sql("""
                    SELECT *
                      FROM intervals
                     WHERE chrom = '{chrom}'
                       AND end > {start}
                       AND start < {end}
                       """.format(chrom=chrom, start=start, end=end))

                for interval_id in overlap:
                    # delete the old intervals
                    dbc.delete_records('sample_reads', {'interval_id': interval_id})
                    dbc.delete_records('intervals_snps', {'interval_id': interval_id})
                    dbc.delete_records('intervals', {'id': interval_id})

                    removed += len(overlap)

                # keep track of what we've done
                added += 1
                num_sites += end - start

                # save the new interval
                dbc.save_record('intervals', record)

        with self.output().open('w') as fout:
            fout.write("INFO: Added {:,} intervals ({:,} bp), removed {:,} intervals".format(added, num_sites, removed))


class PopulateIntervalSNPs(PipelineTask):
    """
    Now we have ascertained all the modern SNPs, let's find those that intersect with the unique intervals.

    :type species: str
    :type population: str
    :type interval: str
    """
    species = luigi.Parameter()
    population = luigi.Parameter()
    interval = luigi.IntParameter()

    def requires(self):
        yield PopulateIntervals(self.species)
        yield ModernSNPsPipeline(self.species)
        yield EnsemblPipeline(self.species)

    def output(self):
        return luigi.LocalTarget('db/{}-{}.log'.format(self.basename, self.classname))

    def run(self):
        dbc = self.db_conn()

        # insert linking records to make future queries much quicker
        exec_time = dbc.execute_sql("""
            INSERT INTO intervals_snps (interval_id, modsnp_id)
                 SELECT i.id, ms.id
                   FROM intervals i
                   JOIN modern_snps ms 
                     ON ms.population = '{population}'
                    AND ms.chrom = i.chrom
                    AND ms.site BETWEEN i.start AND i.end
                  WHERE i.id = {id}
                    AND ms.daf >= {daf}""".format(population=self.population, id=self.interval, daf=MIN_DAF))

        with self.output().open('w') as fout:
            fout.write('INFO: Execution took {}'.format(exec_time))


class ProcessInterval(PipelineTask):
    """
    Scan all the intervals across this chromosome and add the covered bases to the DB, as long as they are variable in
    the modern data.

    # TODO in table `sample_reads` replace (interval_id ?, chrom, site) with modsnp_id
    # TODO increase hard threshold to >= 30 / or / delete where called = 0

    :type species: str
    """
    species = luigi.Parameter()
    population = luigi.Parameter()
    interval = luigi.IntParameter()

    def requires(self):
        yield ExternalFASTA(self.species)
        yield PopulateIntervalSNPs(self.species, self.population, self.interval)
        yield LoadSamples(self.species)

    def output(self):
        return luigi.LocalTarget('db/{}-{}.log'.format(self.basename, self.classname))

    def run(self):
        # open a db connection
        dbc = self.db_conn()

        # get the reference genome
        ref_file = self.input()[0]

        interval = dbc.get_record('intervals', {'id': self.interval})

        # unpack the interval
        interval_id, chrom, start, end = interval['id'], interval['chrom'], interval['start'], interval['end']

        # get all the valid samples
        samples = dbc.get_records_sql("""
            SELECT s.*, GROUP_CONCAT(sf.path) paths
             FROM samples s
             JOIN sample_files sf
               ON sf.sample_id = s.id
            WHERE s.valid = 1
         GROUP BY s.id""")

        # get all the modern SNPs in this interval
        snps = dbc.get_records_sql("""
            SELECT ms.site, ev.ref, ev.alt, ms.ancestral, ms.derived
              FROM intervals i
              JOIN intervals_snps s
                ON s.interval_id = i.id
              JOIN modern_snps ms
                ON ms.id = s.modsnp_id
         LEFT JOIN ensembl_variants ev
                ON ev.id = ms.variant_id
             WHERE i.id = {id}""".format(id=interval_id), key='site')

        with self.output().open('w') as log:
            log.write("INFO: Scanning interval chr{}:{}-{} for {:,} SNPs".format(chrom, start, end, len(snps)))

            # handle chr1 vs. 1 chromosome names
            contig = 'chr' + chrom if self.species == 'horse' else chrom

            # not all SNPs have a dbsnp entry, so we need to scan the reference to find which alleles are REF/ALT
            # because bcftools needs this info to constrain the diploid genotype calls
            with pysam.FastaFile(ref_file.path) as fasta_file:

                for site in snps.keys():
                    if not snps[site]['ref']:
                        # fetch the ref and snp alleles
                        ref_allele = fasta_file.fetch(contig, site-1, site)
                        snp_alleles = [snps[site]['ancestral'], snps[site]['derived']]

                        if ref_allele not in snp_alleles:
                            log.write("WARNING: chr{}:{} REF allele {} not found in SNP alleles {}"
                                      .format(chrom, site, ref_allele, snp_alleles))

                            # remove this site
                            snps.pop(site)

                            continue

                        snp_alleles.remove(ref_allele)

                        # set the REF/ALT alleles
                        snps[site]['ref'] = ref_allele
                        snps[site]['alt'] = snp_alleles.pop()

            # the column headers for batch inserting into the db
            fields = ('interval_id', 'sample_id', 'chrom', 'site', 'base', 'mapq', 'baseq', 'dist')

            num_reads = 0

            # randomise the order of samples to reduce disk I/O for parallel jobs
            ids = samples.keys()
            shuffle(ids)

            # check all the samples for coverage in this interval
            for sample_id in ids:

                # get the sample
                sample = samples[sample_id]

                log.write("INFO: Scanning interval chr{}:{}-{} in sample {}"
                          .format(chrom, start, end, sample['accession']))

                # buffer the reads so we can bulk insert them into the db
                reads = defaultdict(list)

                # there may be multiple BAM files for each sample
                for path in sample['paths'].split(','):

                    # open the BAM file for reading
                    with pysam.AlignmentFile(path, 'rb') as bamfile:

                        # get the full interval
                        for pileupcolumn in bamfile.pileup(contig, start, end + 1):

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
                    log.write("INFO: Calling diploid bases in {:,} sites for sample {}"
                              .format(len(diploid), sample_id))

                    # TODO make these temp files
                    pos_file = 'vcf/diploid-int{}-sample{}.tsv'.format(interval_id, sample_id)
                    vcf_file = 'vcf/diploid-int{}-sample{}.vcf'.format(interval_id, sample_id)

                    # sort the diploid positions
                    diploid.sort()

                    # save all the callable positions to a file
                    with open(pos_file, 'w') as fout:
                        fout.write("\n".join("{}\t{}\t{},{}".format(contig, site, snps[site]['ref'], snps[site]['alt'])
                                             for (chrom, site) in diploid))

                    targets = "{}.gz".format(pos_file)

                    # bgzip and index the target file
                    run_cmd(["bgzip -c {} > {}".format(pos_file, targets)], shell=True)
                    run_cmd(["tabix -s1 -b2 -e2 {}".format(targets)], shell=True)

                    # restrict the callable region using the interval start and end
                    region = "{}:{}-{}".format(contig, start, end)

                    # use all the BAM files
                    bam_files = " ".join(sample['paths'].split(','))

                    params = {'region': region, 'targets': targets, 'ref': ref_file.path, 'bams': bam_files,
                              'vcf': vcf_file}

                    # call bases with bcftools (and drop indels and other junk)
                    # uses both --region (random access) and --targets (streaming) for optimal speed
                    # see https://samtools.github.io/bcftools/bcftools.html#mpileup
                    # TODO none of --samples-file, --ploidy or --ploidy-file given, assuming all sites are diploid
                    cmd = "bcftools mpileup --region {region} --targets-file {targets} --fasta-ref {ref} {bams} " \
                          "| bcftools call --multiallelic-caller --targets-file {targets} --constrain alleles -O v " \
                          "| bcftools view --exclude-types indels,bnd,other --exclude INFO/INDEL=1 --output-file {vcf}"\
                          .format(**params)

                    # run the base calling
                    run_cmd([cmd], shell=True)

                    # parse the results with pysam
                    for rec in pysam.VariantFile(vcf_file).fetch():

                        # get the genotype call for this site (without having to know the @SM code used in the BAM file)
                        geno = rec.samples.values().pop()['GT']

                        # get the genotype quality
                        genoq = int(rec.qual)

                        if genoq >= MIN_GENO_QUAL:
                            # the genotype is good, so drop the raw reads
                            reads.pop((chrom, rec.pos))

                        # decode the GT notation into allele calls (e.g. 0/0, 0/1, 1/1)
                        alleles = [rec.alleles[idx] for idx in geno if idx is not None]

                        for allele in alleles:
                            # compose the read records
                            read = {
                                'interval_id': interval_id,
                                'sample_id':   sample_id,
                                'chrom':       chrom,
                                'site':        rec.pos,
                                'genoq':       genoq,
                                'base':        allele
                            }

                            dbc.save_record('sample_reads', read)

                    # delete the used VCF file
                    os.remove(vcf_file)

                # apply hard filters before inserting (otherwise we swamp the DB with too many low quality reads)
                reads = [read for (chrom, site) in reads for read in reads[(chrom, site)]
                         if read[fields.index('mapq')] >= HARD_MAPQ_CUTOFF
                         and read[fields.index('baseq')] >= HARD_BASEQ_CUTOFF]

                # count the total number of reads
                num_reads += len(reads)

                # bulk insert all the reads for this sample
                if reads:
                    dbc.save_records('sample_reads', fields, reads)

            log.write("INFO: Found {:,} reads for interval chr{}:{}-{}".format(num_reads, chrom, start, end))


class SampleReadsPipeline(PipelineTask):
    """
    Load all the samples reads.

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):
        yield PopulateIntervals(self.species)

    def output(self):
        return luigi.LocalTarget('db/{}-{}.log'.format(self.basename, self.classname))

    def run(self):
        dbc = self.db_conn()

        # get all the intervals
        intervals = dbc.get_records('intervals')

        for population in self.populations:
            for interval in intervals:
                # process all the intervals
                yield ProcessInterval(self.species, population, interval)

        with self.output().open('w') as fout:
            fout.write('INFO: Processed {} intervals'.format(len(intervals)))


if __name__ == '__main__':
    luigi.run()
