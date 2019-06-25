#!/usr/bin/env python
# -*- coding: utf-8 -*-

# standard modules
import glob
import itertools
import os
import random
from collections import defaultdict

# third party modules
import luigi
import pysam

# local modules
from alleletraj import utils
from alleletraj.ancient.vcf import BiallelicSNPsAncientVCF
from alleletraj.bam import SampleBAM, DepthOfCoveragePipeline
from alleletraj.db.conn import Database
from alleletraj.ensembl.link import EnsemblLinkPipeline
from alleletraj.modern.snps import ModernSNPsPipeline
from alleletraj.modern.vcf import ReferencePloidy, MIN_GENO_QUAL
from alleletraj.qtl.load import MIN_DAF, MergeAllLoci, PopulateQTLSNPs
from alleletraj.ref import ReferenceFASTA

# number of bases to hard clip
HARD_CLIP_DIST = 5

# minimum mapping quality (hard filtered)
HARD_MAPQ_CUTOFF = 30

# minimum base quality (hard filtered)
HARD_BASEQ_CUTOFF = 30

# diagnose mispolarized SNPs by sampling N snps and thresholding the DAF
MISPOLAR_SNPS = 10
MISPOLAR_DAF = 0.8


def list_snps_in_locus(species, chrom, start, end):
    """
    Get all the modern SNPs in this locus that fall within a QTL.
    """
    dbc = Database(species)

    snps = dbc.get_records_sql("""
        SELECT DISTINCT ms.site, ms.ancestral, ms.derived
          FROM modern_snps ms
          JOIN qtl_snps qs 
            ON qs.modsnp_id = ms.id
         WHERE ms.chrom = '{chrom}'
           AND ms.site BETWEEN {start} AND {end}
           """.format(chrom=chrom, start=start, end=end), key='site')

    return snps


class LoadAncientDiploidSNPs(utils.MySQLTask):
    """
    Load all the ancient data for SNPs with sufficient coverage to make a high quality diploid call.

    :type species: str
    :type chrom: str
    """
    species = luigi.Parameter()
    chrom = luigi.Parameter()

    def requires(self):
        yield MergeAllLoci(self.species, self.chrom)
        yield PopulateQTLSNPs(self.species, self.chrom)

        for pop, sample in self.list_samples(ancient=True):
            yield BiallelicSNPsAncientVCF(self.species, pop, sample)

    def run(self):
        # unpack the params
        bed_file, _ = self.input()[:2]
        vcf_files = [vcf_file for vcf_file, _ in self.input()[2:]]

        with bed_file.open('r') as qtl_loci, self.output().open('w') as log:

            # iterate over the loci in the BED file
            for locus in qtl_loci:
                chrom, start, end = locus.split()

                snps = list_snps_in_locus(self.species, chrom, start, end)

                log.write("INFO: Scanning locus chr{}:{}-{} for {:,} SNPs".format(chrom, start, end, len(snps)) + '\n')

                # the column headers for batch inserting into the db
                fields = ('sample_id', 'chrom', 'site', 'base', 'genoq')

                num_reads = 0

                samples = self.list_samples(ancient=True)

                # check every sample for reads in this locus
                for (pop, sample), vcf_file in itertools.izip(samples, vcf_files):

                    log.write("INFO: Scanning locus chr{}:{}-{} in sample {}".format(chrom, start, end, sample) + '\n')

                    sample_id = samples[(pop, sample)]['id']

                    reads = []

                    # parse the vcf with pysam
                    for rec in pysam.VariantFile(vcf_file).fetch(chrom, int(start), int(end)):

                        # NOTE VariantRecord.pos is 1 based
                        # https://pysam.readthedocs.io/en/latest/api.html#pysam.VariantRecord.pos
                        site = rec.pos

                        # skip all sites not ascertained in the modern samples
                        if site not in snps:
                            continue

                        # get the genotype quality
                        genoq = int(rec.qual)

                        for base in rec.samples[sample].alleles:
                            # setup the record to insert, in this order
                            read = (sample_id, chrom, site, base, genoq)

                            # store the read so we can batch insert later
                            reads.append(read)

                    # count the number of reads
                    num_reads += len(reads)

                    # bulk insert all the reads for this sample
                    if reads:
                        self.dbc.save_records('sample_reads', fields, reads)

                log.write("INFO: Found {:,} reads for locus chr{}:{}-{}".format(num_reads, chrom, start, end) + '\n')


class LoadAncientSNPs(utils.MySQLTask):
    """
    Load all the ancient data for SNPs that fall within the loci of interest.

    :type species: str
    :type chrom: str
    """
    species = luigi.Parameter()
    chrom = luigi.Parameter()

    def requires(self):
        yield ReferenceFASTA(self.species)
        yield ReferencePloidy(self.species)
        yield MergeAllLoci(self.species, self.chrom)
        yield ModernSNPsPipeline(self.species)
        yield EnsemblLinkPipeline(self.species)

        for pop, sample in self.list_samples(ancient=True):
            yield SampleBAM(self.species, pop, sample)

    def run(self):
        # unpack the params
        (ref_file, _), pld_file, bed_file = self.input()[0:3]
        bam_files = [bam_file for bam_file, _ in self.input()[5:]]

        with bed_file.open('r') as qtl_loci, self.output().open('w') as log:

            # iterate over the loci in the BED file
            for locus in qtl_loci:
                chrom, start, end = locus.split()

                snps = list_snps_in_locus(self.species, chrom, start, end)

                log.write("INFO: Scanning locus chr{}:{}-{} for {:,} SNPs".format(chrom, start, end, len(snps)) + '\n')

                # the column headers for batch inserting into the db
                fields = ('sample_id', 'chrom', 'site', 'base', 'mapq', 'baseq', 'dist')

                num_reads = 0

                samples = self.list_samples(ancient=True)

                # check every sample for reads in this locus
                for (pop, sample), bam_file in itertools.izip(samples, bam_files):

                    log.write("INFO: Scanning locus chr{}:{}-{} in sample {}".format(chrom, start, end, sample) + '\n')

                    sample_id = samples[(pop, sample)]['id']

                    # buffer the reads so we can bulk insert them into the db
                    reads = []

                    # open the BAM file for reading
                    with pysam.AlignmentFile(bam_file.path, 'rb') as pysam_align:

                        for pileup_column in pysam_align.pileup(chrom, int(start), int(end)):

                            # NOTE PileupColumn.reference_pos is 0 based
                            # see http://pysam.readthedocs.io/en/latest/api.html#pysam.PileupColumn.reference_pos
                            site = pileup_column.reference_pos + 1

                            # skip all sites not ascertained in the modern samples
                            if site not in snps:
                                continue

                            quality = []

                            # iterate over all the reads for this site
                            for pileup_read in pileup_column.pileups:

                                # skip alignments that don't have a base at this site (i.e. indels)
                                if pileup_read.is_del or pileup_read.is_refskip:
                                    continue

                                # get the read position
                                read_pos = pileup_read.query_position

                                # get the aligned base for this read
                                base = pileup_read.alignment.query_sequence[read_pos]

                                # get the map quality
                                mapq = pileup_read.alignment.mapping_quality

                                # get the base quality
                                baseq = pileup_read.alignment.query_qualities[read_pos]

                                # get the overall length of the read
                                read_length = len(pileup_read.alignment.query_sequence)

                                # how close is the base to the edge of the read
                                dist = min(read_pos, read_length - read_pos)

                                # store the read so we can batch insert later
                                read = (sample_id, chrom, site, base, mapq, baseq, dist)

                                # apply the hard filters
                                if mapq >= HARD_MAPQ_CUTOFF and baseq >= HARD_BASEQ_CUTOFF and dist >= HARD_CLIP_DIST:
                                    quality.append(read)

                            # call a pseudo-haploid genotype at random
                            if quality:
                                reads.append(random.choice(quality))

                    # count the total number of reads
                    num_reads += len(reads)

                    # bulk insert all the reads for this sample
                    if reads:
                        self.dbc.save_records('sample_reads', fields, reads)

                log.write("INFO: Found {:,} reads for locus chr{}:{}-{}".format(num_reads, chrom, start, end) + '\n')


class FlagMispolarizedSNPs(utils.MySQLTask):
    """
    Flag SNPs which may be mispolarized.

    Get the first MISPOLAR_SNPS snps and calculate the derived allele frequency. If daf > MISPOLAR_DAF then flag the SNP
    as mispolarized.

    :type species: str
    :type chrom: str
    """
    species = luigi.Parameter()
    chrom = luigi.Parameter()

    def requires(self):
        return LoadAncientSNPs(self.species, self.chrom)

    # noinspection SqlAggregates
    def queries(self):

        mispolar = self.dbc.get_records_sql("""
            SELECT ms.id,
                   LENGTH(REPLACE(SUBSTRING(
                        GROUP_CONCAT(
                            IF(sr.base = ms.derived, 1, 0)
                            ORDER BY s.bp_median DESC
                            SEPARATOR ''
                        ), 1, {snps}), 0, '')) / {snps} AS daf
              FROM modern_snps ms
              JOIN sample_reads sr
                ON sr.chrom = ms.chrom
               AND sr.site = ms.site
              JOIN samples s
                ON s.id = sr.sample_id
            WHERE ms.chrom = '{chrom}'
          GROUP BY ms.id
            HAVING daf >= {daf}
               """.format(chrom=self.chrom, snps=MISPOLAR_SNPS, daf=MISPOLAR_DAF))

        for modsnp_id in mispolar:
            self.dbc.save_record('modern_snps', {'id': modsnp_id, 'mispolar': 1})


class AncientSNPsPipeline(utils.PipelineWrapperTask):
    """
    Populate the modern_snps table.

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):
        # calculate the depth of coverage for all ancient samples
        yield DepthOfCoveragePipeline(self.species, ancient=True)

        # process SNPs for all chromosomes
        for chrom in self.chromosomes:
            yield FlagMispolarizedSNPs(self.species, chrom)


if __name__ == '__main__':
    luigi.run()
