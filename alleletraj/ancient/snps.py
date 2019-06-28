#!/usr/bin/env python
# -*- coding: utf-8 -*-

# standard modules
import itertools
import os
import random
from collections import defaultdict

# third party modules
import luigi
import pysam

# local modules
from alleletraj import utils
from alleletraj.ancient.vcf import BCFToolsFilterAncientVCF
from alleletraj.bam import DepthOfCoveragePipeline, SampleBAM
from alleletraj.db.conn import Database
from alleletraj.qtl.load import MergeAllLoci, PopulateQTLSNPs

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
          FROM qtls q
          JOIN qtl_snps qs 
            ON qs.qtl_id = q.id
          JOIN modern_snps ms
            ON ms.id = qs.modsnp_id
         WHERE q.chrom = '{chrom}'
           AND q.start BETWEEN {start} AND {end}
           AND q.valid = 1
           """.format(chrom=chrom, start=start, end=end), key='site')

    return snps


def list_genotypes_in_locus(species, chrom, start, end):
    """
    Get all the sites in this locus that have a diploid call.
    """
    dbc = Database(species)

    records = dbc.get_records_sql("""
       SELECT DISTINCT sr.sample_id, ms.site
          FROM qtls q
          JOIN qtl_snps qs 
            ON qs.qtl_id = q.id
          JOIN modern_snps ms
            ON ms.id = qs.modsnp_id
          JOIN sample_reads sr
            ON sr.chrom = ms.chrom
           AND sr.site = ms.site
           AND sr.genoq IS NOT NULL  # only diploid calls have genotype qualities
         WHERE q.chrom = '{chrom}'
           AND q.start BETWEEN {start} AND {end}
           AND q.valid = 1
           """.format(chrom=chrom, start=start, end=end), key=None)

    diploids = defaultdict(set)

    # nest by sample
    for rec in records:
        diploids[rec['sample_id']].add(rec['site'])

    return diploids


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
            yield BCFToolsFilterAncientVCF(self.species, pop, sample)

    def run(self):
        # unpack the params
        bed_file, _ = self.input()[:2]
        vcf_files = [vcf_file for vcf_file, _ in self.input()[2:]]

        with bed_file.open('r') as qtl_loci, self.output().open('w') as log:

            samples = self.list_samples(ancient=True)

            # iterate over the loci in the BED file
            for locus in qtl_loci:
                chrom, start, end = locus.split()

                snps = list_snps_in_locus(self.species, chrom, start, end)
                qtl_sites = set(snps.keys())  # convert to set for much faster lookups

                log.write("INFO: Scanning locus chr{}:{}-{} for {:,} SNPs".format(chrom, start, end, len(snps)) + '\n')

                # the column headers for batch inserting into the db
                fields = ('sample_id', 'chrom', 'site', 'base', 'genoq')

                num_reads = 0

                # check every sample for reads in this locus
                for (pop, sample), vcf_file in itertools.izip(samples, vcf_files):

                    log.write("INFO: Scanning locus chr{}:{}-{} in sample {}".format(chrom, start, end, sample) + '\n')

                    sample_id = samples[(pop, sample)]['id']

                    reads = []

                    try:
                        # parse the vcf with pysam
                        records = pysam.VariantFile(vcf_file.path).fetch(chrom, int(start), int(end))

                    except ValueError as error:
                        if os.path.isfile(vcf_file.path + '.tbi'):
                            continue  # VCF is empty
                        else:
                            raise error

                    for rec in records:

                        # NOTE VariantRecord.pos is 1 based
                        # https://pysam.readthedocs.io/en/latest/api.html#pysam.VariantRecord.pos
                        site = rec.pos

                        # skip all sites not ascertained in the modern samples
                        if site not in qtl_sites:
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


class LoadAncientHaploidSNPs(utils.MySQLTask):
    """
    Load all the ancient data for SNPs where we cannot make a diploid call.

    :type species: str
    :type chrom: str
    """
    species = luigi.Parameter()
    chrom = luigi.Parameter()

    def requires(self):
        yield MergeAllLoci(self.species, self.chrom)
        yield LoadAncientDiploidSNPs(self.species, self.chrom)

        for pop, sample in self.list_samples(ancient=True):
            yield SampleBAM(self.species, pop, sample)

    def run(self):
        # unpack the params
        bed_file, _ = self.input()[:2]
        bam_files = [bam_file for bam_file, _ in self.input()[2:]]

        with bed_file.open('r') as qtl_loci, self.output().open('w') as log:

            samples = self.list_samples(ancient=True)

            # iterate over the loci in the BED file
            for locus in qtl_loci:
                chrom, start, end = locus.split()

                snps = list_snps_in_locus(self.species, chrom, start, end)
                qtl_sites = set(snps.keys())  # convert to set for much faster lookups

                # get the list of diploid calls for all samples in this locus
                diploids = list_genotypes_in_locus(self.species, chrom, start, end)

                log.write("INFO: Scanning locus chr{}:{}-{} for {:,} SNPs".format(chrom, start, end, len(snps)) + '\n')

                # the column headers for batch inserting into the db
                fields = ('sample_id', 'chrom', 'site', 'base', 'mapq', 'baseq', 'dist')

                num_reads = 0

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
                            if site not in qtl_sites:
                                continue

                            # also skip sites where we have a diploid call for this sample
                            if site in diploids.get(sample_id, {}):
                                continue

                            alleles = [snps[site]['ancestral'], snps[site]['derived']]

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
                                if mapq >= HARD_MAPQ_CUTOFF and baseq >= HARD_BASEQ_CUTOFF \
                                        and dist >= HARD_CLIP_DIST and base in alleles:
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


class ValidateSampleReads(utils.MySQLTask):
    """
    Validate that all the ancient data was loaded correctly.

    Specifically, check that there are exactly 1 haploid read, or 2 diploid reads, per sample per site in the db.

    :type species: str
    :type chrom: str
    """
    species = luigi.Parameter()
    chrom = luigi.Parameter()

    def requires(self):
        return LoadAncientHaploidSNPs(self.species, self.chrom)

    def queries(self):

        if self.chrom in self.autosomes:
            malformed = self.dbc.get_records_sql("""
                SELECT chrom, site, sample_id, COUNT(baseq) haploid, COUNT(genoq) diploid
                  FROM sample_reads
                 WHERE chrom = '{chrom}'
              GROUP BY sample_id, site
            HAVING NOT ((haploid = 1 AND diploid = 0)
                    OR  (haploid = 0 AND diploid = 2))
                   """.format(chrom=self.chrom), key=None)
        else:
            # handle sex chromosomes
            sex1, sex2 = ('M', 'F') if self.chrom == 'X' else ('F', 'M')

            malformed = self.dbc.get_records_sql("""
                SELECT sr.chrom, sr.site, sr.sample_id, s.sex, COUNT(sr.baseq) haploid, COUNT(sr.genoq) diploid
                  FROM samples s
                  JOIN sample_reads sr
                    ON sr.sample_id = s.id
                 WHERE chrom = '{chrom}'
              GROUP BY sample_id, site
            HAVING NOT ((haploid = 1 AND diploid = 0)
                    OR  (haploid = 0 AND diploid = 1 AND sex = '{sex1}')
                    OR  (haploid = 0 AND diploid = 2 AND sex = '{sex2}'))
                   """.format(chrom=self.chrom, sex1=sex1, sex2=sex2), key=None)


        if malformed:
            raise Exception("ERROR: Bad data in sample_reads table - {}".format(malformed))


class FlagMispolarizedSNPs(utils.MySQLTask):
    """
    Flag SNPs which may be mispolarized.

    Calculate the DAF for the oldest and youngest 10 observations, and flag as possibly mispolarised if the DAF in the
    oldest bin is > 80% and the DAF in the youngest bin is the same or lower.

    :type species: str
    :type chrom: str
    """
    species = luigi.Parameter()
    chrom = luigi.Parameter()

    db_lock_tables = ['modern_snps_{chrom}']

    def requires(self):
        return ValidateSampleReads(self.species, self.chrom)

    def queries(self):

        # noinspection SqlAggregates, SqlWithoutWhere
        self.dbc.execute_sql("""
            UPDATE modern_snps ms
              JOIN (

            SELECT ms.id,
                   LENGTH(REPLACE(SUBSTRING(
                        GROUP_CONCAT(
                            IF(sr.base = ms.derived, 1, 0)
                            ORDER BY s.bp_median DESC
                            SEPARATOR ''
                        ), 1, {snps}), 0, '')) / {snps} AS daf_oldest,

                   LENGTH(REPLACE(SUBSTRING(
                        GROUP_CONCAT(
                            IF(sr.base = ms.derived, 1, 0)
                            ORDER BY s.bp_median ASC
                            SEPARATOR ''
                        ), 1, {snps}), 0, '')) / {snps} AS daf_youngest
              FROM modern_snps ms
              JOIN sample_reads sr
                ON sr.chrom = ms.chrom
               AND sr.site = ms.site
              JOIN samples s
                ON s.id = sr.sample_id
            WHERE ms.chrom = '{chrom}'
          GROUP BY ms.id
            HAVING daf_oldest >= {daf}
               AND daf_oldest >= daf_youngest

              ) AS mis 
                ON mis.id = ms.id
               SET ms.mispolar = 1
             WHERE ms.chrom = '{chrom}'
               """.format(chrom=self.chrom, snps=MISPOLAR_SNPS, daf=MISPOLAR_DAF))


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
