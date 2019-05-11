#!/usr/bin/env python
# -*- coding: utf-8 -*-

# standard modules
import random
import re
import struct
from collections import OrderedDict, Counter

# third party modules
import luigi
import pysam

# local modules
from alleletraj import utils
from alleletraj.bam import AlignedBAM
from alleletraj.const import CHROMOSOMES, REF_ASSEMBLY

# column delimiter for PED files
PLINK_COL_DELIM = ' '

# status codes for PED files
PLINK_UNKNOWN = 0
PLINK_MISSING_PHENO = -9
PLINK_MISSING_GENO = 0

# min genotyping rate
PLINK_MIN_GENO = 90

MIN_MAP_QUAL = 30  # TODO
MIN_BASE_QUAL = 30  # TODO
SOFT_CLIP_DIST = 5  # TODO


def call_ancient_genotypes(species, population, sample, bam_path, bed_prefix, sex=PLINK_UNKNOWN, min_mapq=MIN_MAP_QUAL,
                           min_baseq=MIN_BASE_QUAL, soft_clip=SOFT_CLIP_DIST, rand_call=True):
    """
    Generates a plink BED file (binary PED) with genotypes drawn from a given BAM file.

    The BAM file is converted to a pileup, quality filtered, then a pseudo-haploid genotype is called by either
    choosing a random base or by calling the major allele.

    N.B. PLINK does not support more than 2^31 - 3 variants :-(

    :param species:    Species name
    :param population: Population code
    :param sample:     Sample code
    :param bam_path:   Path to the BAM file
    :param bed_prefix: Prefix for the BED file to output
    :param sex:        Optional, sex of the sample (1 = male, 2 = female, 0 = unknown)
    :param min_mapq:   Optional, min MapQ to accept
    :param min_baseq:  Optional, min BaseQ to accept
    :param soft_clip:  Optional, prefer bases away from read ends
    :param rand_call:  Optional, call a random base, else call major allele
    """

    # make a FAM file for this sample (see https://www.cog-genomics.org/plink2/formats#fam)
    fam = OrderedDict([
        ('FID',       population),
        ('IID',       sample),
        ('father',    PLINK_UNKNOWN),
        ('mother',    PLINK_UNKNOWN),
        ('sex',       sex),
        ('phenotype', PLINK_MISSING_PHENO)

    ])

    # save the FAM file
    with open(bed_prefix + '.fam', 'w') as fam_fout:
        fam_fout.write(PLINK_COL_DELIM.join(str(val) for val in fam.values()) + '\n')

    # open the BED and BIM files for writing
    with open(bed_prefix + '.bed', 'wb') as bed_fout, open(bed_prefix + '.bim', 'w') as bim_fout:

        # initialise the BED file (see https://www.cog-genomics.org/plink2/formats#bed)
        bed_fout.write(struct.pack('b', 0x6c))
        bed_fout.write(struct.pack('b', 0x1b))
        bed_fout.write(struct.pack('b', 0x01))

        # open the BAM file for reading
        with pysam.AlignmentFile(bam_path, 'rb') as bam_file:

            # iterate over the whole genome, one chrom at a time
            for contig in CHROMOSOMES[REF_ASSEMBLY[species]]:

                for pileup_column in bam_file.pileup(contig):

                    # get the chrom and pos of the current site in the BAM file
                    chrom, pos = pileup_column.reference_name, \
                                 pileup_column.reference_pos + 1  # NB. reference_pos is 0 based

                    quality_bases = []
                    clipped_bases = []

                    # iterate over all the reads for this site
                    for pileup_read in pileup_column.pileups:

                        # skip alignments that don't have a base at this site (i.e. indels)
                        if pileup_read.is_del or pileup_read.is_refskip:
                            continue

                        # skip bad alignments
                        if pileup_read.alignment.mapping_quality < min_mapq:
                            continue

                        # skip bad base calls
                        if pileup_read.alignment.query_qualities[pileup_read.query_position] < min_baseq:
                            continue

                        # get the aligned base for this read
                        base = pileup_read.alignment.query_sequence[pileup_read.query_position]

                        # get the overall length of the read
                        read_length = len(pileup_read.alignment.query_sequence)

                        # get distance to the nearest end of the read
                        read_pos = min(pileup_read.query_position, read_length-pileup_read.query_position)

                        # soft clip bases near the edges of the read
                        if read_pos <= soft_clip:
                            clipped_bases.append(base)
                        else:
                            quality_bases.append(base)

                    # do we have any
                    if len(quality_bases) > 0 or len(clipped_bases) > 0:

                        if rand_call:
                            # choose a random allele (preferring those outside the soft clip distance)
                            allele = random.choice(quality_bases) if len(quality_bases) > 0 \
                                else random.choice(clipped_bases)
                        else:
                            # call most common allele
                            allele = Counter(quality_bases).most_common(1)[0][0] if len(quality_bases) > 0 \
                                else Counter(clipped_bases).most_common(1)[0][0]

                        # strip the 'chr' prefix
                        chrom = re.sub('^chr', '', contig)

                        # make a locus entry for the BIM file (see https://www.cog-genomics.org/plink2/formats#bim)
                        bim = OrderedDict([
                            ('chromosome',   chrom),
                            ('identifier',   '{}-{}'.format(chrom, pos)),  # make up a unique SNP name
                            ('centimorgans', 0),                           # safe to use dummy value of '0'
                            ('coordinate',   pos),
                            ('allele1',      PLINK_UNKNOWN),
                            ('allele2',      allele)
                        ])

                        # write the locus to the BIM file
                        bim_fout.write(PLINK_COL_DELIM.join(str(val) for val in bim.values()) + '\n')

                        # write the genotype to the BED file (see https://www.cog-genomics.org/plink2/formats#bed)
                        bed_fout.write(struct.pack('b', 3))  # 3 = homozygous for second allele in .bim file


class CallAncientGenotypes(utils.PipelineTask):
    """
    Call ancient genotypes by randomly selecting a base.

    TODO ANGSD doHaploCall sampling a random base
    TODO The haploid data was pseudo-diploidized and converted to plink format using ANGSD

    :type species: str
    :type population: str
    :type sample: str
    """
    species = luigi.Parameter()
    population = luigi.Parameter()
    sample = luigi.Parameter()

    def requires(self):
        return AlignedBAM(self.species, self.population, self.sample)

    def output(self):
        return [luigi.LocalTarget('data/plink/{}.{}'.format(self.basename, ext)) for ext in
                ['bed', 'bim', 'fam', 'log']]

    def run(self):
        # unpack the inputs/outputs
        bam_file, _ = self.input()
        bed_file, _, _, log_file = self.output()

        # TODO add the sex calls
        # random call the ancient genotypes
        call_ancient_genotypes(self.species, self.population, self.sample, bam_file.path, utils.trim_ext(bed_file.path))

        # write a log file to show this task finished
        with log_file.open('w') as log_fout:
            log_fout.write('Done!')
