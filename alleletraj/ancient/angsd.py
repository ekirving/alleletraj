#!/usr/bin/env python
# -*- coding: utf-8 -*-

# third party modules
import luigi

# local modules
from alleletraj import utils
from alleletraj.bam import SampleBAM
from alleletraj.const import CPU_CORES_LOW
from alleletraj.modern.vcf import WholeGenomeSNPsVCF

# number of bases to hard clip
HARD_CLIP_DIST = 5

# minimum mapping quality (hard filtered)
HARD_MAPQ_CUTOFF = 30

# minimum base quality (hard filtered)
HARD_BASEQ_CUTOFF = 30


class CallAncientGenotypes(utils.DatabaseTask):
    """
    Call ancient genotypes with angsd by randomly selecting a base.

    :type species: str
    :type chrom: str
    """
    species = luigi.Parameter()
    chrom = luigi.Parameter()

    resources = {'cpu-cores': 2, 'ram-gb': 4}

    def requires(self):
        for pop, sample in self.list_samples(ancient=True):
            yield SampleBAM(self.species, pop, sample)

    def output(self):
        return [luigi.LocalTarget('data/angsd/{}-ancient.{}'.format(self.basename, ext)) for ext in
                ['haplo.gz', 'arg', 'list', 'log']]

    def run(self):
        # unpack the params
        bam_files = self.input()
        _, _, list_file, log_file = self.output()

        # make a list of all the BAM files
        with list_file.open('w') as fout:
            for bam_file, _ in bam_files:
                fout.write(bam_file.path + '\n')

        # NOTE using a regions file is crazy slow, so call all sites and extract just the VCF SNPs in PlinkTpedToBed()
        cmd = ['angsd',
               '-nThreads',    self.resources['cpu-cores'],
               '-out',         utils.trim_ext(list_file.path),
               '-dohaplocall', 1,
               '-doCounts',    1,
               '-bam',         list_file.path,
               '-r',           '{}:'.format(self.chrom),
               '-minMapQ',     HARD_MAPQ_CUTOFF,
               '-minQ',        HARD_BASEQ_CUTOFF,
               '-trim',        HARD_CLIP_DIST]

        log = utils.run_cmd(cmd)

        with log_file.open('w') as fout:
            fout.write(log)
