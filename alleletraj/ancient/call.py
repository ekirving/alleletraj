#!/usr/bin/env python
# -*- coding: utf-8 -*-

# third party modules
import luigi

# local modules
from alleletraj import utils
from alleletraj.bam import AlignedBAM
from alleletraj.const import CPU_CORES_MED
from alleletraj.modern.vcf import WholeGenomeSNPsVCF

# number of bases to hard clip
HARD_CLIP_DIST = 5

# minimum mapping quality (hard filtered)
HARD_MAPQ_CUTOFF = 30

# minimum base quality (hard filtered)
HARD_BASEQ_CUTOFF = 30


class CallAncientGenotypes(utils.PipelineTask):
    """
    Call ancient genotypes with angsd by randomly selecting a base.

    :type species: str
    """
    species = luigi.Parameter()

    resources = {'cpu-cores': CPU_CORES_MED}

    def requires(self):
        yield WholeGenomeSNPsVCF(self.species)

        for pop, sample in self.all_modern_samples:
            yield AlignedBAM(self.species, pop, sample)

    def output(self):
        return [luigi.LocalTarget('data/angsd/{}-ancient.{}'.format(self.basename, ext)) for ext in
                ['arg', 'haplo.gz', 'pos', 'list']]

    def run(self):
        # unpack the inputs/outputs
        vcf_file, bam_files = self.input()[0], self.input()[1:]
        _, _, pos_file, list_file = self.output()

        # extract the list of sites from the modern VCF
        utils.run_cmd(['bcftools query -f "%CHROM:%POS-%POS\\n" {} > {}'.format(vcf_file.path, pos_file.path)],
                      shell=True)

        # make a list of all the BAM files
        with list_file.open('w') as fout:
            for bam_file, _ in bam_files:
                fout.write(bam_file.path + '\n')

        # NOTE using a regions file is ~15x slower, but we only want the sites that overlap with VCF (inc. non-variant)
        cmd = ['angsd',
               '-nThreads',    self.resources['cpu-cores'],
               '-out',         utils.trim_ext(pos_file),
               '-dohaplocall', 1,
               '-doCounts',    1,
               '-bam',         list_file.path,
               '-rf',          pos_file.path,
               '-minMapQ',     HARD_MAPQ_CUTOFF,
               '-minQ',        HARD_BASEQ_CUTOFF,
               '-trim',        HARD_CLIP_DIST]

        utils.run_cmd([cmd])
