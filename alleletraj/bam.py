#!/usr/bin/env python
# -*- coding: utf-8 -*-

# third party modules
import luigi

# local modules
from alleletraj import utils
from alleletraj.ancient.rescale import MapDamageRescale
from alleletraj.gatk import GATKIndelRealigner


class ExternalBAM(utils.PipelineExternalTask):
    """
    External task dependency for an aligned BAM file.

    N.B. These have been created outside the workflow of this pipeline.

    :type path: str
    """
    path = luigi.Parameter()

    def output(self):
        yield luigi.LocalTarget(self.path)
        yield luigi.LocalTarget(self.path + '.bai')


class ValidateBamFile(utils.PipelineTask):
    """
    Validate an external BAM file using Picard.

    This tool reports on the validity of a BAM file relative to the SAM format specification.

    https://software.broadinstitute.org/gatk/documentation/tooldocs/4.1.2.0/picard_sam_ValidateSamFile.php
    https://software.broadinstitute.org/gatk/documentation/article.php?id=7571

    :type species: str
    :type population: str
    :type sample: str
    """
    species = luigi.Parameter()
    population = luigi.Parameter()
    sample = luigi.Parameter()

    resources = {'cpu-cores': 1, 'ram-gb': 8}

    # do not retry after failure, as this just chews CPU cycles
    retry_count = 0

    def requires(self):
        return ExternalBAM(self.all_modern_data[self.population][self.sample]['path'])

    def output(self):
        # included the bam and bai files in the output
        return list(self.input()) + [luigi.LocalTarget('data/bam/{}.{}'.format(self.sample, ext)) for ext in
                                     ['log', 'err']]

    def run(self):
        bam_file, _ = self.input()
        _, _, log_file, err_file = self.output()

        # validate the BAM file
        with log_file.temporary_path() as log_path:
            utils.run_cmd(['java', self.java_mem,
                           '-jar', 'jar/picard.jar',
                           'ValidateSamFile',
                           'MODE=SUMMARY',
                           'IGNORE=MATE_NOT_FOUND',
                           'INPUT=' + bam_file.path,
                           'OUTPUT=' + err_file.path], stderr=open(log_path, 'w'))


class SAMToolsMerge(utils.PipelineTask):
    """
    Merge multiple libraries into a single BAM file

    :type species: str
    :type population: str
    :type sample: str
    """
    species = luigi.Parameter()
    population = luigi.Parameter()
    sample = luigi.Parameter()

    resources = {'cpu-cores': 1, 'ram-gb': 8}

    def requires(self):
        for accession in self.all_modern_data[self.population][self.sample]['accessions']:
            if self.ancient:
                yield MapDamageRescale(self.species, self.sample, accession)
            else:
                yield GATKIndelRealigner(self.species, self.sample, accession)

    def output(self):
        suffix = 'sort.rmdup.realign.rescale' if self.ancient else 'sort.rmdup.realign'
        return [luigi.LocalTarget('data/bam/{}.{}.{}'.format(self.sample, suffix, ext)) for ext in
                ['bam', 'bam.bai']]

    def run(self):
        bam_inputs = [bam_file.path for bam_file, _, _ in self.input()]
        bam_out, _ = self.output()

        with bam_out.temporary_path() as bam_path:
            utils.run_cmd(['samtools',
                           'merge',
                           bam_path
                           ] + bam_inputs)

        # index the BAM file
        utils.run_cmd(['samtools', 'index', '-b', bam_out.path])


class AlignedBAM(utils.PipelineTask):
    """
    Align raw reads to a reference genome, deduplicate, realign indels, rescale bases, and merge multiple libraries.

    :type species: str
    :type population: str
    :type sample: str
    """
    species = luigi.Parameter()
    population = luigi.Parameter()
    sample = luigi.Parameter()

    def complete(self):
        # this task is complete if our input task is complete
        return self.requires().complete()

    def requires(self):
        if self.all_modern_data[self.population][self.sample].get('path', '') != '':
            # we need to validate any external BAM files before using them
            return ValidateBamFile(self.species, self.population, self.sample)
        else:
            if self.ancient and self.species == 'goat':
                # TODO ancient goat libraries were built from the same PCR reactions, so they need deduping after merge
                pass
            else:
                return SAMToolsMerge(self.species, self.population, self.sample)

    def output(self):
        # only pass on the bam and bai files (i.e. trim off any .log and .err files from ValidateBamFile)
        return self.input()[:2]


class ValidateModernBAMs(utils.PipelineWrapperTask):
    """
    Wrapper taks to validate all the external BAM files.

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):
        for pop, sample in self.all_modern_samples:
            if self.all_modern_data[pop][sample].get('path', '') != '':
                yield ValidateBamFile(self.species, pop, sample)


if __name__ == '__main__':
    luigi.run()
