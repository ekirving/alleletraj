#!/usr/bin/env python
# -*- coding: utf-8 -*-

# third party modules
import luigi

# local modules
from alleletraj import utils
from alleletraj.gatk import GATKIndelRealigner
from alleletraj.ancient.rescale import MapDamageRescale


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

    :type species: str
    :type population: str
    :type sample: str
    """
    species = luigi.Parameter()
    population = luigi.Parameter()
    sample = luigi.Parameter()

    resources = {'cpu-cores': 1, 'ram-gb': 4}

    # do not retry after failure, as this just chews CPU cycles
    retry_count = 0

    def requires(self):
        return ExternalBAM(self.all_populations[self.population][self.sample]['path'])

    def output(self):
        return [luigi.LocalTarget('data/bam/{}.{}'.format(self.sample, ext)) for ext in ['log', 'errs']]

    def run(self):
        bam_file, _ = self.input()
        log_file, errs_file = self.output()

        # validate the BAM file
        with log_file.temporary_path() as log_path:
            utils.run_cmd(['java', self.java_mem,
                           '-jar', 'picard.jar',
                           'ValidateSamFile',
                           'MODE=SUMMARY',
                           'IGNORE=MATE_NOT_FOUND',
                           'INPUT=' + bam_file.path,
                           'OUTPUT=' + errs_file.path], stderr=open(log_path, 'w'))


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
        for accession in self.all_populations[self.population][self.sample]['accessions']:
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

    def requires(self):
        # is there a path defined in the CSV file
        if 'path' in self.all_populations[self.population][self.sample]:
            # we need to check that the provided file is valid
            return ValidateBamFile(self.species, self.population, self.sample)
        else:
            # if not, then we need to align our own BAM file
            return SAMToolsMerge(self.species, self.population, self.sample)

    def output(self):
        return self.input()


if __name__ == '__main__':
    luigi.run()
