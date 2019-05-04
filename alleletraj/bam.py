#!/usr/bin/env python
# -*- coding: utf-8 -*-

import luigi

# import my custom modules
from alleletraj import utils
from alleletraj.gatk import GATKIndelRealigner


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
            yield GATKIndelRealigner(self.species, self.sample, accession)

    def output(self):
        return [luigi.LocalTarget('data/bam/{}.sort.rmdup.realign.{}'.format(self.sample, ext)) for ext in
                ['bam', 'bam.bai']]

    def run(self):
        bam_inputs = [bam_file.path for bam_file, bai_file, log_file in self.input()]
        bam_out, _ = self.output()

        with bam_out.temporary_path() as bam_path:

            utils.run_cmd(['samtools',
                           'merge',
                           bam_path
                           ] + bam_inputs)

        # index the BAM file
        utils.run_cmd(['samtools', 'index', '-b', bam_out.path])


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


class AlignedBAM(utils.PipelineTask):
    """
    Align raw reads to a reference genome, deduplicate, realign indels, and merge multiple libraries.

    :type species: str
    :type population: str
    :type sample: str
    :type ancient: bool
    """
    species = luigi.Parameter()
    population = luigi.Parameter()
    sample = luigi.Parameter()

    def requires(self):
        try:
            # use the provided BAM file
            return ExternalBAM(self.all_populations[self.population][self.sample]['path'])

        except IndexError:
            # align our own BAM file
            return SAMToolsMerge(self.species, self.population, self.sample)

    def output(self):
        return self.input()