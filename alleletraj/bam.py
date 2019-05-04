#!/usr/bin/env python
# -*- coding: utf-8 -*-

# standard modules
import os

# third party modules
import luigi

# local modules
from alleletraj import utils
from alleletraj.ancient.align import FilterUniqueSAMCons
from alleletraj.modern.align import PicardMarkDuplicates, GATK
from alleletraj.ref import PicardSequenceDictionary


class DeduplicatedBAM(utils.PipelineTask):
    """
    Wrapper task to return a deduplicated BAM file.

    We use different trimming, alignment and deduplication methods depending on ancient/modern status of sample.

    :type species: str
    :type population: str
    :type sample: str
    """
    species = luigi.Parameter()
    population = luigi.Parameter()
    sample = luigi.Parameter()

    def requires(self):
        if self.ancient:
            # use the aDNA pipeline
            return FilterUniqueSAMCons(self.species, self.sample, self.accession)
        else:
            # use the modern pipeline
            return PicardMarkDuplicates(self.species, self.sample, self.accession)

    def output(self):
        return self.input()


class GATKRealignerTargetCreator(utils.PipelineTask):
    """
    Use the GATK RealignerTargetCreator to find suspicious indels to target for local realignment

    :type species: str
    :type sample: str
    :type accession: str
    """
    species = luigi.Parameter()
    sample = luigi.Parameter()
    accession = luigi.Parameter()

    resources = {'cpu-cores': 1, 'ram-gb': 4}

    def requires(self):
        yield DeduplicatedBAM(self.species, self.sample, self.accession)
        yield PicardSequenceDictionary(self.species)

    def output(self):
        return [luigi.LocalTarget('data/bam/{}.sort.rmdup.realign.{}'.format(self.accession, ext)) for ext in
                ['intervals', 'intervals.log']]

    def run(self):
        # unpack the inputs/outputs
        (bam_in, _, _), (ref_file, _, _) = self.input()
        itv_file, log_file = self.output()

        with itv_file.temporary_path() as itv_path, open(log_file.path, 'w') as log_fout:
            utils.run_cmd(['java', self.java_mem,
                           '-jar', GATK,
                           '--analysis_type', 'RealignerTargetCreator',
                           '--reference_sequence', ref_file.path,
                           '--input_file', bam_in.path,
                           '--out', itv_path],
                          stdout=log_fout)


class GATKIndelRealigner(utils.PipelineTask):
    """
    Use the GATK IndelRealigner to perform local realignment of reads around indels

    :type species: str
    :type sample: str
    :type accession: str
    """
    species = luigi.Parameter()
    sample = luigi.Parameter()
    accession = luigi.Parameter()

    resources = {'cpu-cores': 1, 'ram-gb': 4}

    def requires(self):
        yield DeduplicatedBAM(self.species, self.sample, self.accession)
        yield GATKRealignerTargetCreator(self.species, self.sample, self.accession)
        yield PicardSequenceDictionary(self.species)

    def output(self):
        return [luigi.LocalTarget('data/bam/{}.sort.rmdup.realign.{}'.format(self.accession, ext)) for ext in
                ['bam', 'bam.bai', 'log']]

    def run(self):
        # unpack the inputs/outputs
        (bam_in, _, _), (itv_file, _), (ref_file, _, _) = self.input()
        bam_out, bai_out, log_file = self.output()

        with bam_out.temporary_path() as bam_path, open(log_file.path, 'w') as log_fout:
            utils.run_cmd(['java', self.java_mem,
                           '-jar', GATK,
                           '--analysis_type', 'IndelRealigner',
                           '--reference_sequence', ref_file.path,
                           '--input_file', bam_in.path,
                           '--targetIntervals', itv_file.path,
                           '--out', bam_path],
                          stdout=log_fout)

            # GATK automatically creates an index for us, but we need to rename the temp file
            os.rename('{}.bai'.format(bam_path), bai_out)


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
    """
    species = luigi.Parameter()
    population = luigi.Parameter()
    sample = luigi.Parameter()

    def requires(self):

        if self.ancient:
            try:
                # TODO use the provided BAM file
                pass

            except IndexError:
                # TODO rescale baseq for ancient samples
                pass
        else:
            try:
                # use the provided BAM file
                return ExternalBAM(self.all_populations[self.population][self.sample]['path'])

            except IndexError:
                # align our own BAM file
                return SAMToolsMerge(self.species, self.population, self.sample)

    def output(self):
        return self.input()


if __name__ == '__main__':
    luigi.run()
