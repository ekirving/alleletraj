#!/usr/bin/env python
# -*- coding: utf-8 -*-

import luigi
import os

# import my custom modules
from alleletraj import utils
from alleletraj.ancient.alignment import FilterUniqueSAMCons
from alleletraj.modern.alignment import PicardMarkDuplicates, GATK
from alleletraj.reference import PicardSequenceDictionary


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


if __name__ == '__main__':
    luigi.run()
