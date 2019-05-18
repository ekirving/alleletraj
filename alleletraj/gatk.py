#!/usr/bin/env python
# -*- coding: utf-8 -*-

# standard modules
import os

# third party modules
import luigi

# local modules
from alleletraj import utils
from alleletraj.ref import PicardSequenceDictionary


class GATKRealignerTargetCreator(utils.PipelineTask):
    """
    Use the GATK RealignerTargetCreator to find suspicious indels to target for local realignment.

    :type species: str
    :type sample: str
    :type accession: str
    """
    species = luigi.Parameter()
    sample = luigi.Parameter()
    accession = luigi.Parameter()

    resources = {'cpu-cores': 1, 'ram-gb': 4}

    def requires(self):
        # import locally to avoid recursive import error
        from alleletraj.bam import DeduplicatedBAM

        yield DeduplicatedBAM(self.species, self.sample, self.accession)
        yield PicardSequenceDictionary(self.species)

    def output(self):
        return [luigi.LocalTarget('data/bam/{}.sort.rmdup.realign.{}'.format(self.accession, ext)) for ext in
                ['intervals', 'intervals.log']]

    def run(self):
        # unpack the inputs/outputs
        (bam_in, _), (ref_file, _, _) = self.input()
        itv_file, log_file = self.output()

        with itv_file.temporary_path() as itv_path, open(log_file.path, 'w') as log_fout:
            utils.run_cmd(['java', self.java_mem,
                           '-jar', 'jar/GenomeAnalysisTK.jar',
                           '--analysis_type', 'RealignerTargetCreator',
                           '--reference_sequence', ref_file.path,
                           '--input_file', bam_in.path,
                           '--out', itv_path],
                          stdout=log_fout)


class GATKIndelRealigner(utils.PipelineTask):
    """
    Use the GATK IndelRealigner to perform local realignment of reads around indels.

    :type species: str
    :type sample: str
    :type accession: str
    """
    species = luigi.Parameter()
    sample = luigi.Parameter()
    accession = luigi.Parameter()

    resources = {'cpu-cores': 1, 'ram-gb': 4}

    def requires(self):
        # import locally to avoid recursive import error
        from alleletraj.bam import DeduplicatedBAM

        yield DeduplicatedBAM(self.species, self.sample, self.accession)
        yield GATKRealignerTargetCreator(self.species, self.sample, self.accession)
        yield PicardSequenceDictionary(self.species)

    def output(self):
        return [luigi.LocalTarget('data/bam/{}.sort.rmdup.realign.{}'.format(self.accession, ext)) for ext in
                ['bam', 'bam.bai', 'log']]

    def run(self):
        # unpack the inputs/outputs
        (bam_in, _), (itv_file, _), (ref_file, _, _) = self.input()
        bam_out, bai_out, log_file = self.output()

        with bam_out.temporary_path() as bam_path, open(log_file.path, 'w') as log_fout:
            utils.run_cmd(['java', self.java_mem,
                           '-jar', 'jar/GenomeAnalysisTK.jar',
                           '--analysis_type', 'IndelRealigner',
                           '--reference_sequence', ref_file.path,
                           '--input_file', bam_in.path,
                           '--targetIntervals', itv_file.path,
                           '--out', bam_path],
                          stdout=log_fout)

            # GATK automatically creates an index for us, but we need to rename the temp file
            os.rename('{}.bai'.format(bam_path), bai_out.path)
