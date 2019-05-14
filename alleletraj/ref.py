#!/usr/bin/env python
# -*- coding: utf-8 -*-

# third party modules
import luigi

# local modules
from alleletraj import utils
from alleletraj.ensembl.load import DownloadEnsemblData


class ReferenceFASTA(utils.PipelineTask):
    """
    Get the reference genome and index it.

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):
        return DownloadEnsemblData(self.species, 'fasta', bgzip=True)

    def output(self):
        return [luigi.LocalTarget('data/ensembl/{}.{}.dna.toplevel.{}'.format(self.binomial, self.assembly, ext))
                for ext in ['fa.gz', 'fa.gz.fai']]

    def run(self):
        # get the downloaded reference assembly
        ref_file = self.input()

        # build an index
        utils.run_cmd(['samtools', 'faidx', ref_file.path])


class BwaIndexBWTSW(utils.PipelineTask):
    """
    Builds the BWA index for the reference genome, needed for performing alignments

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):
        return ReferenceFASTA(self.species)

    def output(self):
        ref_file, _ = self.input()
        return [luigi.LocalTarget('data/ensembl/{}.{}'.format(ref_file.path, ext)) for ext in
                ['amb', 'ann', 'bwt', 'pac', 'sa']]

    def run(self):
        ref_file, _ = self.input()

        utils.run_cmd(['bwa',
                       'index',        # index needed for bwa alignment
                       '-a', 'bwtsw',  # algorithm suitable for mammals
                       ref_file.path])


class PicardSequenceDictionary(utils.PipelineTask):
    """
    Unzip the reference genome and create the sequence dictionary, because GATK is stupid and cannot handle bgzip.

    :type species: str
    """
    species = luigi.Parameter()

    resources = {'cpu-cores': 1, 'ram-gb': 8}

    def requires(self):
        return ReferenceFASTA(self.species)

    def output(self):
        return [luigi.LocalTarget('data/ensembl/{}.{}.dna.toplevel.{}'.format(self.binomial, self.assembly, ext))
                for ext in ['fa', 'fa.fai', 'dict']]

    def run(self):
        ref_in, _ = self.input()
        ref_out, _, dict_file = self.output()

        # unzip the reference genome
        utils.run_cmd(['gunzip', '--keep', '--force', ref_in.path])

        # build a regular index
        utils.run_cmd(['samtools', 'faidx', ref_out.path])

        # create the sequence dictionary
        with dict_file.temporary_path() as dict_path:
            utils.run_cmd(['java', self.java_mem,
                           '-jar', 'jar/picard.jar',
                           'CreateSequenceDictionary',
                           'R=' + ref_out.path,
                           'O=' + dict_path])


if __name__ == '__main__':
    luigi.run()
