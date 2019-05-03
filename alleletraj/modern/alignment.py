#!/usr/bin/env python
# -*- coding: utf-8 -*-

import luigi
import os

from multiprocessing import cpu_count

# import my custom modules
from alleletraj import utils
from alleletraj.ensembl.load import DownloadEnsemblData

# how many CPU cores does this machine have
TOTAL_CORES = cpu_count()

# set how many cores a single working can use
CPU_CORES_LOW = int(TOTAL_CORES * 0.1)   # 10%
CPU_CORES_MED = int(TOTAL_CORES * 0.25)  # 25%
CPU_CORES_HIGH = int(TOTAL_CORES * 0.5)  # 50%
CPU_CORES_MAX = int(TOTAL_CORES * 0.9)   # 90%

# hard filters for TrimGalore!
TRIM_MIN_BASEQ = 20
TRIM_MIN_LENGTH = 25

# location of software tools
GATK = "/usr/local/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar"
PICARD = "/usr/local/picard-tools-2.5.0/picard.jar"


class SraToolsFastqDump(utils.PipelineTask):
    """
    Fetches the paired-end and single end FASTQ files for a given accession code, using SRA Tools

    :type accession: str
    :type paired: bool
    """
    accession = luigi.Parameter()
    paired = luigi.BoolParameter()

    resources = {'cpu-cores': CPU_CORES_LOW}

    def output(self):
        if self.paired:
            return [luigi.LocalTarget('fastq/{}_{}.fastq.gz'.format(self.accession, pair)) for pair in [1, 2]]
        else:
            return [luigi.LocalTarget('fastq/{}.fastq.gz'.format(self.accession))]

    def run(self):
        # use the NCBI SRA toolkit to fetch the fastq files
        utils.run_cmd(['fasterq-dump',
                       '--threads', self.resources['cpu-cores'],
                       '--outdir', './fastq',
                       self.accession])

        # fasterq-dump does not support the old --gzip flag, so we need to do it manually
        for fastq in self.output():
            utils.run_cmd(['gzip', utils.trim_ext(fastq.path)])


class TrimGalore(utils.PipelineTask):
    """
    Trim adapters, low-quality bases from the 3' end, and drop any resulting reads below a min-length threshold.

    :type accession: str
    :type paired: bool
    """
    accession = luigi.Parameter()
    paired = luigi.BoolParameter()

    def requires(self):
        return SraToolsFastqDump(self.accession, self.paired)

    def output(self):
        if self.paired:
            return [luigi.LocalTarget('fastq/{}_{}_trim.fq.gz'.format(self.accession, pair)) for pair in [1, 2]]
        else:
            return [luigi.LocalTarget('fastq/{}_trimmed.fq.gz'.format(self.accession))]

    def run(self):
        fastq_files = self.input()

        cmd = ['trim_galore',
               '--quality', TRIM_MIN_BASEQ,   # trim low-quality ends from reads in addition to adapter removal
               '--length',  TRIM_MIN_LENGTH,  # discard reads that became shorter than length because of trimming
               '--fastqc',                    # run FastQC once trimming is complete
               '--output_dir', 'fastq']

        if self.paired:
            cmd.append('--paired')

        # add the fastq files
        for fastq in fastq_files:
            cmd.append(fastq.path)

        # perform the trimming
        utils.run_cmd(cmd)

        # trim_galore does not let us name the outputs so rename the files
        utils.run_cmd(["rename 's/_val_[12]/_trim/' fastq/{}_*".format(self.accession)], shell=True)


class ReferenceFASTA(utils.PipelineTask):
    """
    Get the reference genome and index it.

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):
        return DownloadEnsemblData(self.species, 'fasta', bgzip=True)

    def output(self):
        return [luigi.LocalTarget('ensembl/{}.{}.dna.toplevel.{}'.format(self.binomial, self.assembly, ext)) for ext in
                ['fa.gz', 'fa.gz.fai']]

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
        return [luigi.LocalTarget('{}.{}'.format(ref_file.path, ext)) for ext in ['amb', 'ann', 'bwt', 'pac', 'sa']]

    def run(self):
        ref_file, _ = self.input()

        utils.run_cmd(['bwa',
                       'index',        # index needed for bwa alignment
                       '-a', 'bwtsw',  # algorithm suitable for mammals
                       ref_file.path])


class BwaMem(utils.PipelineTask):
    """
    Align a fastq file to the reference genome

    :type species: str
    :type sample: str
    :type accession: str
    :type paired: bool
    """
    species = luigi.Parameter()
    sample = luigi.Parameter()
    accession = luigi.Parameter()
    paired = luigi.BoolParameter(default=True)

    resources = {'cpu-cores': CPU_CORES_MED, 'ram-gb': 8}

    def requires(self):
        yield TrimGalore(self.accession, self.paired)
        yield ReferenceFASTA(self.species)
        yield BwaIndexBWTSW(self.species)

    def output(self):
        return [luigi.LocalTarget('bam/{}.sort.{}'.format(self.accession, ext)) for ext in ['bam', 'bam.bai']]

    def run(self):
        # unpack the input params
        fastq_input, (ref_file, _), _ = self.input()
        bam_out, _ = self.output()

        params = {
            # CPU threads to use
            'threads': self.resources['cpu-cores'],

            # compose the read group metadata
            'readgroup': r'@RG\tID:{sample}\tSM:{sample}'.format(sample=self.sample),

            # reference genome
            'reference': ref_file.path,

            # get the fastq file(s) to align
            'fastq': ' '.join(fastq.path for fastq in fastq_input),
        }

        with bam_out.temporary_path() as bam_path:

            # get the temporary path for the bam file
            params['bam'] = bam_path

            # align using the mem algorithm, and convert to a sorted BAM file
            cmd = "bwa mem -t {threads} -R '{readgroup}' {reference} {fastq} " \
                  " | samtools sort -@ {threads} -O bam -o {bam} -".format(**params)

            # perform the alignment
            utils.run_cmd([cmd], shell=True)

        # index the BAM file
        utils.run_cmd(['samtools', 'index', '-b', bam_out.path])


class PicardMarkDuplicates(utils.PipelineTask):
    """
    Remove PCR duplicates, so we don't overestimate coverage

    :type species: str
    :type sample: str
    :type accession: str
    """
    species = luigi.Parameter()
    sample = luigi.Parameter()
    accession = luigi.Parameter()

    resources = {'cpu-cores': 1, 'ram-gb': 8}

    def requires(self):
        return BwaMem(self.species, self.sample, self.accession)

    def output(self):
        return [luigi.LocalTarget('bam/{}.sort.rmdup.{}'.format(self.accession, ext)) for ext in
                ['bam', 'bam.bai', 'log']]

    def run(self):
        # unpack the params
        bam_in, _ = self.input()
        bam_out, _, log_file = self.output()

        with bam_out.temporary_path() as bam_path:

            utils.run_cmd(['java', self.java_mem,
                           '-jar', PICARD,
                           'MarkDuplicates',
                           'INPUT=' + bam_in.path,
                           'OUTPUT=' + bam_path,
                           'METRICS_FILE=' + log_file.path,
                           'REMOVE_DUPLICATES=true',
                           'QUIET=true'])

        # index the BAM file
        utils.run_cmd(['samtools', 'index', '-b', bam_out.path])


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
        return [luigi.LocalTarget('ensembl/{}.{}.dna.toplevel.{}'.format(self.binomial, self.assembly, ext)) for ext in
                ['fa', 'fa.fai', 'dict']]

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
                           '-jar', PICARD,
                           'CreateSequenceDictionary',
                           'R=' + ref_out.path,
                           'O=' + dict_path])


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
        yield PicardMarkDuplicates(self.species, self.sample, self.accession)
        yield PicardSequenceDictionary(self.species)

    def output(self):
        return [luigi.LocalTarget('bam/{}.sort.rmdup.realign.{}'.format(self.accession, ext)) for ext in
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
        yield PicardMarkDuplicates(self.species, self.sample, self.accession)
        yield GATKRealignerTargetCreator(self.species, self.sample, self.accession)
        yield PicardSequenceDictionary(self.species)

    def output(self):
        return [luigi.LocalTarget('bam/{}.sort.rmdup.realign.{}'.format(self.accession, ext)) for ext in
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
        return [luigi.LocalTarget('bam/{}.sort.rmdup.realign.{}'.format(self.sample, ext)) for ext in
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
        try:
            # use the provided BAM file
            return ExternalBAM(self.all_populations[self.population][self.sample]['path'])

        except IndexError:
            # align our own BAM file
            return SAMToolsMerge(self.species, self.population, self.sample)

    def output(self):
        return self.input()


class AlignmentPipeline(utils.PipelineWrapperTask):
    """
    Get BAM files for all samples

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):
        for pop, sample in self.all_samples:
            yield AlignedBAM(self.species, pop, sample)


if __name__ == '__main__':
    luigi.run()
