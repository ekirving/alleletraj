#!/usr/bin/env python
# -*- coding: utf-8 -*-

# third party modules
import luigi

# local modules
from alleletraj import utils
from alleletraj.const import CPU_CORES_MED
from alleletraj.ref import ReferenceFASTA, BwaIndexBWTSW
from alleletraj.sra import SraToolsFastqDump

# hard filters for TrimGalore!
TRIM_MIN_BASEQ = 20
TRIM_MIN_LENGTH = 25


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
            return [luigi.LocalTarget('data/fastq/{}_{}-trim.fq.gz'.format(self.accession, pair)) for pair in [1, 2]]
        else:
            return [luigi.LocalTarget('data/fastq/{}_trimmed.fq.gz'.format(self.accession))]

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
        utils.run_cmd(["rename 's/_val_[12]/-trim/' fastq/{}_*".format(self.accession)], shell=True)


class BwaMem(utils.PipelineTask):
    """
    Align the fastq file(s) to the reference genome using the BWA-MEM algorithm.

    :type species: str
    :type sample: str
    :type accession: str
    :type paired: bool
    """
    species = luigi.Parameter()
    sample = luigi.Parameter()
    accession = luigi.Parameter()
    paired = luigi.BoolParameter(default=True)  # TODO read from CSV

    resources = {'cpu-cores': CPU_CORES_MED, 'ram-gb': 8}

    def requires(self):
        yield TrimGalore(self.accession, self.paired)
        yield ReferenceFASTA(self.species)
        yield BwaIndexBWTSW(self.species)

    def output(self):
        return [luigi.LocalTarget('data/bam/{}.sort.{}'.format(self.accession, ext)) for ext in ['bam', 'bam.bai']]

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
        return [luigi.LocalTarget('data/bam/{}.sort.rmdup.{}'.format(self.accession, ext)) for ext in
                ['bam', 'bam.bai', 'log']]

    def run(self):
        # unpack the params
        bam_in, _ = self.input()
        bam_out, _, log_file = self.output()

        with bam_out.temporary_path() as bam_path:
            utils.run_cmd(['java', self.java_mem,
                           '-jar', 'picard.jar',
                           'MarkDuplicates',
                           'INPUT=' + bam_in.path,
                           'OUTPUT=' + bam_path,
                           'METRICS_FILE=' + log_file.path,
                           'REMOVE_DUPLICATES=true',
                           'QUIET=true'])

        # index the BAM file
        utils.run_cmd(['samtools', 'index', '-b', bam_out.path])


if __name__ == '__main__':
    luigi.run()
