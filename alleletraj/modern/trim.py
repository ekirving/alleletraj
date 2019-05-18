#!/usr/bin/env python
# -*- coding: utf-8 -*-

# third party modules
import luigi

# local modules
from alleletraj import utils
from alleletraj.sra import SraToolsFasterqDump

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
        return SraToolsFasterqDump(self.accession, self.paired)

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
               '--output_dir', 'data/fastq']

        if self.paired:
            cmd.append('--paired')

        # add the fastq files
        for fastq in fastq_files:
            cmd.append(fastq.path)

        # perform the trimming
        utils.run_cmd(cmd)

        # trim_galore does not let us name the outputs so rename the files
        utils.run_cmd(["rename 's/_val_[12]/-trim/' data/fastq/{}_*".format(self.accession)], shell=True)


if __name__ == '__main__':
    luigi.run()
