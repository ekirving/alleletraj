#!/usr/bin/env python
# -*- coding: utf-8 -*-

# third party modules
import luigi

# local modules
from alleletraj import utils
from alleletraj.sra import SraToolsFasterqDump

# hard filters for AdapterRemoval
TRIM_MIN_BASEQ = 20
TRIM_MIN_LENGTH = 25


class AdapterRemoval(utils.PipelineTask):
    """
    Trim adapters, low-quality bases from the 5' and 3' ends, drop any resulting reads below a min-length threshold,
    and collapse overlapping read pairs.

    Becase we collapse paired-end reads, this task always produces a single fastq file.

    :type accession: str
    :type paired: bool
    """
    accession = luigi.Parameter()
    paired = luigi.BoolParameter()

    def requires(self):
        return SraToolsFasterqDump(self.accession, self.paired)

    def output(self):
        yield luigi.LocalTarget('data/fastq/{}-trim.fq.gz'.format(self.accession))
        yield luigi.LocalTarget('data/fastq/{}-trim-discard.fq.gz'.format(self.accession))
        yield luigi.LocalTarget('data/fastq/{}-trim.log'.format(self.accession))

    def run(self):
        fastq_files = self.input()
        fastq_keep, fastq_drop, log_file = self.output()

        cmd = [
            'AdapterRemoval',
            '--discarded',  fastq_drop.path,
            '--settings',   log_file.path,
            '--trimns',                       # trim ambiguous bases (N) at both 5'/3' termini
            '--trimqualities',                # trim bases at both 5'/3' termini with quality scores...
            '--minquality', TRIM_MIN_BASEQ,   # below this value
            '--minlength',  TRIM_MIN_LENGTH,  # discard reads shorter than this length following trimming
            '--gzip'
        ]

        if self.paired:
            # TODO still need to test paired end mode works correctly
            cmd += [
                '--file1', fastq_files[0].path,
                '--file2', fastq_files[1].path,
                '--collapse',
                '--outputcollapsed', fastq_keep.path
            ]
        else:
            cmd += [
                '--basename', utils.trim_ext(log_file.path),
                '--file1',    fastq_files[0].path,
                '--output1',  fastq_keep.path
            ]

        # perform the trimming
        utils.run_cmd(cmd)


if __name__ == '__main__':
    luigi.run()
