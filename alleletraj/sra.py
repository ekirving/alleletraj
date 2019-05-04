#!/usr/bin/env python
# -*- coding: utf-8 -*-

# third party modules
import luigi

# local modules
from alleletraj import utils
from alleletraj.consts import CPU_CORES_LOW


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
            return [luigi.LocalTarget('data/fastq/{}_{}.fastq.gz'.format(self.accession, pair)) for pair in [1, 2]]
        else:
            return [luigi.LocalTarget('data/fastq/{}.fastq.gz'.format(self.accession))]

    def run(self):
        # use the NCBI SRA toolkit to fetch the fastq files
        utils.run_cmd(['fasterq-dump',
                       '--threads', self.resources['cpu-cores'],
                       '--outdir', './fastq',
                       self.accession])

        # fasterq-dump does not support the old --gzip flag, so we need to do it manually
        for fastq in self.output():
            utils.run_cmd(['gzip', utils.trim_ext(fastq.path)])


if __name__ == '__main__':
    luigi.run()
