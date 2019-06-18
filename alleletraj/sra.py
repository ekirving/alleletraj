#!/usr/bin/env python
# -*- coding: utf-8 -*-

# standard modules
import os

# third party modules
import luigi

# local modules
from alleletraj import utils
from alleletraj.const import CPU_CORES_LOW


def entrez_direct_esearch(biosample):
    """
    Get the list of SRA run accessions for a given BioSample code.

    Returns: [(biosample, run_accession, library_layout), ...]

    https://www.ncbi.nlm.nih.gov/books/NBK179288/
    """
    cmd = "esearch -db sra -query 'SAMEA2821680' | " \
          "efetch -format xml | " \
          "xtract -pattern EXPERIMENT_PACKAGE -element SAMPLE@alias RUN@accession " \
          "-block LIBRARY_LAYOUT -if PAIRED -lbl 'paired' -else -lbl 'single'".format(biosample)

    runs = utils.run_cmd([cmd], shell=True)

    return [run.split('\t') for run in runs.strip().split('\n')]


class SraToolsFasterqDump(utils.PipelineTask):
    """
    Fetches the paired-end and single end FASTQ files for a given accession code, using SRA Tools.

    https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump

    :type accession: str
    :type paired: bool
    """
    accession = luigi.Parameter()
    paired = luigi.BoolParameter()

    resources = {'cpu-cores': CPU_CORES_LOW}

    # do not retry after failure, as this just chews CPU cycles
    retry_count = 0

    def output(self):
        if self.paired:
            return [luigi.LocalTarget('data/fastq/{}_{}.fastq.gz'.format(self.accession, pair)) for pair in [1, 2]]
        else:
            return [luigi.LocalTarget('data/fastq/{}.fastq.gz'.format(self.accession))]

    def run(self):
        # use the NCBI SRA toolkit to fetch the fastq files
        cmd = ['fasterq-dump',
               '--threads', self.resources['cpu-cores'],
               '--outdir', 'data/fastq',
               '--temp', 'data/fastq',
               self.accession]

        stderr = utils.run_cmd(cmd)

        for gz_file in self.output():
            fastq_path = utils.trim_ext(gz_file.path)

            if not os.path.isfile(fastq_path):
                # fasterq-dump does not return a valid exit code when there is an error
                raise RuntimeError(stderr)
            else:
                # fasterq-dump does not support the old --gzip flag, so we need to do it manually
                utils.run_cmd(['gzip', fastq_path])


if __name__ == '__main__':
    luigi.run()
