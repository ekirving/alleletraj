#!/usr/bin/env python
# -*- coding: utf-8 -*-

# standard modules
import os

# third party modules
import luigi

# local modules
from alleletraj import utils
from alleletraj.const import CPU_CORES_LOW


def entrez_direct_esearch(bioproject, biosample, fastq_only=True, libname_filter=''):
    """
    Get the list of SRA run accessions for a given BioSample code.

    Returns: [(bioproject, biosample, run_accessions, library_name, library_layout, file_format), ...]

    https://www.ncbi.nlm.nih.gov/books/NBK179288/
    """
    cmd = "esearch -db sra -query '{}' | " \
          "efetch -format xml | " \
          "xtract -pattern EXPERIMENT_PACKAGE " \
          "-block STUDY -def '-' -first EXTERNAL_ID " \
          "-block SAMPLE -def '-' -first EXTERNAL_ID " \
          "-block RUN_SET -def '-' -sep ';' -element RUN@accession " \
          "-block LIBRARY_DESCRIPTOR -def '-' -first LIBRARY_NAME " \
          "-block LIBRARY_LAYOUT -if PAIRED@NOMINAL_LENGTH -lbl 'paired' -else -lbl 'single' " \
          "-block RUN_SET -if AlignInfo -lbl 'bam' -else -lbl 'fastq'".format(biosample)

    result = utils.run_cmd([cmd], shell=True, verbose=False).strip()

    records = []

    if result == '':
        print('WARNING: No data returned from entrez ({}, {})'.format(bioproject, biosample))

    else:
        for line in result.split('\n'):
            record = line.split('\t')

            if len(record) != 6:
                raise RuntimeError('Malformed data from entrez ({}, {}) - {}'.format(bioproject, biosample, record))

            # elif bioproject != record[0] or biosample != record[1]:
            elif biosample != record[1]:
                # make sure the data belongs to the right biosample
                print('WARNING: Record contains incorrect data ({}, {}) - {}'.format(bioproject, biosample, record))

            elif fastq_only and record[5] != 'fastq':
                print('WARNING: Record is not a fastq ({}, {}) - {}'.format(bioproject, biosample, record))

            elif libname_filter and libname_filter in record[3]:
                print('WARNING: Record matches library name filter ({}, {}) - {}'.format(bioproject, biosample, record))
            else:
                # split multiple run accessions
                record[2] = record[2].split(';')

                records.append([None if field == '-' else field for field in record])

    return records


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
