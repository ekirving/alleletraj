#!/usr/bin/env python
# -*- coding: utf-8 -*-

# third party modules
import luigi

# local modules
from alleletraj import utils
from alleletraj.modern.align import BwaMem
from alleletraj.samtools import SAMToolsMerge


class PicardMarkDuplicates(utils.DatabaseTask):
    """
    Remove PCR duplicates, so we don't overestimate coverage.

    MarkDuplicates groups duplicate reads by their 5' alignment and keeps the read with the highest total base qaulity.

    https://software.broadinstitute.org/gatk/documentation/tooldocs/4.1.2.0/picard_sam_markduplicates_MarkDuplicates.php

    NOTE When `accession` is not set, this deduplicates a merged BAM file containing all the sample accessions.

    :type species: str
    :type sample: str
    :type accession: str
    """
    species = luigi.Parameter()
    sample = luigi.Parameter()
    accession = luigi.Parameter(default=None)

    resources = {'cpu-cores': 1, 'ram-gb': 8}

    def requires(self):
        if self.accession:
            return BwaMem(self.species, self.sample, self.accession, self.accession_data['paired'])
        else:
            return SAMToolsMerge(self.species, self.sample)

    def output(self):
        if self.accession:
            return [luigi.LocalTarget('data/bam/{}.sort.rmdup.{}'.format(self.accession, ext)) for ext in
                    ['bam', 'bam.bai', 'log']]
        else:
            return [luigi.LocalTarget('data/bam/{}.merged.rmdup.{}'.format(self.sample, ext)) for ext in
                    ['bam', 'bam.bai', 'log']]

    def run(self):
        # unpack the params
        bam_in, _ = self.input()
        bam_out, _, log_file = self.output()

        # TODO consider switching with MarkDuplicatesWithMateCigar
        # https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.4.0/picard_sam_markduplicates_MarkDuplicatesWithMateCigar.php
        with bam_out.temporary_path() as bam_path:
            utils.run_cmd(['java', self.java_mem,
                           '-jar', 'jar/picard.jar',
                           'MarkDuplicates',
                           'INPUT=' + bam_in.path,
                           'OUTPUT=' + bam_path,
                           'METRICS_FILE=' + log_file.path,
                           'REMOVE_DUPLICATES=true',
                           'QUIET=true'])

        utils.run_cmd(['samtools', 'index', '-b', bam_out.path])


if __name__ == '__main__':
    luigi.run()
