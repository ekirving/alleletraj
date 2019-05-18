#!/usr/bin/env python
# -*- coding: utf-8 -*-

# third party modules
import luigi

# local modules
from alleletraj import utils
from alleletraj.ancient.align import BwaSamSe
from alleletraj.samtools import SAMToolsMerge


class FilterUniqueSAMCons(utils.DatabaseTask):
    """
    Remove PCR duplicates, so we don't overestimate coverage.

    FilterUniqueSAMCons calls consensus on bases in duplicate reads, rather than simply keeping the read with best baseq

    https://bioinf.eva.mpg.de/fastqProcessing/

    NOTE When `accession` is not set, this deduplicates a merged BAM file containing all the sample accessions.

    :type species: str
    :type sample: str
    :type accession: str
    """
    species = luigi.Parameter()
    sample = luigi.Parameter()
    accession = luigi.Parameter(default=None)

    def requires(self):
        if self.accession:
            return BwaSamSe(self.species, self.sample, self.accession, self.accession_data['paired'])
        else:
            return SAMToolsMerge(self.species, self.sample)

    def output(self):
        if self.accession:
            return [luigi.LocalTarget('data/bam/{}.sort.rmdup.{}'.format(self.accession, ext)) for ext in
                    ['bam', 'bam.bai']]
        else:
            return [luigi.LocalTarget('data/bam/{}.merged.rmdup.{}'.format(self.sample, ext)) for ext in
                    ['bam', 'bam.bai']]

    def run(self):
        # unpack the params
        bam_in, _ = self.input()
        bam_out, _ = self.output()

        with bam_out.temporary_path() as bam_path:
            # filter duplicates
            cmd = 'samtools view -h {} | FilterUniqueSAMCons.py | samtools view -b -o {}'.format(bam_in.path, bam_path)

            utils.run_cmd([cmd], shell=True)

        # index the BAM file
        utils.run_cmd(['samtools', 'index', '-b', bam_out.path])


if __name__ == '__main__':
    luigi.run()
