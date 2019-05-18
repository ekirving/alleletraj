#!/usr/bin/env python
# -*- coding: utf-8 -*-

# third party modules
import luigi

# local modules
from alleletraj import utils
from alleletraj.ancient.align import BwaSamSe


class FilterUniqueSAMCons(utils.DatabaseTask):
    """
    Remove PCR duplicates, so we don't overestimate coverage.

    FilterUniqueSAMCons calls consensus on bases in duplicate reads, rather than simply keeping the read with best baseq

    https://bioinf.eva.mpg.de/fastqProcessing/

    :type species: str
    :type sample: str
    :type accession: str
    """
    species = luigi.Parameter()
    sample = luigi.Parameter()
    accession = luigi.Parameter()

    def requires(self):
        return BwaSamSe(self.species, self.sample, self.accession, self.accession_data['paired'])

    def output(self):
        return [luigi.LocalTarget('data/bam/{}.sort.rmdup.{}'.format(self.accession, ext)) for ext in
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
