#!/usr/bin/env python
# -*- coding: utf-8 -*-

# third party modules
import luigi

# local modules
from alleletraj import utils


class SAMToolsMerge(utils.DatabaseTask):
    """
    Merge multiple libraries into a single BAM file.

    :type species: str
    :type sample: str
    """
    species = luigi.Parameter()
    sample = luigi.Parameter()

    resources = {'cpu-cores': 1, 'ram-gb': 8}

    def requires(self):
        # import locally to avoid recursive import error
        from alleletraj.bam import AccessionBAM

        for accession in self.list_accessions():
            yield AccessionBAM(self.species, self.sample, accession)

    def output(self):
        return [luigi.LocalTarget('data/bam/{}.merged.{}'.format(self.sample, ext)) for ext in ['bam', 'bam.bai']]

    def run(self):
        bam_inputs = [bam_file.path for bam_file, _, _ in self.input()]
        bam_out, _ = self.output()

        with bam_out.temporary_path() as bam_path:
            utils.run_cmd(['samtools', 'merge', bam_path] + bam_inputs)

        # index the BAM file
        utils.run_cmd(['samtools', 'index', '-b', bam_out.path])


if __name__ == '__main__':
    luigi.run()
