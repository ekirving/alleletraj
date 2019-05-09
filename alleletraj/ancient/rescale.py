#!/usr/bin/env python
# -*- coding: utf-8 -*-

# third party modules
import luigi

# local modules
from alleletraj import utils
from alleletraj.gatk import GATKIndelRealigner
from alleletraj.ref import ReferenceFASTA

# how many reads to use to estimate the damage frequency
MAPDAMAGE_DOWNSAMPLE = 100000


class MapDamageRescale(utils.PipelineTask):
    """
    Run MapDamage to check for aDNA damage patterns, and rescale baseq scores at read ends.

    :type species: str
    :type sample: str
    :type accession: str
    """
    species = luigi.Parameter()
    sample = luigi.Parameter()
    accession = luigi.Parameter()

    def requires(self):
        yield ReferenceFASTA(self.species)
        yield GATKIndelRealigner(self.species, self.sample, self.accession)

    def output(self):
        return [luigi.LocalTarget('data/bam/{}.sort.rmdup.realign.rescale.{}'.format(self.accession, ext)) for ext in
                ['bam', 'bam.bai', 'log']]

    def run(self):
        # unpack the inputs/outputs
        (ref_file, _), (bam_in, _, _) = self.input()
        bam_out, _, log_file = self.output()

        with bam_out.temporary_path() as bam_path:

            log = utils.run_cmd(['mapDamage',
                                 '-i', bam_in.path,
                                 '-r', ref_file.path,
                                 '-d', 'data/mapdamage/{}'.format(self.sample),
                                 '-n', MAPDAMAGE_DOWNSAMPLE,
                                 '--merge-reference-sequences',
                                 '--rescale',
                                 '--rescale-out={}'.format(bam_path)])

        # index the BAM file
        utils.run_cmd(['samtools', 'index', '-b', bam_out.path])

        with log_file.open('w') as fout:
            fout.write(log)
