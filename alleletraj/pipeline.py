#!/usr/bin/env python
# -*- coding: utf-8 -*-

# third party modules
import luigi

# local modules
from alleletraj import utils
from alleletraj.ancient.snps import AncientSNPsPipeline
from alleletraj.ensembl.link import EnsemblLinkPipeline
from alleletraj.modern.demog import DadiPipeline
from alleletraj.modern.snps import ModernSNPsPipeline
from alleletraj.plink import PlinkMergeBeds
from alleletraj.qtl.analyse import AnalyseQTLsPipeline
from alleletraj.qtl.load import QTLPipeline
from alleletraj.snpchip.link import SNPChipLinkPipeline


class RunAll(utils.PipelineWrapperTask):
    """
    Process all the pipeline tasks.

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):

        # ascertain SNPs in modern whole genome data
        yield ModernSNPsPipeline(self.species)

        # link all the Ensembl data to the modern SNPs
        yield EnsemblLinkPipeline(self.species)

        # link the SNPChip data to the modern SNPs
        yield SNPChipLinkPipeline(self.species)

        # find the best fitting ∂a∂i model
        yield DadiPipeline(self.species)

        # link QTLs and modern snps
        yield QTLPipeline(self.species)

        # load the ancient sample reads for each ascertained SNP
        yield AncientSNPsPipeline(self.species)

        # analyse the coverage and quality for SNPs in each QTLs
        yield AnalyseQTLsPipeline(self.species)

        # TODO replace with analysis pipeline when fully tested
        yield PlinkMergeBeds(self.species)

        # TODO run `selection` on all the 'best' QTL SNPs

        # TODO make all the plots


if __name__ == '__main__':
    luigi.run()
