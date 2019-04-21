#!/usr/bin/env python
# -*- coding: utf-8 -*-

import luigi

from pipeline_utils import PipelineWrapperTask

from pipeline_modern_snps import ModernSNPsPipeline
from pipeline_ensembl import EnsemblPipeline
from pipeline_snpchip import SNPChipPipeline
from pipeline_qtls import QTLPipeline
from pipeline_sample_reads import SampleReadsPipeline

# from pipeline_sample_reads import *
from pipeline_discover_snps import DiscoverSNPs
from pipeline_analyse_qtls import AnalyseQTLs
from populate_samples import populate_pig_samples, populate_horse_samples
from pipeline_ascertainment import PerformAscertainment
from graph_derived import graph_derived

# TODO confirm InnoDB tweaks / https://dev.mysql.com/doc/refman/8.0/en/converting-tables-to-innodb.html
# TODO refactor to use snakemake / parallelise the SNP calling part
# TODO chrom is top level param / use for threading pipeline

# TODO in table `sample_reads` replace (interval_id ?, chrom, site) with modsnp_id
# TODO increase hard threshold to >= 30 / or / delete where called = 0

# TODO map samples / indicate time / don't want to confuse geography for time
# TODO have a look at the "pulse" model in Loog paper


class BuildDatabase(PipelineWrapperTask):
    """
    Build the species database

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):

        # ascertain SNPs in modern whole genome data
        yield ModernSNPsPipeline(self.species)

        # link all the Ensembl data
        yield EnsemblPipeline(self.species)

        # link the SNPChip data
        yield SNPChipPipeline(self.species)

        # load the QTLs from the AnimalQTL database, and other regions of interest
        yield QTLPipeline(self.species)

        # load the sample metadata
        if self.species == 'pig':
            populate_pig_samples()

        elif self.species == 'horse':
            populate_horse_samples()

        # load the sample reads for each ascertained SNP
        yield SampleReadsPipeline()

        # apply quality filters to the sample reads
        for pop in self.populations:
            DiscoverSNPs(self.species, pop)

        # analyse the coverage and quality for SNPs in each QTLs
        yield AnalyseQTLs(self.species)

        if self.species == 'pig':
            # pick the best SNPs to target for a capture array
            yield PerformAscertainment(self.species)

        # graph the age of derived alleles
        graph_derived()
