#!/usr/bin/env python
# -*- coding: utf-8 -*-

import luigi

from pipeline_utils import PipelineWrapperTask

from pipeline_modern_snps import ModernSNPsPipeline
from pipeline_ensembl import EnsemblPipeline
from pipeline_snpchip import SNPChipPipeline
from pipeline_qtls import QTLPipeline
from pipeline_sample_reads import PopulateIntervals, PopulateIntervalSNPs

# from pipeline_sample_reads import *
from discover_snps import discover_snps
from analyse_qtls import analyse_qtls
from populate_samples import populate_pig_samples, populate_horse_samples
from ascertainment import perform_ascertainment
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

        # calculate the unique set of non-overlapping genomic loci from the QTLs
        for chrom in self.chromosomes:
            yield PopulateIntervals(self.species, chrom)

            for pop in self.populations:
                yield PopulateIntervalSNPs(self.species, pop, chrom)

        # load the sample metadata
        if self.species == 'pig':
            populate_pig_samples()

        elif self.species == 'horse':
            populate_horse_samples()

        # TODO make samples file w/ gender call for ploidy
        # load the sample reads for each ascertained SNP
        yield PopulateSampleReads()

        # apply quality filters to the sample reads
        for pop in self.populations:
            yield discover_snps(self.species, pop)

        # analyse the coverage and quality for SNPs in each QTLs
        analyse_qtls()

        if self.species == 'pig':
            # pick the best SNPs to target for a capture array
            perform_ascertainment()

        # graph the age of derived alleles
        graph_derived()
