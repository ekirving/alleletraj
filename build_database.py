#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pipeline_modern_snps import ModernSNPsPipeline
from pipeline_ensembl import EnsemblPipeline
from pipeline_snpchip import SNPChipPipeline

from populate_coverage import *
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


class BuildDatabase(luigi.WrapperTask):
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
        populate_intervals()
        populate_interval_snps(POPULATION)

        # load the sample metadata
        if SPECIES == 'pig':
            populate_pig_samples()

        elif SPECIES == 'horse':
            populate_horse_samples()

        # TODO make samples file w/ gender call for ploidy
        # load the sample reads for each ascertained SNP
        populate_sample_reads()

        # apply quality filters to the sample reads
        discover_snps(POPULATION)

        # analyse the coverage and quality for SNPs in each QTLs
        analyse_qtls()

        if SPECIES == 'pig':
            # pick the best SNPs to target for a capture array
            perform_ascertainment()

        # graph the age of derived alleles
        graph_derived()
