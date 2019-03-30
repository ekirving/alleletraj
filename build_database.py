#!/usr/bin/env python
# -*- coding: utf-8 -*-

from populate_coverage import *
from discover_modern_snps import discover_modern_snps, links_modern_snps
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

# load all the Ensembl data for this species
load_ensembl_genes()
load_ensembl_variants()

# load the SNP Chip data
load_snpchip_variants()

# ascertain SNPs in modern whole genome data
discover_modern_snps()

# link modern SNPs to their dbsnp, gene and snpchip records
links_modern_snps()

# load the QTLs from the AnimalQTL database
populate_qtls()

# load pseudo-QTLs from other sources
populate_sweeps()
populate_mc1r_locus()

if SPECIES == 'pig':
    populate_pig_mummies_loci()

populate_neutral_loci()

# TODO make this work with DOM and DOM2
# link each QTL to the ascertained modern SNPs
populate_qtl_snps(POPULATION)

# flag the modern SNPs which fall into "neutral" regions
mark_neutral_snps()

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
