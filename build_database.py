#!/usr/bin/env python
# -*- coding: utf-8 -*-

from populate_coverage import *
from discover_modern_snps import discover_modern_snps
from discover_snps import discover_snps
from analyse_qtls import analyse_qtls
from populate_samples import populate_samples
from ascertainment import perform_ascertainment

# load all the Ensembl data for this species
load_ensembl_genes()
load_ensembl_variants()

# load the SNP Chip data
load_snpchip_variants()

# ascertain SNPs in modern whole genome data
discover_modern_snps()

# load the QTLs from the AnimalQTL database
populate_qtls()

# load psudo-QTLs from other sources
populate_sweeps()
populate_mc1r_locus()
populate_neutral_loci()

# link each QTL to the ascertained modern SNPs
populate_qtl_snps(POPULATION)

# flag the modern SNPs which fall into "neutral" regions
mark_neutral_snps()

# calculate the unique set of non-overlapping genomic loci from the QTLs
populate_intervals()
populate_interval_snps(POPULATION)

# load the sample metadata
populate_samples()

# load the sample reads for each ascertained SNP
populate_sample_reads()

# apply quality filters to the sample reads
discover_snps(POPULATION)

# analyse the coverage and quality for SNPs in each QTLs
analyse_qtls()

# pick the best SNPs to target for a capture array
perform_ascertainment()
