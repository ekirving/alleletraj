#!/usr/bin/env python
# -*- coding: utf-8 -*-

from db_conn import db_conn

from populate_qtls import *
from populate_samples import populate_samples
from populate_coverage import *
from discover_modern_snps import *
from discover_snps import discover_snps
from analyse_qtls import analyse_qtls
from ascertainment import perform_ascertainment

SPECIES = ['pig']  # , 'horse', 'goat']

for species in SPECIES:

    # load_ensembl_genes(species)
    # load_ensembl_variants(species)

    # load_snpchip_variants(species)

    # discover_modern_snps(species)
    mark_neutral_snps()

    # populate_qtls(species)
    # compute_qtl_windows(species)
    # populate_sweeps(species)

    # TODO load MC1R as a QTL region
    # TODO load neutral regions as QTLs

    # populate_intervals(species)
    # populate_interval_snps(species)

    # populate_samples(species)
    # populate_coverage(species)

    # discover_snps(species)
    # analyse_qtls(species)

    # perform_ascertainment(species)