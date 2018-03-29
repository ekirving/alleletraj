#!/usr/bin/env python
# -*- coding: utf-8 -*-

from db_conn import db_conn

from populate_qtls import *
from populate_samples import populate_samples
from populate_coverage import *
from estimate_allele_freq import *
from discover_snps import discover_snps

SPECIES = ['pig'] #, 'cattle', 'horse']

for species in SPECIES:
    # populate_qtls(species)
    # compute_qtl_windows(species)
    # populate_intervals(species)
    # estimate_allele_freq(species)
    # populate_interval_snps(species)
    # populate_samples(species)
    populate_coverage(species)
    discover_snps(species)
