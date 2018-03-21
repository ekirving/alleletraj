#!/usr/bin/env python
# -*- coding: utf-8 -*-

from db_conn import db_conn

from populate_qtls import populate_qtls
from populate_samples import populate_samples
from populate_coverage import *
from discover_snps import discover_snps

SPECIES = ['pig'] #, 'cattle', 'horse']

for species in SPECIES:
    # populate_qtls(species)
    # populate_intervals(species)
    # populate_samples(species)
    # populate_coverage(species)

    for dataset in ['sample_reads_innodb',
                    'sample_reads_innodb_part',
                    'sample_reads_mysiam',
                    'sample_reads_mysiam_part']:

        discover_snps(species, dataset, dataset + '_q30', min_baseq=30, min_mapq=30, min_dist=5, max_qtl=100000)
