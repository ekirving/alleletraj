#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function

import itertools
import multiprocessing as mp

from selection import *

population = 'DOM2'

dbc = db_conn()

# get the best QTL SNPs
# modsnps = dbc.get_records_sql("""
#     SELECT DISTINCT ms.id
#       FROM qtl_snps qs
#       JOIN modern_snps ms
#         ON ms.id = qs.modsnp_id
#      WHERE qs.num_reads > 30
#        AND ms.neutral IS NULL""").keys()

# hand picked horse SNPs of interest
modsnps = [5333105,   # White markings
           6892292,   # White markings
           16832731,  # Temperament / heterospecific interaction trait
           7733000,   # Altitude adaptation
           6101266,   # Racing performance
           11908332,  # Sperm concentration
           13569570,  # Sperm count
           12492837,  # Semen volume
           ]

if MULTI_THREADED:
    # process the SNPs with multi-threading to make this faster
    pool = mp.Pool(MAX_CPU_CORES)
    pool.map(model_selection, itertools.izip(itertools.repeat(population), modsnps))
else:
    # process the SNPs without multi-threading
    for modsnp_id in modsnps:
        model_selection((population, modsnp_id))

# plot singly to avoid exceeding memory capacity of server
for modsnp_id in modsnps:
    # plot the allele trajectory
    plot_selection(population, modsnp_id)