#!/usr/bin/env python
# -*- coding: utf-8 -*-

import itertools
import pathos.multiprocessing as mp

from alleletraj.db.conn import Database
from alleletraj.analysis.selection import selection_fetch_neutral_snps

# how many CPU cores does this machine have
TOTAL_CORES = mp.cpu_count()

SPECIES = 'cattle'
POPULATION = 'Dom'

dbc = Database(SPECIES)

# get the modsnp id for every GWAS hit
records = dbc.get_records_sql("""
    SELECT DISTINCT ms.id, ifnull(ms.mispolar, 0) mispolar
      FROM qtls q
      JOIN modern_snps ms
        ON ms.chrom = q.chrom
       AND ms.site = q.site
     WHERE q.associationType = 'Association'
       AND q.valid = 1
       """, key=None)

modsnps = []
mispolar = []

for record in records:
    modsnps.append(record['id'])
    mispolar.append(0)

    if record['mispolar']:
        modsnps.append(record['id'])
        mispolar.append(1)


# compute the model likelihoods
pool = mp.ProcessingPool(TOTAL_CORES)
pool.map(selection_fetch_neutral_snps, itertools.repeat(SPECIES), itertools.repeat(POPULATION), modsnps, mispolar)

