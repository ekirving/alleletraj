#!/usr/bin/env python
# -*- coding: utf-8 -*-

from db_conn import db_conn
from populate_coverage import *

species = 'pig'

interval_ids = [
    2328,  # 12048 snps
    2017,  # 5637 snps
    2368,  # 4005 snps
    2554,  # 3052 snps
    263,   # 2005 snps
    1170,  # 1000 snps
    2383,  # 500 snps
    153    # 100 snps
]

dbc = db_conn()

intervals = dbc.get_records_sql("""
    SELECT *
      FROM intervals
     WHERE id IN ({list})
    """.format(list=",".join(interval_ids)))

# get all the valid samples
samples = dbc.get_records_sql(
    """SELECT s.*, GROUP_CONCAT(sf.path) paths
         FROM samples s
         JOIN sample_files sf
           ON sf.sample_id = s.id
        WHERE s.species = '{species}'
          AND s.valid = 1
     GROUP BY s.id""".format(species=species))

print "INFO: Processing {:,} intervals in {:,} {} samples".format(len(intervals), len(samples), species)

began = time()

# process the chromosomes without multi-threading
for interval in intervals.values():
    start = time()
    process_interval((interval, samples))
    print "(%s)." % timedelta(seconds=time() - start)
