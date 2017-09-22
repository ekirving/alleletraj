#!/usr/bin/env python
# -*- coding: utf-8 -*-

from db_conn import db_conn

from populate_qtls import populate_qtls
from populate_samples import populate_samples
from populate_coverage import *
from discover_snps import discover_snps

# populate_qtls()
# populate_samples()
# populate_intervals()
# populate_coverage()

# discover_snps('sample_reads_q20', min_baseq=20, min_mapq=20, min_dist=3, max_qtl=100000)
discover_snps('sample_reads_q30', min_baseq=30, min_mapq=30, min_dist=5, max_qtl=100000)
