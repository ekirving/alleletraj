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

# load_ensembl_genes()
# load_ensembl_variants()

# load_snpchip_variants()

discover_modern_snps()

# populate_qtls()
# compute_qtl_windows()
populate_sweeps()
# populate_mc1r_locus()
# populate_neutral_loci()

# populate_intervals()
# populate_interval_snps()

# populate_samples()
# populate_coverage()

# populate_qtl_snps()
# mark_neutral_snps()

# discover_snps()
# analyse_qtls()

# perform_ascertainment()

