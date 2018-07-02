#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function

import unicodecsv as csv

from pipeline_utils import *


def generate_sample_input(population, modsnp_id):
    """
    Generate the 4-column sample input file for `selection` (Schraiber et al., 2016)

    See https://github.com/Schraiber/selection

    The sample input is a 4-column, white-space-separated file. Each line corresponds to a single sample, which can be
    just one individual or many individuals from approximately the same time period pooled together.

    For each sample, each column is
    1. the number of derived alleles
    2. the sample size (in haploid genomes)
    3. the most ancient end of the possible age of the sample (i.e. the oldest it could be)
    4. the most recent end of the possible age of the sample (i.e. the youngest it could be)
    """

    dbc = db_conn()

    print("INFO: Generating sample input file for SNP #{}".format(modsnp_id))

    # TODO filter for population
    samples = dbc.get_records_sql("""
        # get the ancient frequencies in each bin
        SELECT SUM(sr.base = ms.derived) AS derived_count,
               COUNT(sr.id) AS sample_size,
                -CAST(SUBSTRING_INDEX(sb.bin, ' - ',  1) AS SIGNED INTEGER) AS bin_high,
                -CAST(SUBSTRING_INDEX(sb.bin, ' - ', -1) AS SIGNED INTEGER) - 1 AS bin_low
          FROM modern_snps ms
          JOIN sample_reads sr
            ON sr.chrom = ms.chrom
           AND sr.site = ms.site
           AND sr.called = 1
          JOIN samples s
            ON s.id = sr.sample_id
          JOIN sample_bins sb
            ON sb.sample_id = s.id
         WHERE ms.id = {modsnp_id}
      GROUP BY sb.bin

         UNION

        # add the modern frequency 
        SELECT derived_count, ancestral_count + derived_count, 0, 0
          FROM modern_snps ms
         WHERE ms.id = {modsnp_id}

      ORDER BY bin_high""".format(modsnp_id=modsnp_id), key=None)

    # write the sample input file
    with open("selection/{}-{}-modsnp_{}.input".format(SPECIES, population, modsnp_id), "wb") as tsv_file:

        fields = ['derived_count', 'sample_size', 'bin_high', 'bin_low']
        writer = csv.DictWriter(tsv_file, fieldnames=fields, delimiter='\t')

        # write the data to disk
        for sample in samples:
            writer.writerow(sample)


def run_selection(population, modsnp_id):
    """
    Run `selection` for the given SNP
    """

    # compose the input and output file paths
    input_file = "selection/{}-{}-modsnp_{}.input".format(SPECIES, population, modsnp_id)
    output_file = "selection/{}-{}-modsnp_{}".format(SPECIES, population, modsnp_id)
    pop_history = "data/selection/{}-{}-constant.pop".format(SPECIES, population)

    # get the generation time and Ne
    gen_time = GENERATION_TIME[SPECIES]
    pop_size = POPULATION_SIZE[SPECIES][population]

    run_cmd(['sr',
             '-D', input_file,        # path to data input file
             '-P', pop_history,       # path to population size history file
             '-o', output_file,       # output file prefix
             '-a',                    # flag to infer allele age
             '-G', gen_time,          # generation time
             '-N', pop_size,          # reference population size
             '-n', MCMC_CYCLES,       # number of MCMC cycles to run
             '-f', MCMC_PRINT,        # frequency of printing output to the screen
             '-s', MCMC_SAMPLE_FREQ,  # frequency of sampling from the posterior
             '-F', MCMC_FRACTION,     # fraction of the allele frequency to update during a trajectory update move
             '-e', MCMC_RANDOM_SEED,  # random number seed
             ])


