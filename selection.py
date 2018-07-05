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


def generate_mc1r_snp_input(population):
    """
    Custom code for handling the PCR data for modsnp #71891
    """

    dbc = db_conn()

    print("INFO: Generating sample input file for SNP #71891")

    samples = dbc.get_records_sql("""
        # get the ancient frequencies in each bin
        SELECT SUM(CASE s.mc1r_snp 
                       WHEN 'A/A' THEN 2 
                       WHEN 'G/A' THEN 1 
                       ELSE 0 
                   END) AS derived_count,
               count(s.id)*2 AS sample_size,
               -CAST(SUBSTRING_INDEX(sb.bin, ' - ',  1) AS SIGNED INTEGER) AS bin_high,
               -CAST(SUBSTRING_INDEX(sb.bin, ' - ', -1) AS SIGNED INTEGER) - 1 AS bin_low
          FROM samples s
          JOIN sample_bins sb
            ON sb.sample_id = s.id
          WHERE s.valid = 1 
            AND s.status IN ('Domestic')
            AND s.mc1r_snp IS NOT NULL 
       GROUP BY sb.bin

          UNION

          # add the modern frequency 
          SELECT derived_count, ancestral_count + derived_count, 0, 0
            FROM modern_snps ms
           WHERE ms.id = 71891

        ORDER BY bin_high""", key=None)

    # write the sample input file
    with open("selection/{}-{}-modsnp_71891.input".format(SPECIES, population), "wb") as tsv_file:

        fields = ['derived_count', 'sample_size', 'bin_high', 'bin_low']
        writer = csv.DictWriter(tsv_file, fieldnames=fields, delimiter='\t')

        # write the data to disk
        for sample in samples:
            writer.writerow(sample)


def run_selection(population, modsnp_id, mcmc_cycles, mcmc_freq):
    """
    Run `selection` for the given SNP.
    """

    begin = time()

    # compose the input and output file paths
    input_file = "selection/{}-{}-modsnp_{}.input".format(SPECIES, population, modsnp_id)
    output_prefix = "selection/{}-{}-modsnp_{}_c{}-f{}".format(SPECIES, population, modsnp_id, mcmc_cycles, mcmc_freq)

    # get the generation time, Ne and population history
    gen_time = GENERATION_TIME[SPECIES]
    pop_size = POPULATION_SIZE[SPECIES][population]
    pop_hist = POPULATION_HISTORY[SPECIES][population]

    print("INFO: Started selection for {} SNP #{}".format(population, modsnp_id))

    log = run_cmd(['sr',
                   '-D', input_file,        # path to data input file
                   '-P', pop_hist,          # path to population size history file
                   '-o', output_prefix,     # output file prefix
                   '-a',                    # flag to infer allele age
                   '-G', gen_time,          # generation time
                   '-N', pop_size,          # reference population size
                   '-n', mcmc_cycles,       # number of MCMC cycles to run
                   '-f', MCMC_PRINT,        # frequency of printing output to the screen
                   '-s', mcmc_freq,         # frequency of sampling from the posterior
                   '-F', MCMC_FRACTION,     # fraction of the allele frequency to update during a trajectory update move
                   '-e', MCMC_RANDOM_SEED,  # random number seed
                   ])

    # save the log file
    with open(output_prefix + '.log', 'w') as fout:
        fout.write(log)

    # TODO measure ESS and enforce threshold
    # https://www.rdocumentation.org/packages/LaplacesDemon/versions/16.1.0/topics/ESS
    # https://cran.r-project.org/web/packages/coda/index.html

    print("INFO: Finished selection for {} SNP #{} ({})".format(population, modsnp_id, timedelta(seconds=time() - begin)))


def plot_selection(population, modsnp_id, mcmc_cycles, mcmc_freq):
    """
    Plot the allele trajectory.

    WARNING: very slow and memory costly!
    """

    begin = time()

    # compose the input and output file paths
    input_file = "selection/{}-{}-modsnp_{}.input".format(SPECIES, population, modsnp_id)
    output_prefix = "selection/{}-{}-modsnp_{}_c{}-f{}".format(SPECIES, population, modsnp_id, mcmc_cycles, mcmc_freq)
    gen_time = GENERATION_TIME[SPECIES]
    pop_size = POPULATION_SIZE[SPECIES][population]

    print("INFO: Started plotting for {} SNP #{}".format(population, modsnp_id))

    # burn in by 20%
    burn_in = int(mcmc_cycles / mcmc_freq * 0.2)

    # plot the allele trajectory
    run_cmd(['Rscript', 'rscript/plot-selection.R', input_file, output_prefix, gen_time, pop_size, burn_in])

    # TODO plot the strength of selection (a1 and a2)

    print("INFO: Finished plotting for {} SNP #{} ({})".format(population, modsnp_id, timedelta(seconds=time() - begin)))


def model_selection(args):
    """
    Model the allele trajectory of the given SNP.
    """

    # extract the nested tuple of arguments (an artifact of using izip to pass args to mp.Pool)
    (population, modsnp_id, (mcmc_cycles, mcmc_freq)) = args

    if modsnp_id == 71891:
        generate_mc1r_snp_input(population)
    else:
        # convert the SNP data into the input format for `selection`
        generate_sample_input(population, modsnp_id)

    # run `selection` for the given SNP
    run_selection(population, modsnp_id, mcmc_cycles, mcmc_freq)

    # plot the allele trajectory
    plot_selection(population, modsnp_id, mcmc_cycles, mcmc_freq)
