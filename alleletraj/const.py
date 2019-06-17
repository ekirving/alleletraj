#!/usr/bin/env python
# -*- coding: utf-8 -*-

# standard modules
from multiprocessing import cpu_count

# QTLdb release number (changing this number will trigger a complete db rebuild)
QTLDB_RELEASE = 'rel38'

# the scientific name of each species (necessary for resolving Ensembl URLs)
BINOMIAL_NAME = {
    'cattle': 'Bos_taurus',
    'goat':   'Capra_hircus',
    'horse':  'Equus_caballus',
    'pig':    'Sus_scrofa',
}

# which reference assembly should we use for each species (changing the assembly will trigger a complete db rebuild)
REF_ASSEMBLY = {
    'cattle': 'UMD3.1',       # TODO replace with ARS-UCD1.2
    'goat':   'ARS1',
    'horse':  'EquCab2',      # TODO replace with EquCab3.0
    'pig':    'Sscrofa10.2',  # TODO replace with Sscrofa11.1
}

# the name of the chromosomes in each assembly (NOTE: MT omitted intentionally)
CHROMOSOMES = {
    # cattle assemblies
    'UMD3.1':     map(str, range(1, 30)) + ['X'],
    'ARS-UCD1.2': map(str, range(1, 30)) + ['X'],

    # goat assemblies
    'ARS1': map(str, range(1, 30)),

    # horse assemblies
    'EquCab2':   map(str, range(1, 32)) + ['X'],
    'EquCab3.0': map(str, range(1, 32)) + ['X'],

    # pig assemblies
    'Sscrofa10.2': map(str, range(1, 19)) + ['X', 'Y'],
    'Sscrofa11.1': map(str, range(1, 19)) + ['X', 'Y'],
}

# genomic regions of selective sweeps ascertained in other papers
SWEEP_DATA = {
    # 'cattle': {}, # TODO add other species
    # 'goat': {},
    # 'pig': {},

    # see https://www.nature.com/articles/ng.3394
    'pig': {
        'loci': 'data/sweep/EUD_Sweep_p001_FINAL_cutoff_MERGED10kb.bed',
        'snps': 'data/sweep/EUD_Sweep_p001_FINAL_cutoff.bed'
    }
}


# the per-site mutation rate mutation rate per generation
MUTATION_RATE = {
    'cattle': 1.10e-8,    # 5 years / see https://www.frontiersin.org/articles/10.3389/fgene.2018.00728/full
    # 'cattle': 1.26e-8,  # 6 years / see https://www.nature.com/articles/s41467-018-04737-0
    # 'cattle': 2.20e-9 / 5,        # see https://www.nature.com/articles/s41559-018-0562-y
    #                     actually citing https://www.pnas.org/content/99/2/803

    'goat':   2.5e-8,    # https://www.nature.com/articles/s41467-018-03206-y

    'horse':  7.242e-9,  # https://www.pnas.org/content/111/52/E5661.full#sec-17
                         # https://www.sciencedirect.com/science/article/pii/S0960982215010039

    'pig':    2.5e-8,    # https://www.nature.com/articles/ng.3197
}

# the average generation time in years
GENERATION_TIME = {
    'cattle': 5,  # https://www.frontiersin.org/articles/10.3389/fgene.2018.00728/full
    'goat':   2,  # https://www.nature.com/articles/s41467-018-03206-y
    'horse':  8,  # https://www.pnas.org/content/111/52/E5661.full#sec-17
    'pig':    5,  # https://www.nature.com/articles/ng.3197
}

# TODO deprecated
# the reference effective population size
# POPULATION_SIZE = {
#     'pig': {  # see https://www.nature.com/articles/ng.3394
#         'EUD': 20563,
#         'EUW': 8497
#     },
#     'horse': {  # see https://www.sciencedirect.com/science/article/pii/S0960982215010039
#         'DOM2': 16000,
#         'DOM2WLD': 16000
#     }
# }



# how many CPU cores does this machine have
TOTAL_CORES = cpu_count()

# set how many cores a single working can use
CPU_CORES_LOW = max(int(TOTAL_CORES * 0.1), 1)
CPU_CORES_MED = int(TOTAL_CORES * 0.25)
CPU_CORES_HIGH = int(TOTAL_CORES * 0.5)
CPU_CORES_MAX = int(TOTAL_CORES * 0.9)

# what level should Treemix group by, pops OR samples
GROUP_BY_POPS = 'grp-pops'
GROUP_BY_SMPL = 'grp-smpl'

# the name of the outgroup population
OUTGROUP_POP = 'OUT'
