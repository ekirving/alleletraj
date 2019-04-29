#!/usr/bin/env python
# -*- coding: utf-8 -*-

from multiprocessing import cpu_count

# QTLdb release number (changing this number will trigger a complete db rebuild)
QTLDB_RELEASE = 'rel37'

# the scientific name of each species (necessary for resolving Ensembl URLs)
BINOMIAL_NAME = {
    'cattle': 'Bos_taurus',
    'goat':   'Capra_hircus',
    'horse':  'Equus_caballus',
    'pig':    'Sus_scrofa',
}

# which reference assembly should we use for each species (changing the assembly will trigger a complete db rebuild)
REF_ASSEMBLY = {
    'cattle': 'UMD3.1',
    'goat':   'ARS1',
    'horse':  'EquCab2',
    'pig':    'Sscrofa10.2',
}

# the name of the chromosomes in each assembly
CHROMOSOMES = {
    # cattle assemblies
    'UMD3.1':     map(str, range(1, 30)) + ['X', 'MT'],
    'ARS-UCD1.2': map(str, range(1, 30)) + ['X', 'MT'],

    # goat assemblies
    'ARS1': map(str, range(1, 30)) + ['MT'],

    # horse assemblies
    'EquCab2':   map(str, range(1, 32)) + ['X'],  # NOTE: MT omitted intentionally, as it is absent from the BAM files
    'EquCab3.0': map(str, range(1, 32)) + ['X'],

    # pig assemblies
    'Sscrofa10.2': map(str, range(1, 19)) + ['X', 'Y', 'MT'],
    'Sscrofa11.1': map(str, range(1, 19)) + ['X', 'Y', 'MT'],
}

# TODO move into spreadsheet
BAM_FILES = {
    'horse': {
        'Arab_0237A': '/home/ludo/inbox/BAMs/modern/Arab_0237A_SAMN02439777.Horse_nuc_wY.realigned.bam',
        'Conn_0004A': '/home/ludo/inbox/BAMs/modern/Conn_0004A_CON15.Horse_nuc_wY.realigned.bam',
        'Duel_0238A': '/home/ludo/inbox/BAMs/modern/Duel_0238A_SAMN02422919.Horse_nuc_wY.realigned.bam',
        'Esom_0226A': '/home/ludo/inbox/BAMs/modern/Esom_0226A_Esomalicus.Horse_nuc_wY.realigned.bam',
        'Frie_0296A': '/home/ludo/inbox/BAMs/modern/Frie_0296A_SAMEA3951218.Horse_nuc_wY.realigned.bam',
        'FrMo_0065A': '/home/ludo/inbox/BAMs/modern/FrMo_0065A_FM1798.Horse_nuc_wY.realigned.bam',
        'Hano_0235A': '/home/ludo/inbox/BAMs/modern/Hano_0235A_SAMN02439779.Horse_nuc_wY.realigned.bam',
        'Heav_0269A': '/home/ludo/inbox/BAMs/modern/Heav_0269A_SAMN03955412.Horse_nuc_wY.realigned.bam',
        'Icel_0144A': '/home/ludo/inbox/BAMs/modern/Icel_0144A_P5782.Horse_nuc_wY.realigned.bam',
        'Icel_0247A': '/home/ludo/inbox/BAMs/modern/Icel_0247A_IS074.Horse_nuc_wY.realigned.bam',
        'Jeju_0275A': '/home/ludo/inbox/BAMs/modern/Jeju_0275A_SAMN01057172.Horse_nuc_wY.realigned.bam',
        'Marw_0239A': '/home/ludo/inbox/BAMs/modern/Marw_0239A_SRR1275408.Horse_nuc_wY.realigned.bam',
        'Mong_0153A': '/home/ludo/inbox/BAMs/modern/Mong_0153A_KB7754.Horse_nuc_wY.realigned.bam',
        'Mong_0215A': '/home/ludo/inbox/BAMs/modern/Mong_0215A_TG1111D2628.Horse_nuc_wY.realigned.bam',
        'Morg_0096A': '/home/ludo/inbox/BAMs/modern/Morg_0096A_EMS595.Horse_nuc_wY.realigned.bam',
        'Prze_0150A': '/home/ludo/inbox/BAMs/modern/Prze_0150A_KB3879.Horse_nuc_wY.realigned.bam',
        'Prze_0151A': '/home/ludo/inbox/BAMs/modern/Prze_0151A_KB7674.Horse_nuc_wY.realigned.bam',
        'Prze_0157A': '/home/ludo/inbox/BAMs/modern/Prze_0157A_SB293.Horse_nuc_wY.realigned.bam',
        'Prze_0158A': '/home/ludo/inbox/BAMs/modern/Prze_0158A_SB339.Horse_nuc_wY.realigned.bam',
        'Prze_0159A': '/home/ludo/inbox/BAMs/modern/Prze_0159A_SB4329.Horse_nuc_wY.realigned.bam',
        'Prze_0160A': '/home/ludo/inbox/BAMs/modern/Prze_0160A_SB533.Horse_nuc_wY.realigned.bam',
        'Quar_0073A': '/home/ludo/inbox/BAMs/modern/Quar_0073A_A2085.Horse_nuc_wY.realigned.bam',
        'Shet_0249A': '/home/ludo/inbox/BAMs/modern/Shet_0249A_SPH020.Horse_nuc_wY.realigned.bam',
        'Shet_0250A': '/home/ludo/inbox/BAMs/modern/Shet_0250A_SPH041.Horse_nuc_wY.realigned.bam',
        'Sorr_0236A': '/home/ludo/inbox/BAMs/modern/Sorr_0236A_SAMN02439778.Horse_nuc_wY.realigned.bam',
        'Stan_0081A': '/home/ludo/inbox/BAMs/modern/Stan_0081A_M5256.Horse_nuc_wY.realigned.bam',
        'Thor_0145A': '/home/ludo/inbox/BAMs/modern/Thor_0145A_Twilight.Horse_nuc_wY.realigned.bam',
        'Thor_0290A': '/home/ludo/inbox/BAMs/modern/Thor_0290A_SAMN01047706.Horse_nuc_wY.realigned.bam',
        'Yaku_0163A': '/home/ludo/inbox/BAMs/modern/Yaku_0163A_Yak1.Horse_nuc_wY.realigned.bam',
        'Yaku_0170A': '/home/ludo/inbox/BAMs/modern/Yaku_0170A_Yak8.Horse_nuc_wY.realigned.bam',
        'Yaku_0171A': '/home/ludo/inbox/BAMs/modern/Yaku_0171A_Yak9.Horse_nuc_wY.realigned.bam'
    }
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


# the per-site mutation rate
MUTATION_RATE = {
    # TODO cattle, pig and goat
    'horse': 7.242e-9  # https://www.pnas.org/content/111/52/E5661.full#sec-17
}

# the average generation time in years
GENERATION_TIME = {
    # TODO cattle and goat
    'pig': 5,  # see https://www.nature.com/articles/ng.3197
    'horse': 8  # see https://www.sciencedirect.com/science/article/pii/S0960982215010039
}

# TODO deprecated
# the reference effective population size
POPULATION_SIZE = {
    'pig': {  # see https://www.nature.com/articles/ng.3394
        'EUD': 20563,
        'EUW': 8497
    },
    'horse': {  # see https://www.sciencedirect.com/science/article/pii/S0960982215010039
        'DOM2': 16000,
        'DOM2WLD': 16000
    }
}

# how many CPU cores does this machine have
TOTAL_CORES = cpu_count()

# set how many cores a single working can use
CPU_CORES_LOW = int(TOTAL_CORES * 0.1)   # 10%
CPU_CORES_MED = int(TOTAL_CORES * 0.25)  # 25%
CPU_CORES_HIGH = int(TOTAL_CORES * 0.5)  # 50%
CPU_CORES_MAX = int(TOTAL_CORES * 0.9)   # 90%
