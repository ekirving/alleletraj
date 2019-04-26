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

# the name of the outgroup sample
OUTGROUP = {
    # 'cattle': '',  # TODO add me
    # 'goat':   '',  # TODO add me
    'horse': 'Esom_0226A',  # Equus africanus somaliensis / Somali wild ass
    'pig': 'SVSV01U01_Sverrucosus_rh',  # Sus verrucosus / Javan warty pig
}

# TODO move into spreadsheet
SAMPLES = {

    # https://www.nature.com/articles/ng.3394
    'pig': {

        # the 81 European domestic pigs
        'EUD': ['AS01F01_AnglerSattleschwein_rh', 'AS01F09_Angler_Sattelsw_rh', 'BB01M47_Bunte_Bentheimer_rh',
                'BK01F10_Berkshire_rh', 'BK01M20_Berkshire_rh', 'BS01F10_British_Saddle_rh',
                'BS01F35_British_Saddleback_rh', 'CA01F14_Calabrese_rh', 'CM01F17_Chato_Murciano_rh',
                'CM01F18_Chato_Murciano_rh', 'CS01F02_Cinta_Senese_rh', 'CT01F13_Cassertana_rh',
                'CT01M12_Cassertana_rh', 'DU22M01_Duroc_rh', 'DU22M02_Duroc_rh', 'DU22M03_Duroc_rh', 'DU23M01_Duroc_rh',
                'DU23M02_Duroc_rh', 'DU23M03_Duroc_rh', 'DU23M04_Duroc_rh', 'GO01F04_Gl_Old_Spots_rh',
                'GO01F23_GloucesterOldSpot_rh', 'HA20U01_Hampshire_rh', 'HA20U02_Hampshire_rh', 'HA20U04_Hampshire_rh',
                'HA20U06_Hampshire_rh', 'LB01F49_Large_Black_rh', 'LE01F25_Leicoma_rh', 'LR21M03_rh', 'LR24F01_rh',
                'LR24F08_rh', 'LR24M17_rh', 'LR24M18_rh', 'LR24M19_rh', 'LR24M20_rh', 'LR24M21_rh', 'LR24M22_rh',
                'LR30F02_rh', 'LR30F03_rh', 'LR30F04_Landrace_rh', 'LS01F04_Linderodsvin_rh', 'LW22F01_rh',
                'LW22F02_rh', 'LW22F03_rh', 'LW22F04_rh', 'LW22F06_rh', 'LW22F07_rh', 'LW22F08_LargeWhite_rh',
                'LW22F09_LargeWhite_rh', 'LW22M04_rh', 'LW36F01_rh', 'LW36F02_rh', 'LW36F03_rh', 'LW36F04_rh',
                'LW36F05_rh', 'LW36F06_rh', 'LW37M01_rh', 'LW38MF02_rh', 'LW39M01_rh', 'LW39M02_rh', 'LW39M03_rh',
                'LW39M04_rh', 'LW39M05_rh', 'LW39M07_rh', 'LW39M08_rh', 'MA01F18_Mangalica_rh', 'MA01F20_Mangalica_rh',
                'MW01F29_Middle_White_rh', 'MW01F33_Middle_White_rh', 'NI01U07_Negro_Iberico_rh',
                'NS01F05_Nera_Siciliana_rh', 'PI21F02_rh', 'PI21F06_rh', 'PI21F07_Pietrain_rh', 'PI21F08_Pietrain_rh',
                'PI21F09_Pietrain_rh', 'PI21M17_rh', 'PI21M20_rh', 'RE01F51_Retinto_rh', 'TA01F19_Tamworth_rh',
                'TA01M06_Tamworth_rh'],

        # the 22 Asian domestic pigs
        'ASD': ['JI01U08_Jinhua_rh', 'JI01U10_Jinhua_rh', 'JQ01U02_Jiangquhai_rh', 'JQ01U03_Jiangquahai_rh',
                'JQ01U08_Jiangquahai_rh', 'LSP01U16_LepingSpotted_rh', 'LSP01U18_LepingSpotted_rh',
                'MS20M03_Meishan_rh', 'MS20M05_Meishan_rh', 'MS20U10_Meishan_rh', 'MS20U11_Meishan_rh',
                'MS20U13_Meishan_rh', 'MS21M01_Meishan_rh', 'MS21M05_Meishan_rh', 'MS21M07_Meishan_rh',
                'MS21M08_Meishan_rh', 'MS21M14_Meishan_rh', 'WS01U03_WannanSpotted_rh', 'WS01U13_WannanSpotted_rh',
                'XI01U03_rh', 'XI01U04_Xiang_rh', 'ZA01U02_Zang_rh'],

        # the 2 Sumatran scrofa (for ascertaining ancient alleles)
        'SUM': ['INDO22_Sumatra_rh', 'INDO33_Sumatra_rh'],
    },

    'horse': {

        # the 6 Przewalski horses
        'DOM': ['Prze_0150A', 'Prze_0151A', 'Prze_0157A', 'Prze_0158A', 'Prze_0159A', 'Prze_0160A'],

        # the 24 Domestic horses
        'DOM2': ['Arab_0237A', 'Conn_0004A', 'Duel_0238A', 'FrMo_0065A', 'Frie_0296A', 'Hano_0235A', 'Heav_0269A',
                 'Icel_0144A', 'Icel_0247A', 'Jeju_0275A', 'Marw_0239A', 'Mong_0153A', 'Mong_0215A', 'Morg_0096A',
                 'Quar_0073A', 'Shet_0249A', 'Shet_0250A', 'Sorr_0236A', 'Stan_0081A', 'Thor_0145A', 'Thor_0290A',
                 'Yaku_0163A', 'Yaku_0170A', 'Yaku_0171A'],
    }
}

# TODO move into spreadsheet
SAMPLE_SEX = {

    'horse': {
        'Arab_0237A': 'M',
        'Conn_0004A': 'M',
        'Duel_0238A': 'F',
        'Esom_0226A': 'F',
        'Frie_0296A': 'M',
        'FrMo_0065A': 'M',
        'Hano_0235A': 'M',
        'Heav_0269A': 'M',
        'Icel_0144A': 'M',
        'Icel_0247A': 'M',
        'Jeju_0275A': 'M',
        'Marw_0239A': 'M',
        'Mong_0153A': 'M',
        'Mong_0215A': 'M',
        'Morg_0096A': 'F',
        'Prze_0150A': 'F',
        'Prze_0151A': 'M',
        'Prze_0157A': 'M',
        'Prze_0158A': 'F',
        'Prze_0159A': 'M',
        'Prze_0160A': 'M',
        'Quar_0073A': 'M',
        'Shet_0249A': 'F',
        'Shet_0250A': 'F',
        'Sorr_0236A': 'M',
        'Stan_0081A': 'M',
        'Thor_0145A': 'F',
        'Thor_0290A': 'M',
        'Yaku_0163A': 'M',
        'Yaku_0170A': 'M',
        'Yaku_0171A': 'M'
    }
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

# TODO move into a spreadsheet
# TODO testing data only
SRA_ACCESSIONS = {
    'pig': {
        'AS01F01_AnglerSattleschwein_rh': ['ERR977190'],
        'AS01F09_Angler_Sattelsw_rh': ['ERR977191', 'ERR977192', 'ERR977193', 'ERR977194']
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
    'horse': 7.242e-9  # https://www.pnas.org/content/111/52/E5661.full#sec-17
}

# the average generation time in years
GENERATION_TIME = {
    'pig': 5,  # see https://www.nature.com/articles/ng.3197
    'horse': 8  # see https://www.sciencedirect.com/science/article/pii/S0960982215010039
}

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
