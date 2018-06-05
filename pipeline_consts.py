#!/usr/bin/env python
# -*- coding: utf-8 -*-

from collections import OrderedDict
from socket import gethostname
from multiprocessing import cpu_count

# enforce max interval size of 1 Gb
MAX_INTERVAL_SIZE = int(1e6)

# offset to use for the QTL window (+/- 50 Kb)
QTL_WINDOW = 50000

# sizes of each chrom in the given assemblies
CHROM_SIZE = {

    # UMD_3.1.1
    'cattle': {'1': 158337067, '2': 137060424, '3': 121430405, '4': 120829699, '5': 121191424, '6': 119458736,
               '7': 112638659, '8': 113384836, '9': 105708250, '10': 104305016, '11': 107310763, '12': 91163125,
               '13': 84240350, '14': 84648390, '15': 85296676, '16': 81724687, '17': 75158596, '18': 66004023,
               '19': 64057457, '20': 72042655, '21': 71599096, '22': 61435874, '23': 52530062, '24': 62714930,
               '25': 42904170, '26': 51681464, '27': 45407902, '28': 46312546, '29': 51505224, 'X': 148823899},

    # Sscrofa10.2
    'pig':    {'1': 315321322, '2': 162569375, '3': 144787322, '4': 143465943, '5': 111506441, '6': 157765593,
               '7': 134764511, '8': 148491826, '9': 153670197, '10': 79102373, '11': 87690581, '12': 63588571,
               '13': 218635234, '14': 153851969, '15': 157681621, '16': 86898991, '17': 69701581, '18': 61220071,
               'X': 144288218, 'Y': 1637716},

    # EquCab2.0
    'horse':  {'1': 185838109, '2': 120857687, '3': 119479920, '4': 108569075, '5': 99680356, '6': 84719076,
               '7': 98542428, '8': 94057673, '9': 83561422, '10': 83980604, '11': 61308211, '12': 33091231,
               '13': 42578167, '14': 93904894, '15': 91571448, '16': 87365405, '17': 80757907, '18': 82527541,
               '19': 59975221, '20': 64166202, '21': 57723302, '22': 49946797, '23': 55726280, '24': 46749900,
               '25': 39536964, '26': 41866177, '27': 39960074, '28': 46177339, '29': 33672925, '30': 30062385,
               '31': 24984650, 'X': 124114077},

    # CHIR_1.0
    'goat':   {'1': 155011307, '2': 135415751, '3': 116796116, '4': 115961478, '5': 111055201, '6': 114334461,
               '7': 106547263, '8': 111020524, '9': 90293942, '10': 99198151, '11': 105305070, '12': 82535142,
               '13': 80625018, '14': 92306894, '15': 78986926, '16': 77678508, '17': 71877645, '18': 61067880,
               '19': 62130014, '20': 71279863, '21': 66773250, '22': 57956300, '23': 49403180, '24': 61756751,
               '25': 41496684, '26': 50169583, '27': 44118849, '28': 43231948, '29': 48376377, 'X': 121952644},
}

ENSEMBL_DATA = {

    # see ftp://ftp.ensembl.org/pub/release-89/variation/gvf/sus_scrofa/
    #     ftp://ftp.ensembl.org/pub/release-89/gtf/sus_scrofa/
    'pig': {'gtf': 'data/ensembl/release-89/gtf/sus_scrofa/Sus_scrofa.Sscrofa10.2.89.gtf.gz',
            'gvf': 'data/ensembl/release-89/gvf/sus_scrofa/Sus_scrofa.gvf.gz'},
}

SWEEP_DATA = {

    # see https://www.nature.com/articles/ng.3394
    'pig': {'loci': 'data/sweep/EUD_Sweep_p001_FINAL_cutoff_MERGED10kb.bed',
            'snps': 'data/sweep/EUD_Sweep_p001_FINAL_cutoff.bed'}
}

SNP_CHIP_DATA = {

    # see http://bioinformatics.tecnoparco.org/SNPchimp
    'pig': 'data/SNPchimp/SNPchimp_pig.tsv.gz'
}

VERBOSE = True

# the minimum phred scaled genotype quality (30 = 99.9%)
MIN_BASE_QUAL = 30
MIN_MAP_QUAL = 30
MIN_GENO_QUAL = 30

# number of bases to soft clip
SOFT_CLIP_DIST = 3


# the arbitrary +/- age uncertainty for median age dates
MEDIAN_AGE_UNCERT = 100

BIN_WIDTH = 500
BIN_PERCENT = 0.5  # samples must overlap a bin by >= 50%

# Pigs_allTo20042016_shared / old Google sheet
# SHEET_ID = '154wbTEcAUPz4d5v7QLLSwIKTdK8AxXP5US-riCjt2og'
# SHEET_TABS = ['Europe and NE Pigs']  #, 'SE Asian Pigs']

GOOGLE_SHEET = {

    # Pig_Table_Final_05_03_18
    'pig': {
        'id':   '1IWCt8OtTz6USOmN5DO0jcYxZOLnnOVdstTGzRcBZolI',
        'tabs': ['Everything for the paper - updated'],
        'cols': OrderedDict([
                    ('Extract No.',       'accession'),
                    ('Total Reads',       'map_reads'),
                    ('% Mapped',          'map_prcnt'),
                    ('Age',               'age'),
                    ('Period',            'period'),
                    ('Location',          'location'),
                    ('Country',           'country'),
                    ('Wild/Dom Status',   'status'),
                    ('GMM Status',        'gmm_status'),
                    ('Group',             'group'),
                    ('Haplogroup',        'haplogroup'),
                    ('DNA',               'dna')
                ])
    }
}

# list of junk input to mask with NULL
SHEET_NA = ['n/a', 'NA', '-', '?', 'NULL', 'None', '...', '']

AGE_MAP = {

    'pig': {
        'id': '1bH5u_qDaFXJdTyybeahqgF7je17td0FdyOMs_tlECdA',
        'tabs': ['Age Map'],
        'cols': OrderedDict([
            ('Age',         'age'),
            ('Confident',   'confident'),
            ('Lower (BP)',  'lower'),
            ('Upper (BP)',  'upper'),
            ('Median (BP)', 'median'),
        ])
    }
}

C14_SHEET = {

    'pig': {
        'id': '1odoL9hQh87bLLe3yipbo-CKKXLvIgb5n_kfoqSALHi8',
        'tabs': ['All Dates'],
        'cols': OrderedDict([
            ('Extract_No',              'accession'),
            ('From Cal BP (Int Cal13)', 'lower'),
            ('To Cal BP',               'upper'),
        ])
    }
}


# list of permissible countries in Europe
EUROPE = ['Belgium', 'Bulgaria', 'Croatia', 'Czech Rep.', 'Denmark', 'England', 'Estonia', 'Faroes', 'France',
          'Germany', 'Greece', 'Hungary', 'Iceland', 'Italy', 'Macedonia (FYROM)', 'Moldova', 'Netherlands', 'Poland',
          'Portugal', 'Romania', 'Serbia', 'Slovakia', 'Spain', 'Sweden', 'Switzerland', 'Ukraine']

# the number of flanking SNPs (on either side) to include
QTL_FLANK_NUM_SNPS = 2

# the number of selective sweep SNPs to include
SWEEP_NUM_SNPS = 3

# the distance between sweep peaks
SWEEP_PEAK_WIDTH = 1000

# offset all genes by 100 Kb to preclude linkage with our 'neutral' SNPs
GENE_OFFSET = int(1e5)

# the Ensembl gene ID
MC1R_GENE_ID = 'ENSSSCG00000020924'

# the number of "neutral" SNPs to include in the ascertainament
NUM_NEUTRAL_SNPS = 50000

# minimum depth of coverage to call diploid genotypes
MIN_GENO_DEPTH = 10

# TODO what about the others species
# location of reference genome
REF_FILE = "fasta/Sus_scrofa.Sscrofa10.2.dna.toplevel.fa"

# minumum mapping quality (hard filtered)
HARD_MAPQ_CUTOFF = 20

# minimum base quality (hard filtered)
HARD_BASEQ_CUTOFF = 20

# should we use multi-threading to speed up record insertion
MULTI_THREADED = True if gethostname() != 'macbookpro.local' else False

# the minimum minor allele frequency of modern SNPs to include
MIN_MAF = 0.05

# no single worker should use more than 30% of the available cores
MAX_CPU_CORES = int(cpu_count() * 0.3)
