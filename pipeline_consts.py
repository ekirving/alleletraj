#!/usr/bin/env python
# -*- coding: utf-8 -*-

from collections import OrderedDict
from socket import gethostname
from multiprocessing import cpu_count

# TODO this needs refactoring to work with multiple species
# SPECIES = 'pig'
# POPULATION = 'EUD'

# OUTGROUP = 'SVSV01U01_Sverrucosus_rh'  # Sus verrucosus / Javan warty pig

SPECIES = 'horse'
POPULATION = 'DOM'

OUTGROUP = 'Esom_0226A'  # Equus africanus somaliensis / Somali wild ass

# location of reference genome
REF_FILE = {
    'pig':   "fasta/Sus_scrofa.Sscrofa10.2.dna.toplevel.fa",
    'horse': "fasta/Equus_caballus.EquCab2.dna.toplevel.fa",
}

# path to the folder containing the modern fasta files
FASTA_PATH = '/media/jbod/raid1-sdc1/laurent/full_run_results/Pig/modern/FASTA'

SAMPLES = {

    'pig': {

        # the 81 European domestic pigs
        'EUD': ['AS01F01_AnglerSattleschwein_rh', 'AS01F09_Angler_Sattelsw_rh', 'BB01M47_Bunte_Bentheimer_rh',
                'BK01F10_Berkshire_rh', 'BK01M20_Berkshire_rh', 'BS01F10_British_Saddle_rh',
                'BS01F35_British_Saddleback_rh', 'CA01F14_Calabrese_rh', 'CM01F17_Chato_Murciano_rh',
                'CM01F18_Chato_Murciano_rh', 'CS01F02_Cinta_Senese_rh', 'CT01F13_Cassertana_rh', 'CT01M12_Cassertana_rh',
                'DU22M01_Duroc_rh', 'DU22M02_Duroc_rh', 'DU22M03_Duroc_rh', 'DU23M01_Duroc_rh', 'DU23M02_Duroc_rh',
                'DU23M03_Duroc_rh', 'DU23M04_Duroc_rh', 'GO01F04_Gl_Old_Spots_rh', 'GO01F23_GloucesterOldSpot_rh',
                'HA20U01_Hampshire_rh', 'HA20U02_Hampshire_rh', 'HA20U04_Hampshire_rh', 'HA20U06_Hampshire_rh',
                'LB01F49_Large_Black_rh', 'LE01F25_Leicoma_rh', 'LR21M03_rh', 'LR24F01_rh', 'LR24F08_rh', 'LR24M17_rh',
                'LR24M18_rh', 'LR24M19_rh', 'LR24M20_rh', 'LR24M21_rh', 'LR24M22_rh', 'LR30F02_rh', 'LR30F03_rh',
                'LR30F04_Landrace_rh', 'LS01F04_Linderodsvin_rh', 'LW22F01_rh', 'LW22F02_rh', 'LW22F03_rh', 'LW22F04_rh',
                'LW22F06_rh', 'LW22F07_rh', 'LW22F08_LargeWhite_rh', 'LW22F09_LargeWhite_rh', 'LW22M04_rh', 'LW36F01_rh',
                'LW36F02_rh', 'LW36F03_rh', 'LW36F04_rh', 'LW36F05_rh', 'LW36F06_rh', 'LW37M01_rh', 'LW38MF02_rh',
                'LW39M01_rh', 'LW39M02_rh', 'LW39M03_rh', 'LW39M04_rh', 'LW39M05_rh', 'LW39M07_rh', 'LW39M08_rh',
                'MA01F18_Mangalica_rh', 'MA01F20_Mangalica_rh', 'MW01F29_Middle_White_rh', 'MW01F33_Middle_White_rh',
                'NI01U07_Negro_Iberico_rh', 'NS01F05_Nera_Siciliana_rh', 'PI21F02_rh', 'PI21F06_rh', 'PI21F07_Pietrain_rh',
                'PI21F08_Pietrain_rh', 'PI21F09_Pietrain_rh', 'PI21M17_rh', 'PI21M20_rh', 'RE01F51_Retinto_rh',
                'TA01F19_Tamworth_rh', 'TA01M06_Tamworth_rh'],

        # the 22 Asian domestic pigs
        'ASD': ['JI01U08_Jinhua_rh', 'JI01U10_Jinhua_rh', 'JQ01U02_Jiangquhai_rh', 'JQ01U03_Jiangquahai_rh',
                'JQ01U08_Jiangquahai_rh', 'LSP01U16_LepingSpotted_rh', 'LSP01U18_LepingSpotted_rh', 'MS20M03_Meishan_rh',
                'MS20M05_Meishan_rh', 'MS20U10_Meishan_rh', 'MS20U11_Meishan_rh', 'MS20U13_Meishan_rh',
                'MS21M01_Meishan_rh', 'MS21M05_Meishan_rh', 'MS21M07_Meishan_rh', 'MS21M08_Meishan_rh',
                'MS21M14_Meishan_rh', 'WS01U03_WannanSpotted_rh', 'WS01U13_WannanSpotted_rh', 'XI01U03_rh',
                'XI01U04_Xiang_rh', 'ZA01U02_Zang_rh'],

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

# possible values include:
FASTA_MAP = {
    'R': ['A', 'G'],
    'Y': ['C', 'T'],
    'K': ['G', 'T'],
    'M': ['A', 'C'],
    'S': ['C', 'G'],
    'W': ['A', 'T'],
}

# enforce max interval size of 1 Gb
MAX_INTERVAL_SIZE = int(1e6)

# offset to use for the QTL window (+/- 50 Kb)
QTL_WINDOW = 50000

# sizes of each chrom in the given assemblies
CHROM_SIZE = {

    # Sscrofa10.2
    'pig': OrderedDict([
        ('1',  315321322), ('2',  162569375), ('3',  144787322), ('4', 143465943), ('5', 111506441), ('6', 157765593),
        ('7',  134764511), ('8',  148491826), ('9',  153670197), ('10', 79102373), ('11', 87690581), ('12', 63588571),
        ('13', 218635234), ('14', 153851969), ('15', 157681621), ('16', 86898991), ('17', 69701581), ('18', 61220071),
        ('X',  144288218), ('Y',    1637716)
    ]),

    # EquCab2.0
    'horse': OrderedDict([
        ('1', 185838109), ('2', 120857687), ('3', 119479920), ('4', 108569075), ('5',  99680356), ('6',  84719076),
        ('7',  98542428), ('8',  94057673), ('9',  83561422), ('10', 83980604), ('11', 61308211), ('12', 33091231),
        ('13', 42578167), ('14', 93904894), ('15', 91571448), ('16', 87365405), ('17', 80757907), ('18', 82527541),
        ('19', 59975221), ('20', 64166202), ('21', 57723302), ('22', 49946797), ('23', 55726280), ('24', 46749900),
        ('25', 39536964), ('26', 41866177), ('27', 39960074), ('28', 46177339), ('29', 33672925), ('30', 30062385),
        ('31', 24984650), ('X', 124114077)
    ]),

    # CHIR_1.0
    'goat': OrderedDict([
        ('1', 155011307), ('2', 135415751), ('3', 116796116), ('4', 115961478), ('5', 111055201), ('6', 114334461),
        ('7', 106547263), ('8', 111020524), ('9', 90293942), ('10', 99198151), ('11', 105305070), ('12', 82535142),
        ('13', 80625018), ('14', 92306894), ('15', 78986926), ('16', 77678508), ('17', 71877645), ('18', 61067880),
        ('19', 62130014), ('20', 71279863), ('21', 66773250), ('22', 57956300), ('23', 49403180), ('24', 61756751),
        ('25', 41496684), ('26', 50169583), ('27', 44118849), ('28', 43231948), ('29', 48376377), ('X', 121952644)
    ]),

    # UMD_3.1.1
    'cattle': OrderedDict([
        ('1', 158337067), ('2', 137060424), ('3', 121430405), ('4',  120829699), ('5',  121191424), ('6', 119458736),
        ('7', 112638659), ('8', 113384836), ('9', 105708250), ('10', 104305016), ('11', 107310763), ('12', 91163125),
        ('13', 84240350), ('14', 84648390), ('15', 85296676), ('16',  81724687), ('17',  75158596), ('18', 66004023),
        ('19', 64057457), ('20', 72042655), ('21', 71599096), ('22',  61435874), ('23',  52530062), ('24', 62714930),
        ('25', 42904170), ('26', 51681464), ('27', 45407902), ('28',  46312546), ('29',  51505224), ('X', 148823899)
    ]),

}

ENSEMBL_DATA = {

    # see ftp://ftp.ensembl.org/pub/release-89/gtf/sus_scrofa/
    #     ftp://ftp.ensembl.org/pub/release-89/variation/gvf/sus_scrofa/
    'pig': {'gtf': 'data/ensembl/release-89/gtf/sus_scrofa/Sus_scrofa.Sscrofa10.2.89.gtf.gz',
            'gvf': 'data/ensembl/release-89/gvf/sus_scrofa/Sus_scrofa.gvf.gz'},

    # see ftp://ftp.ensembl.org/pub/release-92/gtf/equus_caballus/
    #     ftp://ftp.ensembl.org/pub/release-92/variation/gvf/equus_caballus/
    'horse': {'gtf': 'data/ensembl/release-92/gtf/equus_caballus/Equus_caballus.EquCab2.92.gtf.gz',
              'gvf': 'data/ensembl/release-92/gvf/equus_caballus/equus_caballus.gvf.gz'},

    # see ftp://ftp.ensembl.org/pub/release-92/gtf/capra_hircus/
    #     ftp://ftp.ensembl.org/pub/release-92/variation/gvf/capra_hircus/
    'goat':  {'gtf': 'data/ensembl/release-92/gtf/capra_hircus/Capra_hircus.ARS1.92.gtf.gz',
              'gvf': 'data/ensembl/release-92/gvf/capra_hircus/capra_hircus.gvf.gz'},
}

SWEEP_DATA = {

    # see https://www.nature.com/articles/ng.3394
    'pig': {'loci': 'data/sweep/EUD_Sweep_p001_FINAL_cutoff_MERGED10kb.bed',
            'snps': 'data/sweep/EUD_Sweep_p001_FINAL_cutoff.bed'}
}

SNP_CHIP_DATA = {

    # see http://bioinformatics.tecnoparco.org/SNPchimp
    'pig':   'data/SNPchimp/SNPchimp_pig.tsv.gz',
    'horse': 'data/SNPchimp/SNPchimp_horse.tsv.gz',
    'goat':  'data/SNPchimp/SNPchimp_goat.tsv.gz'
}

# the Ensembl gene ID
MC1R_GENE_ID = {
    'pig':   'ENSSSCG00000020924',
    'horse': 'ENSECAG00000000900'
}

VERBOSE = False

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
    },

    # HorseSelection_LO4EIP-TRANSFERED
    'horse': {
        'id':   '1BMvIwYj-d8t3mpf67rzabrEvDoB8hBZbyS6XfGwcwUU',
        'tabs': ['Ancient'],
        'cols': OrderedDict([
                    ('Name',     'accession'),
                    ('Status',   'status'),
                    ('path',     'path'),
                    ('Age BP',   'age'),
                    ('Age',      'period'),
                    ('Site',     'location'),
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

RADIOCARBON_SHEET = {

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
QTL_FLANK_NUM_SNPS = 3

# the number of selective sweep SNPs to include
SWEEP_NUM_SNPS = 5

# the distance between sweep peaks
SWEEP_PEAK_WIDTH = 1000

# offset all genes by 100 Kb to preclude linkage with our 'neutral' SNPs
GENE_OFFSET = 100000

# the number of "neutral" SNPs to include in the ascertainment
NUM_NEUTRAL_SNPS = 60000

# the number of "ancestral" SNPs to include in the ascertainment
NUM_ANCESTRAL_SNPS = 30000

# minimum distance from an INDEL
INDEL_BUFFER = 10

# minimum depth of coverage to call diploid genotypes
MIN_GENO_DEPTH = 10

# minumum mapping quality (hard filtered)
HARD_MAPQ_CUTOFF = 20

# minimum base quality (hard filtered)
HARD_BASEQ_CUTOFF = 20

# should we use multi-threading to speed up record insertion
MULTI_THREADED = True if gethostname() != 'macbookpro.local' else False

# no single worker should use more than 50% of the available cores
MAX_CPU_CORES = int(cpu_count() * 0.5)

# the minimum derived allele frequency of modern SNPs to include
MIN_DAF = 0.05

# the minimum phred scaled genotype quality (30 = 99.9%)
MIN_BASE_QUAL = 30
MIN_MAP_QUAL = 30
MIN_GENO_QUAL = 30

# number of bases to soft clip
SOFT_CLIP_DIST = 3

# the maximum number of rows to insert in a single operation
MAX_INSERT_SIZE = 50000

# the maximum number of conditions in a single query
MAX_QUERY_SIZE = 5000
