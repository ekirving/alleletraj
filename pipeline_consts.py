#!/usr/bin/env python
# -*- coding: utf-8 -*-

from collections import OrderedDict
from multiprocessing import cpu_count
from socket import gethostname

# TODO move module specific constants into those modules

# QTLdb settings
QTLDB_RELEASE = 'rel37'

# the scientific name of each species
BINOMIAL_NAME = {
    'cattle': 'Bos_taurus',
    'goat':   'Capra_hircus',
    'horse':  'Equus_caballus',
    'pig':    'Sus_scrofa',
}

# which reference assembly should we use for each species
REF_ASSEMBLY = {
    'cattle': 'UMD3.1',
    'goat':   'ARS1',
    'horse':  'EquCab2',
    'pig':    'Sscrofa10.2',
}

# TODO extract these from the fasta file index
# sizes of each chrom in the given assemblies
CHROM_SIZE = {

    'Sscrofa10.2': OrderedDict([
        ('1',  315321322), ('2',  162569375), ('3',  144787322), ('4', 143465943), ('5', 111506441), ('6', 157765593),
        ('7',  134764511), ('8',  148491826), ('9',  153670197), ('10', 79102373), ('11', 87690581), ('12', 63588571),
        ('13', 218635234), ('14', 153851969), ('15', 157681621), ('16', 86898991), ('17', 69701581), ('18', 61220071),
        ('X',  144288218), ('Y',    1637716)
    ]),

    'EquCab2': OrderedDict([
        ('1', 185838109), ('2', 120857687), ('3', 119479920), ('4', 108569075), ('5',  99680356), ('6',  84719076),
        ('7',  98542428), ('8',  94057673), ('9',  83561422), ('10', 83980604), ('11', 61308211), ('12', 33091231),
        ('13', 42578167), ('14', 93904894), ('15', 91571448), ('16', 87365405), ('17', 80757907), ('18', 82527541),
        ('19', 59975221), ('20', 64166202), ('21', 57723302), ('22', 49946797), ('23', 55726280), ('24', 46749900),
        ('25', 39536964), ('26', 41866177), ('27', 39960074), ('28', 46177339), ('29', 33672925), ('30', 30062385),
        ('31', 24984650), ('X', 124114077)
    ]),

    # TODO wrong!
    'CHIR_1.0': OrderedDict([
        ('1', 155011307), ('2', 135415751), ('3', 116796116), ('4', 115961478), ('5', 111055201), ('6', 114334461),
        ('7', 106547263), ('8', 111020524), ('9', 90293942), ('10', 99198151), ('11', 105305070), ('12', 82535142),
        ('13', 80625018), ('14', 92306894), ('15', 78986926), ('16', 77678508), ('17', 71877645), ('18', 61067880),
        ('19', 62130014), ('20', 71279863), ('21', 66773250), ('22', 57956300), ('23', 49403180), ('24', 61756751),
        ('25', 41496684), ('26', 50169583), ('27', 44118849), ('28', 43231948), ('29', 48376377), ('X', 121952644)
    ]),

    'UMD3.1': OrderedDict([
        ('1', 158337067), ('2', 137060424), ('3', 121430405), ('4',  120829699), ('5',  121191424), ('6', 119458736),
        ('7', 112638659), ('8', 113384836), ('9', 105708250), ('10', 104305016), ('11', 107310763), ('12', 91163125),
        ('13', 84240350), ('14', 84648390), ('15', 85296676), ('16',  81724687), ('17',  75158596), ('18', 66004023),
        ('19', 64057457), ('20', 72042655), ('21', 71599096), ('22',  61435874), ('23',  52530062), ('24', 62714930),
        ('25', 42904170), ('26', 51681464), ('27', 45407902), ('28',  46312546), ('29',  51505224), ('X', 148823899)
    ]),

}

OUTGROUP = {
    # 'cattle': '',  # TODO add me
    # 'goat':   '',  # TODO add me
    'horse':  'Esom_0226A',                # Equus africanus somaliensis / Somali wild ass
    'pig':    'SVSV01U01_Sverrucosus_rh',  # Sus verrucosus / Javan warty pig
}

# TODO move into spreadsheet
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


# the per-site mutation rate
MUTATION_RATE = {
    'horse': 7.242e-9  # https://www.pnas.org/content/111/52/E5661.full#sec-17
}

# the average generation time in years
GENERATION_TIME = {
    'pig': 5,   # see https://www.nature.com/articles/ng.3197
    'horse': 8  # see https://www.sciencedirect.com/science/article/pii/S0960982215010039
}

# the reference effective population size
POPULATION_SIZE = {
    'pig': {    # see https://www.nature.com/articles/ng.3394
        'EUD': 20563,
        'EUW': 8497
    },
    'horse': {  # see https://www.sciencedirect.com/science/article/pii/S0960982215010039
        'DOM2': 16000,
        'DOM2WLD': 16000
    }
}

# the minimum derived allele frequency of modern SNPs to include
MIN_DAF = 0.05

# the minimum phred scaled genotype quality (30 = 99.9%)
MIN_GENO_QUAL = 30
