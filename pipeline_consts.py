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
        ('30', 30062385), ('31', 24984650)  # TODO restore when done testing
        # ('1', 185838109), ('2', 120857687), ('3', 119479920), ('4', 108569075), ('5',  99680356), ('6',  84719076),
        # ('7',  98542428), ('8',  94057673), ('9',  83561422), ('10', 83980604), ('11', 61308211), ('12', 33091231),
        # ('13', 42578167), ('14', 93904894), ('15', 91571448), ('16', 87365405), ('17', 80757907), ('18', 82527541),
        # ('19', 59975221), ('20', 64166202), ('21', 57723302), ('22', 49946797), ('23', 55726280), ('24', 46749900),
        # ('25', 39536964), ('26', 41866177), ('27', 39960074), ('28', 46177339), ('29', 33672925), ('30', 30062385),
        # ('31', 24984650), ('X', 124114077)
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
        'A_Ch40_CGG_1_016168': '/home/ludo/inbox/BAMs/ancient/A_Ch40_CGG_1_016168.Horse_nuc_wY.realigned.r.t.bam',
        'AM115_CGG_1_018579': '/home/ludo/inbox/BAMs/ancient/AM115_CGG_1_018579.Horse_nuc_wY.realigned.r.t.m.bam',
        'AM181_CGG_1_018580': '/home/ludo/inbox/BAMs/ancient/AM181_CGG_1_018580.Horse_nuc_wY.realigned.r.t.m.bam',
        'Arab_0237A': '/home/ludo/inbox/BAMs/modern/Arab_0237A_SAMN02439777.Horse_nuc_wY.realigned.bam',
        'ARUS_0222A': '/home/ludo/inbox/BAMs/ancient/ARUS_0222A_CGG101397.Horse_nuc_wY.realigned.r.t.m.bam',
        'ARUS_0223A': '/home/ludo/inbox/BAMs/ancient/ARUS_0223A_Batagai.Horse_nuc_wY.realigned.r.t.m.bam',
        'ARUS_0224A': '/home/ludo/inbox/BAMs/ancient/ARUS_0224A_CGG10022.Horse_nuc_wY.realigned.r.t.m.bam',
        'ARUS_0225A': '/home/ludo/inbox/BAMs/ancient/ARUS_0225A_CGG10023.Horse_nuc_wY.realigned.r.t.m.bam',
        'Arz15_CGG_1_017084': '/home/ludo/inbox/BAMs/ancient/Arz15_CGG_1_017084_U_i4_TGACCA.r.t.3p1.bam',
        'Arz17_CGG_1_017086': '/home/ludo/inbox/BAMs/ancient/Arz17_CGG_1_017086_U_i5_ACAGTG.r.t.3p1.bam',
        'Arz3_CGG_1_017089': '/home/ludo/inbox/BAMs/ancient/Arz3_CGG_1_017089_i9_GATCAG_U.r.t.5p1.3p1.bam',
        'B_Ch24_CGG_1_016169': '/home/ludo/inbox/BAMs/ancient/B_Ch24_CGG_1_016169.Horse_nuc_wY.realigned.r.t.bam',
        'BAPSKA_unregistered': '/home/ludo/inbox/BAMs/ancient/BAPSKA_unregistered_i63_CAACTG_U.r.t.3p1.bam',
        'Belgheis_TrBWBX116_CGG_1_019521': '/home/ludo/inbox/BAMs/ancient/Belgheis_TrBWBX116_CGG_1_019521.Horse_nuc_wY.realigned.r.t.m.bam',
        'Botai1_CGG_1_018173_Extract1': '/home/ludo/inbox/BAMs/ancient/Botai1_CGG_1_018173.Horse_nuc_wY.realigned.r.t.m.bam',
        'Botai2_CGG_1_018174': '/home/ludo/inbox/BAMs/ancient/Botai2_CGG_1_018174.Horse_nuc_wY.realigned.r.t.m.bam',
        'Botai3_CGG_1_018175': '/home/ludo/inbox/BAMs/ancient/Botai3_CGG_1_018175.Horse_nuc_wY.realigned.r.t.m.bam',
        'Botai4_CGG_1_018176': '/home/ludo/inbox/BAMs/ancient/Botai4_CGG_1_018176.Horse_nuc_wY.realigned.r.t.m.bam',
        'Botai5_CGG_1_018177': '/home/ludo/inbox/BAMs/ancient/Botai5_CGG_1_018177.Horse_nuc_wY.realigned.r.t.m.bam',
        'Botai6_CGG_1_018178': '/home/ludo/inbox/BAMs/ancient/Botai6_CGG_1_018178.Horse_nuc_wY.realigned.r.t.m.bam',
        'Bru4_CGG_1_018376': '/home/ludo/inbox/BAMs/ancient/Bru4_CGG_1_018376.Horse_nuc_wY.realigned.r.t.s.m.bam',
        'Cap102_CGG_1_016984': '/home/ludo/inbox/BAMs/ancient/Cap102_CGG_1_016984.Horse_nuc_wY.realigned.r.t.s.m.bam',
        'CdY2_CGG_1_018391': '/home/ludo/inbox/BAMs/ancient/CdY2_CGG_1_018391.Horse_nuc_wY.realigned.r.t.m.bam',
        'Conn_0004A': '/home/ludo/inbox/BAMs/modern/Conn_0004A_CON15.Horse_nuc_wY.realigned.bam',
        'D_Ch47_CGG_1_016171': '/home/ludo/inbox/BAMs/ancient/D_Ch47_CGG_1_016171.Horse_nuc_wY.realigned.r.t.bam',
        'DaAn01_CGG_1_018337': '/home/ludo/inbox/BAMs/ancient/DaAn01_CGG_1_018337_i28_TCTCGC_U.r.t.3p1.bam',
        'Duel_0238A': '/home/ludo/inbox/BAMs/modern/Duel_0238A_SAMN02422919.Horse_nuc_wY.realigned.bam',
        'Duk2_CGG_1_018386': '/home/ludo/inbox/BAMs/ancient/Duk2_CGG_1_018386.Horse_nuc_wY.realigned.r.t.m.bam',
        'E_Ch25_CGG_1_016172': '/home/ludo/inbox/BAMs/ancient/E_Ch25_CGG_1_016172.Horse_nuc_wY.realigned.r.t.t.s.bam',
        'Earb5_CGG_1_018497': '/home/ludo/inbox/BAMs/ancient/Earb5_CGG_1_018497_i80_TTTTGG_U.r.t.3p1.bam',
        'Earb6_CGG_1_018495': '/home/ludo/inbox/BAMs/ancient/Earb6_CGG_1_018495.Horse_nuc_wY.realigned.r.t.s.m.3p1.bam',
        'Esom_0226A': '/home/ludo/inbox/BAMs/modern/Esom_0226A_Esomalicus.Horse_nuc_wY.realigned.bam',
        'F_Ch26_CGG_1_016173': '/home/ludo/inbox/BAMs/ancient/F_Ch26_CGG_1_016173.Horse_nuc_wY.realigned.r.t.bam',
        'Fen4_CGG_1_018396': '/home/ludo/inbox/BAMs/ancient/Fen4_CGG_1_018396.Horse_nuc_wY.realigned.r.t.s.5p2.3p2.bam',
        'Fr1_CGG_1_018151': '/home/ludo/inbox/BAMs/ancient/Fr1_CGG_1_018151.Horse_nuc_wY.realigned.r.t.m.bam',
        'Frie_0296A': '/home/ludo/inbox/BAMs/modern/Frie_0296A_SAMEA3951218.Horse_nuc_wY.realigned.bam',
        'FrMo_0065A': '/home/ludo/inbox/BAMs/modern/FrMo_0065A_FM1798.Horse_nuc_wY.realigned.bam',
        'G_Ch82_CGG_1_016174': '/home/ludo/inbox/BAMs/ancient/G_Ch82_CGG_1_016174.Horse_nuc_wY.realigned.r.t.bam',
        'Gar3_CGG_1_018389': '/home/ludo/inbox/BAMs/ancient/Gar3_CGG_1_018389.Horse_nuc_wY.realigned.r.t.m.bam',
        'Georgia2_CGG_1_020266': '/home/ludo/inbox/BAMs/ancient/Georgia2_CGG_1_020266_i13_CTATCA_U.r.t.s.3p2.bam',
        'GEP13_CGG_1_018049_5': '/home/ludo/inbox/BAMs/ancient/GEP13_CGG_1_018049.Horse_nuc_wY.realigned.r.t.m.bam',
        'GEP14_CGG_1_018050': '/home/ludo/inbox/BAMs/ancient/GEP14_CGG_1_018050.Horse_nuc_wY.realigned.r.t.m.bam',
        'GEP21_CGG_1_018057_Digestion2': '/home/ludo/inbox/BAMs/ancient/GEP21.Horse_nuc_wY.realigned.r.t.m.bam',
        'GVA1_CGG_1_019259': '/home/ludo/inbox/BAMs/ancient/GVA1_CGG_1_019259_i52_U_AAACAA.r.t.bam',
        'GVA111_CGG_1_019369': '/home/ludo/inbox/BAMs/ancient/GVA111_CGG_1_019369_i74_GTGGGG_U.r.t.5p2.3p1.bam',
        'GVA112_CGG_1_019370': '/home/ludo/inbox/BAMs/ancient/GVA112_CGG_1_019370_i11_GGCTAC_U.r.t.bam',
        'GVA115_CGG_1_019373': '/home/ludo/inbox/BAMs/ancient/GVA115_CGG_1_019373_i13_CTATCA_U.r.t.5p1.3p1.bam',
        'GVA122_CGG_1_019380': '/home/ludo/inbox/BAMs/ancient/GVA122_CGG_1_019380.Horse_nuc_wY.realigned.r.t.m.bam',
        'GVA123_CGG_1_019381': '/home/ludo/inbox/BAMs/ancient/GVA123_CGG_1_019381.Horse_nuc_wY.realigned.r.t.m.bam',
        'GVA124_CGG_1_019382': '/home/ludo/inbox/BAMs/ancient/GVA124_CGG_1_019382.Horse_nuc_wY.realigned.r.t.m.bam',
        'GVA126_CGG_1_019384': '/home/ludo/inbox/BAMs/ancient/GVA126_CGG_1_019384.Horse_nuc_wY.realigned.r.t.m.bam',
        'GVA129_CGG_1_019387': '/home/ludo/inbox/BAMs/ancient/GVA129_CGG_1_019387_i16_GTGTAT_U.r.t.3p1.bam',
        'GVA133_CGG_1_019707': '/home/ludo/inbox/BAMs/ancient/GVA133_CGG_1_019707.Horse_nuc_wY.realigned.r.t.m.1.3p.bam',
        'GVA135_CGG_1_019709': '/home/ludo/inbox/BAMs/ancient/GVA135_CGG_1_019709.Horse_nuc_wY.realigned.r.t.m.5p2.31p.bam',
        'GVA140_CGG_1_019714': '/home/ludo/inbox/BAMs/ancient/GVA140_CGG_1_019714.Horse_nuc_wY.realigned.r.t.m.1.3p.bam',
        'GVA191_CGG_1_019765': '/home/ludo/inbox/BAMs/ancient/GVA191_CGG_1_019765.Horse_nuc_wY.realigned.r.t.m.bam',
        'GVA199_CGG_1_019773': '/home/ludo/inbox/BAMs/ancient/GVA199_CGG_1_019773.Horse_nuc_wY.realigned.r.t.m.5p2.3p2.bam',
        'GVA201_CGG_1_019775': '/home/ludo/inbox/BAMs/ancient/GVA201_CGG_1_019775.Horse_nuc_wY.realigned.r.t.m.bam',
        'GVA237_CGG_1_019877': '/home/ludo/inbox/BAMs/ancient/GVA237_CGG_1_019877_i55_AGATCG_U.r.t.3p1.bam',
        'GVA242_CGG_1_019816': '/home/ludo/inbox/BAMs/ancient/GVA242_CGG_1_019816.Horse_nuc_wY.realigned.r.t.m.bam',
        'GVA26_CGG_1_019284': '/home/ludo/inbox/BAMs/ancient/GVA26_CGG_1_019284.Horse_nuc_wY.realigned.r.t.m.bam',
        'GVA28_CGG_1_019286': '/home/ludo/inbox/BAMs/ancient/GVA28_CGG_1_019286_i2_CGATGT_U.r.t.bam',
        'GVA307_CGG_1_019857': '/home/ludo/inbox/BAMs/ancient/GVA307_CGG_1_019857.Horse_nuc_wY.realigned.r.t.m.bam',
        'GVA308_CGG_1_019858': '/home/ludo/inbox/BAMs/ancient/GVA308_CGG_1_019858.Horse_nuc_wY.realigned.r.t.m.bam',
        'GVA309_CGG_1_019859': '/home/ludo/inbox/BAMs/ancient/GVA309_CGG_1_019859.Horse_nuc_wY.realigned.r.t.m.bam',
        'GVA310_CGG_1_019860': '/home/ludo/inbox/BAMs/ancient/GVA310_CGG_1_019860_i56_AGGGGA_U.r.t.bam',
        'GVA311_CGG_1_019861': '/home/ludo/inbox/BAMs/ancient/GVA311_CGG_1_019861.Horse_nuc_wY.realigned.r.t.m.5p2.3p3.bam',
        'GVA321_CGG_1_019871': '/home/ludo/inbox/BAMs/ancient/GVA321_CGG_1_019871_i01_U_ATCACG.r.t.3p1.bam',
        'GVA36_CGG_1_019294': '/home/ludo/inbox/BAMs/ancient/GVA36_CGG_1_019294.Horse_nuc_wY.realigned.r.t.m.1.3p.1.5p.1.3p.bam',
        'GVA375_CGG_1_019925': '/home/ludo/inbox/BAMs/ancient/GVA375_CGG_1_019925.Horse_nuc_wY.realigned.r.t.m.bam',
        'GVA4_CGG_1_019262': '/home/ludo/inbox/BAMs/ancient/GVA4_CGG_1_019262.Horse_nuc_wY.realigned.r.t.m.5p2.3p1.bam',
        'GVA43_CGG_1_019301': '/home/ludo/inbox/BAMs/ancient/GVA43_CGG_1_019301.Horse_nuc_wY.realigned.r.t.m.5p2.3p1.s.bam',
        'GVA47_CGG_1_019305': '/home/ludo/inbox/BAMs/ancient/GVA47_CGG_1_019305_i56_U_AGGGGA.r.t.3p1.bam',
        'GVA48_CGG_1_019306': '/home/ludo/inbox/BAMs/ancient/GVA48_CGG_1_019306_i57_U_ATACCT.r.t.bam',
        'GVA53_CGG_1_019311': '/home/ludo/inbox/BAMs/ancient/GVA53_CGG_1_019311_i71_GCTAGC_U.r.t.5p1.bam',
        'GVA56_CGG_1_019314': '/home/ludo/inbox/BAMs/ancient/GVA56_CGG_1_019314.Horse_nuc_wY.realigned.r.t.m.bam',
        'GVA60_CGG_1_019318': '/home/ludo/inbox/BAMs/ancient/GVA60_CGG_1_019318_i72_GGGCCG_U.r.t.5p1.3p2.bam',
        'GVA75_CGG_1_019333': '/home/ludo/inbox/BAMs/ancient/GVA75_CGG_1_019333_i75_TAGTAA_U.r.t.5p1.bam',
        'GVA81_CGG_1_019339': '/home/ludo/inbox/BAMs/ancient/GVA81_CGG_1_019339.Horse_nuc_wY.realigned.r.t.m.5p2.bam',
        'GVA9_CGG_1_019267': '/home/ludo/inbox/BAMs/ancient/GVA9_CGG_1_019267_i5_ACAGTG_U.r.t.s.bam',
        'H_Ch78_CGG_1_016175': '/home/ludo/inbox/BAMs/ancient/H_Ch78_CGG_1_016175.Horse_nuc_wY.realigned.r.t.t.bam',
        'Hano_0235A': '/home/ludo/inbox/BAMs/modern/Hano_0235A_SAMN02439779.Horse_nuc_wY.realigned.bam',
        'Hasanlu1140_CGG_1_019998': '/home/ludo/inbox/BAMs/ancient/Hasanlu1140_CGG_1_019998_i4_TGACCA_U.r.t.5p2.3p2.bam',
        'Hasanlu2327_CGG_1_019995': '/home/ludo/inbox/BAMs/ancient/Hasanlu2327_CGG_1_019995_i56_AGGGGA_U.r.t.5p2.3p1.bam',
        'Hasanlu2405_CGG_1_019992': '/home/ludo/inbox/BAMs/ancient/Hasanlu2405_CGG_1_019992.Horse_nuc_wY.realigned.r.t.m.bam',
        'Hasanlu2529_CGG_1_019988': '/home/ludo/inbox/BAMs/ancient/Hasanlu2529_CGG_1_019988_i57_ATACCT_U.r.t.3p1.bam',
        'Hasanlu2689_CGG_1_019996': '/home/ludo/inbox/BAMs/ancient/Hasanlu2689_CGG_1_019996_i58_ATGAGC_U.r.t.5p1.3p2.bam',
        'Hasanlu3394_CGG_1_019997': '/home/ludo/inbox/BAMs/ancient/Hasanlu3394_CGG_1_019997_i69_GATGCA_U.r.t.3p1.bam',
        'Hasanlu3398_CGG_1_019986': '/home/ludo/inbox/BAMs/ancient/Hasanlu3398_CGG_1_019986_i47_CATAGA_U.r.t.5p2.3p2.bam',
        'Hasanlu3461_CGG_1_020003': '/home/ludo/inbox/BAMs/ancient/Hasanlu3461_CGG_1_020003_i37_ACTGCC_U.r.t.5p1.3p1.bam',
        'Hasanlu368_CGG_1_019994': '/home/ludo/inbox/BAMs/ancient/Hasanlu368_CGG_1_019994_i52_AAACAA_U.r.5p2.3p1.bam',
        'Haunstetten_CGG_1_017139': '/home/ludo/inbox/BAMs/ancient/Haunstetten_CGG_1_017139.Horse_nuc_wY.realigned.r.t.m.bam',
        'Heav_0269A': '/home/ludo/inbox/BAMs/modern/Heav_0269A_SAMN03955412.Horse_nuc_wY.realigned.bam',
        'I_Ch118_CGG_1_016176': '/home/ludo/inbox/BAMs/ancient/I_Ch118_CGG_1_016176.Horse_nuc_wY.realigned.r.t.bam',
        'I-K2_Arz1_CGG_1_017079': '/home/ludo/inbox/BAMs/ancient/I-K2_Arz1_CGG_1_017079.Horse_nuc_wY.realigned.r.t.t.bam',
        'I-K3_Arz2_CGG_1_017088': '/home/ludo/inbox/BAMs/ancient/I-K3_Arz2_CGG_1_017088.Horse_nuc_wY.realigned.r.t.t.s.bam',
        'Iberia254_CGG_1_020032': '/home/ludo/inbox/BAMs/ancient/Iberia254_CGG_1_020032_i69_GATGCA_U.r.t.5p2.3p3.bam',
        'Icel_0144A': '/home/ludo/inbox/BAMs/modern/Icel_0144A_P5782.Horse_nuc_wY.realigned.bam',
        'Icel_0247A': '/home/ludo/inbox/BAMs/modern/Icel_0247A_IS074.Horse_nuc_wY.realigned.bam',
        'Issyk1_CGG_1_018577': '/home/ludo/inbox/BAMs/ancient/Issyk1_CGG_1_018577.Horse_nuc_wY.realigned.r.t.m.bam',
        'Jeju_0275A': '/home/ludo/inbox/BAMs/modern/Jeju_0275A_SAMN01057172.Horse_nuc_wY.realigned.bam',
        'JG160_CGG_1_019246': '/home/ludo/inbox/BAMs/ancient/JG160_CGG_1_019246.Horse_nuc_wY.realigned.r.t.m.bam',
        'K_Ch137_CGG_1_016177': '/home/ludo/inbox/BAMs/ancient/K_Ch137_CGG_1_016177.Horse_nuc_wY.realigned.r.t.bam',
        'Kha2_t1_CGG_1_018909': '/home/ludo/inbox/BAMs/ancient/Kha2_t1_CGG_1_018909.Horse_nuc_wY.realigned.r.t.m.bam',
        'KSH4_CGG_1_017098': '/home/ludo/inbox/BAMs/ancient/KSH4_CGG_1_017098.Horse_nuc_wY.realigned.r.t.m.5p1.3p1.bam',
        'KSH5_CGG_1_017099': '/home/ludo/inbox/BAMs/ancient/KSH5_CGG_1_017099_i54_ACCATC_U.r.t.s.3p2.bam',
        'KYRH10_CGG_1_018031': '/home/ludo/inbox/BAMs/ancient/KYRH10_CGG_1_018031.Horse_nuc_wY.realigned.r.t.m.bam',
        'KYRH8_CGG_1_018029': '/home/ludo/inbox/BAMs/ancient/KYRH8_CGG_1_018029.Horse_nuc_wY.realigned.r.t.m.5p2.3p1.bam',
        'L_Ch102_CGG_1_016178': '/home/ludo/inbox/BAMs/ancient/L_Ch102_CGG_1_016178.Horse_nuc_wY.realigned.r.t.bam',
        'LIS2_Tr19': '/home/ludo/inbox/BAMs/ancient/LIS2_Tr19_i8_ACTTGA_U.r.t.5p1.bam',
        'LOBOT_8_CGG_1_020179': '/home/ludo/inbox/BAMs/ancient/LOBOT_8_CGG_1_020179.Horse_nuc_wY.realigned.r.t.m.bam',
        'LOBOT_A_CGG_1_020181': '/home/ludo/inbox/BAMs/ancient/LOBOT_A_CGG_1_020181.Horse_nuc_wY.realigned.r.t.m.1.3p.bam',
        'LOBOT_B_CGG_1_020182': '/home/ludo/inbox/BAMs/ancient/LOBOT_B_CGG_1_020182.Horse_nuc_wY.realigned.r.t.m.3p1.bam',
        'LOBOT_C_CGG_1_020183': '/home/ludo/inbox/BAMs/ancient/LOBOT_C_CGG_1_020183.Horse_nuc_wY.realigned.r.t.m.bam',
        'LOBOT_D1_CGG_1_020204': '/home/ludo/inbox/BAMs/ancient/LOBOT_D1_CGG_1_020204.Horse_nuc_wY.realigned.r.t.m.bam',
        'LOBOT_D2_CGG_1_020205': '/home/ludo/inbox/BAMs/ancient/LOBOT_D2_CGG_1_020205_i21_ACATAC_U.r.t.bam',
        'LOBOT_D4_CGG_1_020207': '/home/ludo/inbox/BAMs/ancient/LOBOT_D4_CGG_1_020207.Horse_nuc_wY.realigned.r.t.m.bam',
        'LOBOT_D5_CGG_1_020208': '/home/ludo/inbox/BAMs/ancient/LOBOT_D5_CGG_1_020208.Horse_nuc_wY.realigned.r.t.m.bam',
        'LOBOT_D6_CGG_1_020209': '/home/ludo/inbox/BAMs/ancient/LOBOT_D6_CGG_1_020209.Horse_nuc_wY.realigned.r.t.m.bam',
        'LOBOT_E_CGG_1_020185': '/home/ludo/inbox/BAMs/ancient/LOBOT_E_CGG_1_020185_i27_TGCATA_U.r.t.3p2.bam',
        'LOBOT_F_CGG_1_020186': '/home/ludo/inbox/BAMs/ancient/LOBOT_F_CGG_1_020186.Horse_nuc_wY.realigned.r.t.m.bam',
        'LOBOT_G_CGG_1_020187': '/home/ludo/inbox/BAMs/ancient/LOBOT_G_CGG_1_020187.Horse_nuc_wY.realigned.r.t.m.bam',
        'LOBOT_I_CGG_1_020189': '/home/ludo/inbox/BAMs/ancient/LOBOT_I_CGG_1_020189.Horse_nuc_wY.realigned.r.t.m.bam',
        'LOBOT_K_CGG_1_020191': '/home/ludo/inbox/BAMs/ancient/LOBOT_K_CGG_1_020191.Horse_nuc_wY.realigned.r.t.m.bam',
        'LOBOT_L_CGG_1_020192': '/home/ludo/inbox/BAMs/ancient/LOBOT_L_CGG_1_020192.LOBOT_M_CGG_1_020193.Horse_nuc_wY.realigned.r.t.m.bam',
        'LOBOT_N_CGG_1_020194': '/home/ludo/inbox/BAMs/ancient/LOBOT_N_CGG_1_020194_i9_GATCAG_U.r.t.5p1.3p2.bam',
        'LOBOT_O_CGG_1_020195': '/home/ludo/inbox/BAMs/ancient/LOBOT_O_CGG_1_020195_i10_TAGCTT_U.r.t.5p2.3p2.bam',
        'LOBOT_P_CGG_1_020196': '/home/ludo/inbox/BAMs/ancient/LOBOT_P_CGG_1_020196.Horse_nuc_wY.realigned.r.t.m.bam',
        'LOBOT_Petrous_CGG_1_020210': '/home/ludo/inbox/BAMs/ancient/LOBOT_Petrous_CGG_1_020210.Horse_nuc_wY.realigned.r.t.m.bam',
        'LOBOT_R_CGG_1_020198': '/home/ludo/inbox/BAMs/ancient/LOBOT_R_CGG_1_020198.Horse_nuc_wY.realigned.r.t.m.bam',
        'LOBOT_T_CGG_1_020200': '/home/ludo/inbox/BAMs/ancient/LOBOT_T_CGG_1_020200_i47_U_CATAGA.r.t.bam',
        'M_Ch112_CGG_1_016179': '/home/ludo/inbox/BAMs/ancient/M_Ch112_CGG_1_016179.Horse_nuc_wY.realigned.r.t.bam',
        'Marvele01_CGG_1_019388': '/home/ludo/inbox/BAMs/ancient/Marvele01_CGG_1_019388.Horse_nuc_wY.realigned.r.t.m.bam',
        'Marvele02_CGG_1_019389': '/home/ludo/inbox/BAMs/ancient/Marvele02_CGG_1_019389_i13_CTATCA_U.r.t.3p2.bam',
        'Marvele05_CGG_1_019392': '/home/ludo/inbox/BAMs/ancient/Marvele05_CGG_1_019392_i21_ACATAC_U.r.t.3p2.bam',
        'Marvele16_CGG_1_019403': '/home/ludo/inbox/BAMs/ancient/Marvele16_CGG_1_019403_i28_TCTCGC_U.r.t.5p2.bam',
        'Marvele18_CGG_1_019405': '/home/ludo/inbox/BAMs/ancient/Marvele18_CGG_1_019405.Horse_nuc_wY.realigned.r.t.m.1.3p.bam',
        'Marvele21_CGG_1_019408': '/home/ludo/inbox/BAMs/ancient/Marvele21_CGG_1_019408.Horse_nuc_wY.realigned.r.t.m.5p2.3p1.bam',
        'Marvele22_CGG_1_019409': '/home/ludo/inbox/BAMs/ancient/Marvele22_CGG_1_019409_i46_CACGAA_U.r.t.bam',
        'Marvele27_CGG_1_019414': '/home/ludo/inbox/BAMs/ancient/Marvele27_CGG_1_019414_i47_CATAGA_U.r.t.bam',
        'Marvele32_CGG_1_019419': '/home/ludo/inbox/BAMs/ancient/Marvele32_CGG_1_019419.Horse_nuc_wY.realigned.r.t.m.bam',
        'Marw_0239A': '/home/ludo/inbox/BAMs/modern/Marw_0239A_SRR1275408.Horse_nuc_wY.realigned.bam',
        'Mic2_CGG_1_018388': '/home/ludo/inbox/BAMs/ancient/Mic2_CGG_1_018388_Enriched_i32_ACGCAT_U.r.t.s.bam',
        'Mon23_CGG_1_018059': '/home/ludo/inbox/BAMs/ancient/Mon23_CGG_1_018059.Horse_nuc_wY.realigned.r.t.m.1.3p.bam',
        'Mon24_CGG_1_018060_Extraction1': '/home/ludo/inbox/BAMs/ancient/Mon24_CGG_1_018060.Horse_nuc_wY.realigned.r.t.m.1.3p.bam',
        'Mon25_CGG_1_018061': '/home/ludo/inbox/BAMs/ancient/Mon25_CGG_1_018061.Horse_nuc_wY.realigned.r.t.m.1.3p.bam',
        'Mon26_CGG_1_018062': '/home/ludo/inbox/BAMs/ancient/Mon26_CGG_1_018062.Horse_nuc_wY.realigned.r.t.m.1.3p.bam',
        'Mon27_CGG_1_018063': '/home/ludo/inbox/BAMs/ancient/Mon27_CGG_1_018063.Horse_nuc_wY.realigned.r.t.m.1.3p.bam',
        'Mon28_CGG_1_018064': '/home/ludo/inbox/BAMs/ancient/Mon28_CGG_1_018064.Horse_nuc_wY.realigned.r.t.m.bam',
        'Mon37_CGG_1_018073': '/home/ludo/inbox/BAMs/ancient/Mon37_CGG_1_018073.Horse_nuc_wY.realigned.r.t.m.bam',
        'Mon39': '/home/ludo/inbox/BAMs/ancient/Mon39.Horse_nuc_wY.realigned.r.t.m.1.3p.bam',
        'Mon40_CGG_1_018076': '/home/ludo/inbox/BAMs/ancient/Mon40_CGG_1_018076.Horse_nuc_wY.realigned.r.t.s.m.3p1.bam',
        'Mon41_CGG_1_018077': '/home/ludo/inbox/BAMs/ancient/Mon41_CGG_1_018077.Horse_nuc_wY.realigned.r.t.m.bam',
        'Mon42_CGG_1_018078_Extraction1': '/home/ludo/inbox/BAMs/ancient/Mon42_CGG_1_018078.Horse_nuc_wY.realigned.r.t.m.5p2.3p1.bam',
        'Mon43_CGG_1_018079': '/home/ludo/inbox/BAMs/ancient/Mon43_CGG_1_018079.Horse_nuc_wY.realigned.r.t.m.5p2.3p1.bam',
        'Mon44_CGG_1_018080': '/home/ludo/inbox/BAMs/ancient/Mon44_CGG_1_018080.Horse_nuc_wY.realigned.r.t.m.bam',
        'Mon45_CGG_1_018081': '/home/ludo/inbox/BAMs/ancient/Mon45_CGG_1_018081.Horse_nuc_wY.realigned.r.t.m.5p2.3p1.bam',
        'Mon79_CGG_1_018115': '/home/ludo/inbox/BAMs/ancient/Mon79_CGG_1_018115.Horse_nuc_wY.realigned.r.t.m.bam',
        'Mon84_CGG_1_018120_Extraction2': '/home/ludo/inbox/BAMs/ancient/Mon84_CGG_1_018120.Horse_nuc_wY.realigned.r.t.m.bam',
        'Mon86_CGG_1_018122': '/home/ludo/inbox/BAMs/ancient/Mon86_CGG_1_018122.Horse_nuc_wY.realigned.r.t.m.bam',
        'Mon87_CGG_1_018123': '/home/ludo/inbox/BAMs/ancient/Mon87_CGG_1_018123.Horse_nuc_wY.realigned.r.t.m.bam',
        'Mon89': '/home/ludo/inbox/BAMs/ancient/Mon89_i28_TCTCGC_U.r.t.3p2.bam',
        'Mong_0153A': '/home/ludo/inbox/BAMs/modern/Mong_0153A_KB7754.Horse_nuc_wY.realigned.bam',
        'Mong_0215A': '/home/ludo/inbox/BAMs/modern/Mong_0215A_TG1111D2628.Horse_nuc_wY.realigned.bam',
        'Morg_0096A': '/home/ludo/inbox/BAMs/modern/Morg_0096A_EMS595.Horse_nuc_wY.realigned.bam',
        'MV178_CGG_1_020446': '/home/ludo/inbox/BAMs/ancient/MV178_CGG_1_020446.Horse_nuc_wY.realigned.r.t.m.bam',
        'Mzr1_CGG_1_018150': '/home/ludo/inbox/BAMs/ancient/Mzr1_CGG_1_018150.Horse_nuc_wY.realigned.r.t.m.bam',
        'NB_175_CGG_1_020499': '/home/ludo/inbox/BAMs/ancient/NB_175_CGG_1_020499_i26_TGTGAC_U.r.t.5p1.3p2.bam',
        'NB_63_CGG_1_020512': '/home/ludo/inbox/BAMs/ancient/NB_63_CGG_1_020512_i72_GGGCCG_U.r.t.5p1.3p2.bam',
        'NB_K9279_CGG_1_020505': '/home/ludo/inbox/BAMs/ancient/NB_K9279_CGG_1_020505_i65_CGGCAC_U.r.t.5p2.3p2.bam',
        'NB_P9261_CGG_1_020506': '/home/ludo/inbox/BAMs/ancient/NB_P9261_CGG_1_020506_i66_CTCGGT_U.r.t.5p1.3p2.bam',
        'NB_Ra_8_49_CGG_1_020509': '/home/ludo/inbox/BAMs/ancient/NB_Ra_8_49_CGG_1_020509_i69_GATGCA_U.r.t.5p1.3p2.bam',
        'NewBotai_13_CGG_1_017005': '/home/ludo/inbox/BAMs/ancient/NewBotai_13_CGG_1_017005.Horse_nuc_wY.realigned.r.t.s.m.5p2.3p1.bam',
        'NewBotai_18_CGG_1_017010': '/home/ludo/inbox/BAMs/ancient/NewBotai_18_CGG_1_017010.Horse_nuc_wY.realigned.r.t.s.m.bam',
        'NewBotai_31_CGG_1_017023': '/home/ludo/inbox/BAMs/ancient/NewBotai_31_CGG_1_017023.Horse_nuc_wY.realigned.r.t.s.m.5p3.3p2.bam',
        'NewBotai_4_CGG_1_016996': '/home/ludo/inbox/BAMs/ancient/NewBotai_4_CGG_1_016996.Horse_nuc_wY.realigned.r.t.m.bam',
        'NewBotai_9_CGG_1_017001': '/home/ludo/inbox/BAMs/ancient/NewBotai_9_CGG_1_017001.Horse_nuc_wY.realigned.r.t.s.m.bam',
        'NewBotai10': '/home/ludo/inbox/BAMs/ancient/NewBotai10.Horse_nuc_wY.realigned.r.t.s.m.bam',
        'NewBotai15_CGG_1_017007': '/home/ludo/inbox/BAMs/ancient/NewBotai15_CGG_1_017007.Horse_nuc_wY.realigned.r.t.s.m.5p2.3p2.bam',
        'NewBotai2_CGG_1_016994': '/home/ludo/inbox/BAMs/ancient/NewBotai2_CGG_1_016994.Horse_nuc_wY.realigned.r.t.s.m.bam',
        'NewBotai35_CGG_1_017027': '/home/ludo/inbox/BAMs/ancient/NewBotai35_CGG_1_017027.Horse_nuc_wY.realigned.r.t.m.2.5p.s.bam',
        'NewBotai44_CGG_1_017036': '/home/ludo/inbox/BAMs/ancient/NewBotai44_CGG_1_017036.Horse_nuc_wY.realigned.r.t.s.m.5p2.3p2.bam',
        'NewBotai45_CGG_1_017037': '/home/ludo/inbox/BAMs/ancient/NewBotai45_CGG_1_017037.Horse_nuc_wY.realigned.r.t.s.m.5p2.3p1.bam',
        'NewBotai46_CGG_1_017038': '/home/ludo/inbox/BAMs/ancient/NewBotai46_CGG_1_017038.Horse_nuc_wY.realigned.r.t.m.bam',
        'NUSTAR4_CGG_1_020438': '/home/ludo/inbox/BAMs/ancient/NUSTAR4_CGG_1_020438_i60_ATTAAA_U.r.t.5p2.bam',
        'NUSTAR5_': '/home/ludo/inbox/BAMs/ancient/NUSTAR5_.Horse_nuc_wY.realigned.r.t.m.bam',
        'OKG1_CGG_1_018397': '/home/ludo/inbox/BAMs/ancient/OKG1_CGG_1_018397.Horse_nuc_wY.realigned.r.t.bam',
        'OKG2_CGG_1_018398': '/home/ludo/inbox/BAMs/ancient/OKG2_CGG_1_018398.Horse_nuc_wY.realigned.r.t.m.3p2.bam',
        'Ote2_CGG_1_018473': '/home/ludo/inbox/BAMs/ancient/Ote2_CGG_1_018473.Horse_nuc_wY.realigned.r.t.m.1.3p.s.bam',
        'OTOK16_unregistered': '/home/ludo/inbox/BAMs/ancient/OTOK16_unregistered_i62_CAAAAT_U.r.t.bam',
        'PAVH11_CGG_1_018171': '/home/ludo/inbox/BAMs/ancient/PAVH11_CGG_1_018171.Horse_nuc_wY.realigned.r.t.m.bam',
        'PAVH2_CGG_1_018154': '/home/ludo/inbox/BAMs/ancient/PAVH2_CGG_1_018154.Horse_nuc_wY.realigned.r.t.m.bam',
        'PAVH4_CGG_1_018157': '/home/ludo/inbox/BAMs/ancient/PAVH4_CGG_1_018157.Horse_nuc_wY.realigned.r.t.m.bam',
        'PAVH6_CGG_1_018161': '/home/ludo/inbox/BAMs/ancient/PAVH6_CGG_1_018161.Horse_nuc_wY.realigned.r.t.m.bam',
        'PAVH8_CGG_1_018165': '/home/ludo/inbox/BAMs/ancient/PAVH8_CGG_1_018165.Horse_nuc_wY.realigned.r.t.m.bam',
        'PAVH9_CGG_1_018167': '/home/ludo/inbox/BAMs/ancient/PAVH9_CGG_1_018167.Horse_nuc_wY.realigned.r.t.m.bam',
        'Prze_0150A': '/home/ludo/inbox/BAMs/modern/Prze_0150A_KB3879.Horse_nuc_wY.realigned.bam',
        'Prze_0151A': '/home/ludo/inbox/BAMs/modern/Prze_0151A_KB7674.Horse_nuc_wY.realigned.bam',
        'Prze_0157A': '/home/ludo/inbox/BAMs/modern/Prze_0157A_SB293.Horse_nuc_wY.realigned.bam',
        'Prze_0158A': '/home/ludo/inbox/BAMs/modern/Prze_0158A_SB339.Horse_nuc_wY.realigned.bam',
        'Prze_0159A': '/home/ludo/inbox/BAMs/modern/Prze_0159A_SB4329.Horse_nuc_wY.realigned.bam',
        'Prze_0160A': '/home/ludo/inbox/BAMs/modern/Prze_0160A_SB533.Horse_nuc_wY.realigned.bam',
        'Prze_0213A_Paratype': '/home/ludo/inbox/BAMs/ancient/Prze_0213A_Paratype.Horse_nuc_wY.realigned.r.t.m.bam',
        'Quar_0073A': '/home/ludo/inbox/BAMs/modern/Quar_0073A_A2085.Horse_nuc_wY.realigned.bam',
        'Rid1_CGG_1_018468': '/home/ludo/inbox/BAMs/ancient/Rid1_CGG_1_018468.Horse_nuc_wY.realigned.r.t.s.m.3p1.bam',
        'Rid2_CGG_1_018469': '/home/ludo/inbox/BAMs/ancient/Rid2_CGG_1_018469.Horse_nuc_wY.realigned.r.t.m.5p2.3p2.bam',
        'Rus11_CGG_1_019162': '/home/ludo/inbox/BAMs/ancient/Rus11_CGG_1_019162_i08_U_ACTTGA.r.t.3p2.bam',
        'Rus14_CGG_1_019164': '/home/ludo/inbox/BAMs/ancient/Rus14_CGG_1_019164.Horse_nuc_wY.realigned.r.t.m.bam',
        'Rus16_CGG_1_019166': '/home/ludo/inbox/BAMs/ancient/Rus16_CGG_1_019166.Horse_nuc_wY.realigned.r.t.s.5p0.3p2.s.bam',
        'Rus19_CGG_1_019169': '/home/ludo/inbox/BAMs/ancient/Rus19_CGG_1_019169_i3_TTAGGC_U.r.t.5p1.3p2.bam',
        'Rus3_CGG_1_019154': '/home/ludo/inbox/BAMs/ancient/Rus3_CGG_1_019154_i32_ACGCAT_U.r.t.3p1.bam',
        'Rus37_CGG_1_019185': '/home/ludo/inbox/BAMs/ancient/Rus37_CGG_1_019185.Horse_nuc_wY.realigned.r.t.m.bam',
        'Rus38_CGG_1_019186': '/home/ludo/inbox/BAMs/ancient/Rus38_CGG_1_019186.Horse_nuc_wY.realigned.r.t.s.5p1.3p2.s.bam',
        'Rus41_CGG_1_019189': '/home/ludo/inbox/BAMs/ancient/Rus41_CGG_1_019189.Horse_nuc_wY.realigned.r.t.m.5p2.3p1.bam',
        'Rus45_CGG_1_019192': '/home/ludo/inbox/BAMs/ancient/Rus45_CGG_1_019192.Horse_nuc_wY.realigned.r.t.m.5p1.3p2.bam',
        'Rus48_CGG_1_019195': '/home/ludo/inbox/BAMs/ancient/Rus48_CGG_1_019195_i77_TCCCGG_U.r.t.5p1.3p2.bam',
        'Rus9_CGG_1_019160': '/home/ludo/inbox/BAMs/ancient/Rus9_CGG_1_019160.Horse_nuc_wY.realigned.r.t.m.3p2.s.bam',
        'Saa1_CGG_1_018474': '/home/ludo/inbox/BAMs/ancient/Saa1_CGG_1_018474.Horse_nuc_wY.realigned.r.t.m.bam',
        'SAG_S27_CGG_1_019559': '/home/ludo/inbox/BAMs/ancient/SAG_S27_CGG_1_019559.Horse_nuc_wY.realigned.r.t.m.5p2.3p2.bam',
        'Seb131_CGG_1_018487': '/home/ludo/inbox/BAMs/ancient/Seb131_CGG_1_018487_i1_ATCACG_U.r.t.5p2.3p1.bam',
        'Shet_0249A': '/home/ludo/inbox/BAMs/modern/Shet_0249A_SPH020.Horse_nuc_wY.realigned.bam',
        'Shet_0250A': '/home/ludo/inbox/BAMs/modern/Shet_0250A_SPH041.Horse_nuc_wY.realigned.bam',
        'Sorr_0236A': '/home/ludo/inbox/BAMs/modern/Sorr_0236A_SAMN02439778.Horse_nuc_wY.realigned.bam',
        'Spain38_CGG_1_020484': '/home/ludo/inbox/BAMs/ancient/Spain38_CGG_1_020484.Horse_nuc_wY.realigned.r.t.5p1.3p2.bam',
        'Spain39_CGG_1_020485': '/home/ludo/inbox/BAMs/ancient/Spain39_CGG_1_020485.Horse_nuc_wY.realigned.r.t.5p1.s.bam',
        'Stan_0081A': '/home/ludo/inbox/BAMs/modern/Stan_0081A_M5256.Horse_nuc_wY.realigned.bam',
        'Svi6_CGG_1_018375': '/home/ludo/inbox/BAMs/ancient/Svi6_CGG_1_018375.Horse_nuc_wY.realigned.r.t.m.5p2.3p2.bam',
        'Syr1_t1_c3_CGG_1_018919': '/home/ludo/inbox/BAMs/ancient/Syr1_t1_c3_CGG_1_018919.Horse_nuc_wY.realigned.r.t.m.bam',
        'Syr1_t1_c4_CGG_1_018920': '/home/ludo/inbox/BAMs/ancient/Syr1_t1_c4_CGG_1_018920_i52_AAACAA_U.r.t.5p4.3p4.bam',
        'Thor_0145A': '/home/ludo/inbox/BAMs/modern/Thor_0145A_Twilight.Horse_nuc_wY.realigned.bam',
        'Thor_0290A': '/home/ludo/inbox/BAMs/modern/Thor_0290A_SAMN01047706.Horse_nuc_wY.realigned.bam',
        'TP4_CGG_1_018394': '/home/ludo/inbox/BAMs/ancient/TP4_CGG_1_018394.Horse_nuc_wY.realigned.r.t.m.3p1.bam',
        'Tur140_CGG_1_018706': '/home/ludo/inbox/BAMs/ancient/Tur140_CGG_1_018706.Horse_nuc_wY.realigned.r.t.m.1.3p.bam',
        'Tur141_CGG_1_018707': '/home/ludo/inbox/BAMs/ancient/Tur141_CGG_1_018707.Horse_nuc_wY.realigned.r.t.m.5p2.3p1.bam',
        'Tur142_CGG_1_018708': '/home/ludo/inbox/BAMs/ancient/Tur142_CGG_1_018708.Horse_nuc_wY.realigned.r.t.m.5p2.3p1.s.bam',
        'Tur145_CGG_1_018711': '/home/ludo/inbox/BAMs/ancient/Tur145_CGG_1_018711.Horse_nuc_wY.realigned.r.t.m.1.3p.bam',
        'Tur146_CGG_1_018712': '/home/ludo/inbox/BAMs/ancient/Tur146_CGG_1_018712.Horse_nuc_wY.realigned.r.t.m.5p2.3p1.bam',
        'Tur150_CGG_1_018716': '/home/ludo/inbox/BAMs/ancient/Tur150_CGG_1_018716.Horse_nuc_wY.realigned.r.t.m.1.3p.bam',
        'Tur170_CGG_1_018736': '/home/ludo/inbox/BAMs/ancient/Tur170_CGG_1_018736.Horse_nuc_wY.realigned.r.t.m.3p1.bam',
        'Tur171_CGG_1_018737': '/home/ludo/inbox/BAMs/ancient/Tur171_CGG_1_018737.Horse_nuc_wY.realigned.r.t.m.1.3p.bam',
        'Tur172_CGG_1_018738': '/home/ludo/inbox/BAMs/ancient/Tur172_CGG_1_018738.Horse_nuc_wY.realigned.r.t.m.bam',
        'Tur173_CGG_1_018739': '/home/ludo/inbox/BAMs/ancient/Tur173_CGG_1_018739.Horse_nuc_wY.realigned.r.t.m.1.3p.bam',
        'Tur175_CGG_1_018741': '/home/ludo/inbox/BAMs/ancient/Tur175_CGG_1_018741.Horse_nuc_wY.realigned.r.t.m.bam',
        'Tur176_CGG_1_018742': '/home/ludo/inbox/BAMs/ancient/Tur176_CGG_1_018742.Horse_nuc_wY.realigned.r.t.m.bam',
        'Tur181_CGG_1_018747': '/home/ludo/inbox/BAMs/ancient/Tur181_CGG_1_018747.Horse_nuc_wY.realigned.r.t.m.bam',
        'Tur193_CGG_1_018759': '/home/ludo/inbox/BAMs/ancient/Tur193_CGG_1_018759.Horse_nuc_wY.realigned.r.t.m.1.3p.bam',
        'Tur194_CGG_1_018760': '/home/ludo/inbox/BAMs/ancient/Tur194_CGG_1_018760.Horse_nuc_wY.realigned.r.t.m.1.3p.bam',
        'Tur229_CGG_1_018795_Extraction2': '/home/ludo/inbox/BAMs/ancient/Tur229_CGG_1_018795.Horse_nuc_wY.realigned.r.t.m.1.3p.bam',
        'Tur243_CGG_1_018809': '/home/ludo/inbox/BAMs/ancient/Tur243_CGG_1_018809.Horse_nuc_wY.realigned.r.t.m.bam',
        'UCIE2012_85_': '/home/ludo/inbox/BAMs/ancient/UCIE2012_85_.Horse_nuc_wY.realigned.r.t.m.bam',
        'UE2275x2_CGG_1_020989': '/home/ludo/inbox/BAMs/ancient/UE2275x2_CGG_1_020989.Horse_nuc_wY.realigned.r.t.s.m.5p5.3p5.s.bam',
        'UE4618_CGG_1_020962': '/home/ludo/inbox/BAMs/ancient/UE4618_CGG_1_020962.Horse_nuc_wY.realigned.r.t.s.m.t.bam',
        'UK08_CGG_1_019435': '/home/ludo/inbox/BAMs/ancient/UK08_CGG_1_019435_i16_GTGTAT_U.r.t.5p2.3p1.bam',
        'UK15_CGG_1_019442': '/home/ludo/inbox/BAMs/ancient/UK15_CGG_1_019442.Horse_nuc_wY.realigned.r.t.m.3p1.bam',
        'UK16_CGG_1_019443': '/home/ludo/inbox/BAMs/ancient/UK16_CGG_1_019443.Horse_nuc_wY.realigned.r.t.m.bam',
        'UK17_CGG_1_019444': '/home/ludo/inbox/BAMs/ancient/UK17_CGG_1_019444.Horse_nuc_wY.realigned.r.t.m.bam',
        'UK18_CGG_1_019445': '/home/ludo/inbox/BAMs/ancient/UK18_CGG_1_019445.Horse_nuc_wY.realigned.r.t.m.bam',
        'UK19_CGG_1_019446': '/home/ludo/inbox/BAMs/ancient/UK19_CGG_1_019446.Horse_nuc_wY.realigned.r.t.m.5p2.3p1.bam',
        'UK20_CGG_1_019447': '/home/ludo/inbox/BAMs/ancient/UK20_CGG_1_019447_i11_GGCTAC_U.r.t.5p2.3p1.bam',
        'Upps02_CGG_1_018490': '/home/ludo/inbox/BAMs/ancient/Upps02_CGG_1_018490.Horse_nuc_wY.realigned.r.t.s.m.5p2.3p2.bam',
        'Vert293_CGG_1_018522': '/home/ludo/inbox/BAMs/ancient/Vert293_CGG_1_018522.Horse_nuc_wY.realigned.r.t.m.bam',
        'Vert300_CGG_1_018529': '/home/ludo/inbox/BAMs/ancient/Vert300_CGG_1_018529_i6_GCCAAT_U.r.t.3p1.bam',
        'Vert304_CGG_1_018533': '/home/ludo/inbox/BAMs/ancient/Vert304_CGG_1_018533_i10_TAGCTT_U.r.t.5p1.3p2.bam',
        'Vert311_CGG_1_018540': '/home/ludo/inbox/BAMs/ancient/Vert311_CGG_1_018540.Horse_nuc_wY.realigned.r.t.m.3p1.bam',
        'VHR010_CGG_1_020949': '/home/ludo/inbox/BAMs/ancient/VHR010_CGG_1_020949_i74_GTGGGG_U.r.t.s.5p1.3p2.bam',
        'VHR011_CGG_1_020950': '/home/ludo/inbox/BAMs/ancient/VHR011_CGG_1_020950_i76_TCAGCT_U.r.t.bam',
        'VHR017_CGG_1_020952': '/home/ludo/inbox/BAMs/ancient/VHR017_CGG_1_020952_i62_CAAAAT_U.r.t.s.3p1.bam',
        'VHR031_CGG_1_020955': '/home/ludo/inbox/BAMs/ancient/VHR031_CGG_1_020955_i61_ATTCTC_U.r.t.3p1.bam',
        'VHR037_CGG_1_020957': '/home/ludo/inbox/BAMs/ancient/VHR037_CGG_1_020957_i60_ATTAAA_U.r.t.5p2.3p1.bam',
        'VHR062_CGG_1_020959': '/home/ludo/inbox/BAMs/ancient/VHR062_CGG_1_020959_i75_TAGTAA_U.r.t.3p1.bam',
        'VHR102_CGG_1_020961': '/home/ludo/inbox/BAMs/ancient/VHR102_CGG_1_020961_i54_ACCATC_U.r.t.3p1.bam',
        'VIR175_CGG_1_016987': '/home/ludo/inbox/BAMs/ancient/VIR175_CGG_1_016987_i41_TGTCTG.r.t.3p2.bam',
        'Yaku_0163A': '/home/ludo/inbox/BAMs/modern/Yaku_0163A_Yak1.Horse_nuc_wY.realigned.bam',
        'Yaku_0170A': '/home/ludo/inbox/BAMs/modern/Yaku_0170A_Yak8.Horse_nuc_wY.realigned.bam',
        'Yaku_0171A': '/home/ludo/inbox/BAMs/modern/Yaku_0171A_Yak9.Horse_nuc_wY.realigned.bam',
        'YER28_CGG_1_020254': '/home/ludo/inbox/BAMs/ancient/YER28_CGG_1_020254.Horse_nuc_wY.realigned.r.t.m.1.3p.s.bam',
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


# how many CPU cores does this machine have
TOTAL_CORES = cpu_count()

# set how many cores a single working can use
CPU_CORES_ONE = 1
CPU_CORES_LOW = int(TOTAL_CORES * 0.1)   # 10%
CPU_CORES_MED = int(TOTAL_CORES * 0.25)  # 25%
CPU_CORES_HIGH = int(TOTAL_CORES * 0.5)  # 50%
CPU_CORES_MAX = int(TOTAL_CORES * 0.9)   # 90%

# the minimum derived allele frequency of modern SNPs to include
MIN_DAF = 0.05

# the minimum phred scaled genotype quality (30 = 99.9%)
MIN_GENO_QUAL = 30
