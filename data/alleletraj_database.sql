# ************************************************************
# Sequel Pro SQL dump
# Version 4541
#
# http://www.sequelpro.com/
# https://github.com/sequelpro/sequelpro
#
# Host: 127.0.0.1 (MySQL 5.6.19)
# Database: alleletraj_horse_equcab2_rel38
# Generation Time: 2019-06-18 16:48:50 +0000
# ************************************************************


/*!40101 SET @OLD_CHARACTER_SET_CLIENT=@@CHARACTER_SET_CLIENT */;
/*!40101 SET @OLD_CHARACTER_SET_RESULTS=@@CHARACTER_SET_RESULTS */;
/*!40101 SET @OLD_COLLATION_CONNECTION=@@COLLATION_CONNECTION */;
/*!40101 SET NAMES utf8 */;
/*!40014 SET @OLD_FOREIGN_KEY_CHECKS=@@FOREIGN_KEY_CHECKS, FOREIGN_KEY_CHECKS=0 */;
/*!40101 SET @OLD_SQL_MODE=@@SQL_MODE, SQL_MODE='NO_AUTO_VALUE_ON_ZERO' */;
/*!40111 SET @OLD_SQL_NOTES=@@SQL_NOTES, SQL_NOTES=0 */;


# Dump of table ascertainment
# ------------------------------------------------------------

CREATE TABLE `ascertainment` (
  `id` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `qtl_id` int(11) unsigned DEFAULT NULL,
  `type` varchar(255) NOT NULL,
  `rsnumber` varchar(255) NOT NULL,
  `chrom` char(2) NOT NULL,
  `site` int(11) NOT NULL,
  `ref` char(1) NOT NULL,
  `alt` char(1) NOT NULL,
  `chip_name` varchar(255) DEFAULT '',
  `snp_name` varchar(510) DEFAULT '',
  PRIMARY KEY (`id`),
  UNIQUE KEY `qtl_chrom_site` (`qtl_id`,`chrom`,`site`),
  KEY `chrom_site` (`chrom`,`site`),
  KEY `rsnumber` (`rsnumber`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8;



# Dump of table ensembl_genes
# ------------------------------------------------------------

CREATE TABLE `ensembl_genes` (
  `id` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `source` varchar(255) NOT NULL,
  `gene_id` varchar(255) NOT NULL,
  `gene_name` varchar(255) DEFAULT '',
  `version` int(11) NOT NULL,
  `biotype` varchar(255) NOT NULL,
  `chrom` varchar(15) NOT NULL,
  `start` int(11) NOT NULL,
  `end` int(11) NOT NULL,
  PRIMARY KEY (`id`),
  UNIQUE KEY `gene_id` (`gene_id`),
  KEY `chrom_start_end` (`chrom`,`start`,`end`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8;



# Dump of table ensembl_variants
# ------------------------------------------------------------

CREATE TABLE `ensembl_variants` (
  `id` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `dbxref` varchar(255) NOT NULL,
  `rsnumber` varchar(255) NOT NULL,
  `type` char(20) NOT NULL,
  `chrom` varchar(10) NOT NULL,
  `start` int(11) NOT NULL,
  `end` int(11) NOT NULL,
  `ref` varchar(100) NOT NULL,
  `alt` varchar(250) NOT NULL,
  `indel` tinyint(1) DEFAULT NULL,
  KEY `id` (`id`),
  KEY `rsnumber` (`rsnumber`),
  KEY `chrom_start_end` (`chrom`,`start`,`end`),
  KEY `type` (`type`),
  KEY `indel` (`indel`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8
/*!50500 PARTITION BY LIST  COLUMNS(chrom)
(PARTITION p1 VALUES IN ('1') ENGINE = InnoDB,
 PARTITION p2 VALUES IN ('2') ENGINE = InnoDB,
 PARTITION p3 VALUES IN ('3') ENGINE = InnoDB,
 PARTITION p4 VALUES IN ('4') ENGINE = InnoDB,
 PARTITION p5 VALUES IN ('5') ENGINE = InnoDB,
 PARTITION p6 VALUES IN ('6') ENGINE = InnoDB,
 PARTITION p7 VALUES IN ('7') ENGINE = InnoDB,
 PARTITION p8 VALUES IN ('8') ENGINE = InnoDB,
 PARTITION p9 VALUES IN ('9') ENGINE = InnoDB,
 PARTITION p10 VALUES IN ('10') ENGINE = InnoDB,
 PARTITION p11 VALUES IN ('11') ENGINE = InnoDB,
 PARTITION p12 VALUES IN ('12') ENGINE = InnoDB,
 PARTITION p13 VALUES IN ('13') ENGINE = InnoDB,
 PARTITION p14 VALUES IN ('14') ENGINE = InnoDB,
 PARTITION p15 VALUES IN ('15') ENGINE = InnoDB,
 PARTITION p16 VALUES IN ('16') ENGINE = InnoDB,
 PARTITION p17 VALUES IN ('17') ENGINE = InnoDB,
 PARTITION p18 VALUES IN ('18') ENGINE = InnoDB,
 PARTITION p19 VALUES IN ('19') ENGINE = InnoDB,
 PARTITION p20 VALUES IN ('20') ENGINE = InnoDB,
 PARTITION p21 VALUES IN ('21') ENGINE = InnoDB,
 PARTITION p22 VALUES IN ('22') ENGINE = InnoDB,
 PARTITION p23 VALUES IN ('23') ENGINE = InnoDB,
 PARTITION p24 VALUES IN ('24') ENGINE = InnoDB,
 PARTITION p25 VALUES IN ('25') ENGINE = InnoDB,
 PARTITION p26 VALUES IN ('26') ENGINE = InnoDB,
 PARTITION p27 VALUES IN ('27') ENGINE = InnoDB,
 PARTITION p28 VALUES IN ('28') ENGINE = InnoDB,
 PARTITION p29 VALUES IN ('29') ENGINE = InnoDB,
 PARTITION p30 VALUES IN ('30') ENGINE = InnoDB,
 PARTITION p31 VALUES IN ('31') ENGINE = InnoDB,
 PARTITION pX VALUES IN ('X') ENGINE = InnoDB,
 PARTITION pY VALUES IN ('Y') ENGINE = InnoDB,
 PARTITION pMT VALUES IN ('MT') ENGINE = InnoDB) */;



# Dump of table modern_snp_daf
# ------------------------------------------------------------

CREATE TABLE `modern_snp_daf` (
  `id` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `modsnp_id` int(11) unsigned NOT NULL,
  `population` varchar(255) NOT NULL,
  `ancestral_count` int(11) unsigned NOT NULL,
  `derived_count` int(11) unsigned NOT NULL,
  `daf` float NOT NULL,
  PRIMARY KEY (`id`),
  UNIQUE KEY `modsnp_population` (`modsnp_id`,`population`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8;



# Dump of table modern_snps
# ------------------------------------------------------------

CREATE TABLE `modern_snps` (
  `id` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `chrom` char(2) NOT NULL,
  `site` int(11) NOT NULL,
  `ancestral` char(1) NOT NULL,
  `derived` char(1) NOT NULL,
  `type` char(2) NOT NULL,
  `variant_id` int(11) DEFAULT NULL,
  `snpchip_id` int(11) DEFAULT NULL,
  `gene_id` int(11) unsigned DEFAULT NULL,
  `neutral` tinyint(1) DEFAULT NULL,
  `mispolar` tinyint(1) DEFAULT NULL,
  KEY `id` (`id`),
  KEY `chrom_site` (`chrom`,`site`),
  KEY `variant_id` (`variant_id`),
  KEY `snpchip_id` (`snpchip_id`),
  KEY `gene_id` (`gene_id`),
  KEY `neutral` (`neutral`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8
/*!50500 PARTITION BY LIST  COLUMNS(chrom)
(PARTITION p1 VALUES IN ('1') ENGINE = InnoDB,
 PARTITION p2 VALUES IN ('2') ENGINE = InnoDB,
 PARTITION p3 VALUES IN ('3') ENGINE = InnoDB,
 PARTITION p4 VALUES IN ('4') ENGINE = InnoDB,
 PARTITION p5 VALUES IN ('5') ENGINE = InnoDB,
 PARTITION p6 VALUES IN ('6') ENGINE = InnoDB,
 PARTITION p7 VALUES IN ('7') ENGINE = InnoDB,
 PARTITION p8 VALUES IN ('8') ENGINE = InnoDB,
 PARTITION p9 VALUES IN ('9') ENGINE = InnoDB,
 PARTITION p10 VALUES IN ('10') ENGINE = InnoDB,
 PARTITION p11 VALUES IN ('11') ENGINE = InnoDB,
 PARTITION p12 VALUES IN ('12') ENGINE = InnoDB,
 PARTITION p13 VALUES IN ('13') ENGINE = InnoDB,
 PARTITION p14 VALUES IN ('14') ENGINE = InnoDB,
 PARTITION p15 VALUES IN ('15') ENGINE = InnoDB,
 PARTITION p16 VALUES IN ('16') ENGINE = InnoDB,
 PARTITION p17 VALUES IN ('17') ENGINE = InnoDB,
 PARTITION p18 VALUES IN ('18') ENGINE = InnoDB,
 PARTITION p19 VALUES IN ('19') ENGINE = InnoDB,
 PARTITION p20 VALUES IN ('20') ENGINE = InnoDB,
 PARTITION p21 VALUES IN ('21') ENGINE = InnoDB,
 PARTITION p22 VALUES IN ('22') ENGINE = InnoDB,
 PARTITION p23 VALUES IN ('23') ENGINE = InnoDB,
 PARTITION p24 VALUES IN ('24') ENGINE = InnoDB,
 PARTITION p25 VALUES IN ('25') ENGINE = InnoDB,
 PARTITION p26 VALUES IN ('26') ENGINE = InnoDB,
 PARTITION p27 VALUES IN ('27') ENGINE = InnoDB,
 PARTITION p28 VALUES IN ('28') ENGINE = InnoDB,
 PARTITION p29 VALUES IN ('29') ENGINE = InnoDB,
 PARTITION p30 VALUES IN ('30') ENGINE = InnoDB,
 PARTITION p31 VALUES IN ('31') ENGINE = InnoDB,
 PARTITION pX VALUES IN ('X') ENGINE = InnoDB,
 PARTITION pY VALUES IN ('Y') ENGINE = InnoDB,
 PARTITION pMT VALUES IN ('MT') ENGINE = InnoDB) */;



# Dump of table pubmeds
# ------------------------------------------------------------

CREATE TABLE `pubmeds` (
  `id` int(11) unsigned NOT NULL,
  `authors` text,
  `year` char(255) DEFAULT NULL,
  `title` text,
  `journal` char(255) DEFAULT NULL,
  PRIMARY KEY (`id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8;



# Dump of table qtl_snps
# ------------------------------------------------------------

CREATE TABLE `qtl_snps` (
  `id` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `qtl_id` int(11) unsigned NOT NULL,
  `modsnp_id` int(11) unsigned NOT NULL,
  `num_reads` int(11) DEFAULT NULL,
  `best` tinyint(1) DEFAULT NULL,
  PRIMARY KEY (`id`),
  UNIQUE KEY `qtl_modsnp` (`qtl_id`,`modsnp_id`),
  KEY `best` (`best`),
  KEY `num_reads` (`num_reads`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8;



# Dump of table qtl_stats
# ------------------------------------------------------------

CREATE TABLE `qtl_stats` (
  `id` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `qtl_id` int(11) unsigned DEFAULT NULL,
  `chrom` char(2) DEFAULT NULL,
  `class` varchar(255) DEFAULT NULL,
  `type` varchar(255) DEFAULT NULL,
  `name` varchar(255) DEFAULT NULL,
  `pvalue` varchar(255) DEFAULT NULL,
  `significance` varchar(255) DEFAULT NULL,
  `snps` int(11) NOT NULL DEFAULT '0',
  `max_samples` int(11) DEFAULT NULL,
  `avg_samples` decimal(6,2) DEFAULT NULL,
  `max_reads` int(11) DEFAULT NULL,
  `avg_reads` decimal(6,2) DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `qtl_id` (`qtl_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8;



# Dump of table qtls
# ------------------------------------------------------------

CREATE TABLE `qtls` (
  `id` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `qtldb_id` int(11) unsigned DEFAULT NULL,
  `associationType` char(255) DEFAULT NULL,
  `symbol` char(255) DEFAULT NULL,
  `pubmed_id` int(11) unsigned DEFAULT NULL,
  `trait_id` int(11) unsigned DEFAULT NULL,
  `chrom` char(2) DEFAULT NULL,
  `A1` char(255) DEFAULT NULL,
  `A2` char(255) DEFAULT NULL,
  `peak` char(255) DEFAULT NULL,
  `B1` char(255) DEFAULT NULL,
  `B2` char(255) DEFAULT NULL,
  `linkageLoc_peak` float DEFAULT NULL,
  `linkageLoc_start` float DEFAULT NULL,
  `linkageLoc_end` float DEFAULT NULL,
  `genomeLoc_start` int(11) DEFAULT NULL,
  `genomeLoc_end` int(11) DEFAULT NULL,
  `testBase` char(255) DEFAULT NULL,
  `mapType` char(255) DEFAULT NULL,
  `testModel` char(255) DEFAULT NULL,
  `LSmean` char(255) DEFAULT NULL,
  `Pvalue` char(255) DEFAULT NULL,
  `LODscore` char(255) DEFAULT NULL,
  `Fstat` char(255) DEFAULT NULL,
  `varance` char(255) DEFAULT NULL,
  `LR` char(255) DEFAULT NULL,
  `bayes` char(255) DEFAULT NULL,
  `significance` char(255) DEFAULT NULL,
  `domEffect` char(255) DEFAULT NULL,
  `addEffect` char(255) DEFAULT NULL,
  `gene_id` char(255) DEFAULT NULL,
  `geneSource` char(255) DEFAULT NULL,
  `valid` tinyint(1) DEFAULT NULL,
  `site` int(11) DEFAULT NULL,
  `start` int(11) DEFAULT NULL,
  `end` int(11) DEFAULT NULL,
  PRIMARY KEY (`id`),
  UNIQUE KEY `qtldb_id` (`qtldb_id`),
  KEY `trait_fk` (`trait_id`),
  KEY `pubmed_fk` (`pubmed_id`),
  KEY `peak` (`peak`),
  KEY `chrom_start_end` (`chrom`,`start`,`end`),
  KEY `valid` (`valid`),
  KEY `associationtype_valid` (`associationType`,`valid`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8;



# Dump of table sample_bins
# ------------------------------------------------------------

CREATE TABLE `sample_bins` (
  `id` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `name` varchar(255) NOT NULL,
  `max` int(11) NOT NULL,
  `min` int(11) NOT NULL,
  `num_samples` int(11) DEFAULT NULL,
  PRIMARY KEY (`id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8;



# Dump of table sample_reads
# ------------------------------------------------------------

CREATE TABLE `sample_reads` (
  `id` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `sample_id` int(11) unsigned NOT NULL,
  `chrom` varchar(2) NOT NULL,
  `site` int(11) NOT NULL,
  `base` varchar(1) NOT NULL,
  `dist` int(11) DEFAULT NULL,
  `mapq` int(11) DEFAULT NULL,
  `baseq` int(11) DEFAULT NULL,
  `genoq` int(11) DEFAULT NULL,
  KEY `id` (`id`),
  KEY `sampleID` (`sample_id`),
  KEY `chrom_site` (`chrom`,`site`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8
/*!50500 PARTITION BY LIST  COLUMNS(chrom)
(PARTITION p1 VALUES IN ('1') ENGINE = InnoDB,
 PARTITION p2 VALUES IN ('2') ENGINE = InnoDB,
 PARTITION p3 VALUES IN ('3') ENGINE = InnoDB,
 PARTITION p4 VALUES IN ('4') ENGINE = InnoDB,
 PARTITION p5 VALUES IN ('5') ENGINE = InnoDB,
 PARTITION p6 VALUES IN ('6') ENGINE = InnoDB,
 PARTITION p7 VALUES IN ('7') ENGINE = InnoDB,
 PARTITION p8 VALUES IN ('8') ENGINE = InnoDB,
 PARTITION p9 VALUES IN ('9') ENGINE = InnoDB,
 PARTITION p10 VALUES IN ('10') ENGINE = InnoDB,
 PARTITION p11 VALUES IN ('11') ENGINE = InnoDB,
 PARTITION p12 VALUES IN ('12') ENGINE = InnoDB,
 PARTITION p13 VALUES IN ('13') ENGINE = InnoDB,
 PARTITION p14 VALUES IN ('14') ENGINE = InnoDB,
 PARTITION p15 VALUES IN ('15') ENGINE = InnoDB,
 PARTITION p16 VALUES IN ('16') ENGINE = InnoDB,
 PARTITION p17 VALUES IN ('17') ENGINE = InnoDB,
 PARTITION p18 VALUES IN ('18') ENGINE = InnoDB,
 PARTITION p19 VALUES IN ('19') ENGINE = InnoDB,
 PARTITION p20 VALUES IN ('20') ENGINE = InnoDB,
 PARTITION p21 VALUES IN ('21') ENGINE = InnoDB,
 PARTITION p22 VALUES IN ('22') ENGINE = InnoDB,
 PARTITION p23 VALUES IN ('23') ENGINE = InnoDB,
 PARTITION p24 VALUES IN ('24') ENGINE = InnoDB,
 PARTITION p25 VALUES IN ('25') ENGINE = InnoDB,
 PARTITION p26 VALUES IN ('26') ENGINE = InnoDB,
 PARTITION p27 VALUES IN ('27') ENGINE = InnoDB,
 PARTITION p28 VALUES IN ('28') ENGINE = InnoDB,
 PARTITION p29 VALUES IN ('29') ENGINE = InnoDB,
 PARTITION p30 VALUES IN ('30') ENGINE = InnoDB,
 PARTITION p31 VALUES IN ('31') ENGINE = InnoDB,
 PARTITION pX VALUES IN ('X') ENGINE = InnoDB,
 PARTITION pY VALUES IN ('Y') ENGINE = InnoDB,
 PARTITION pMT VALUES IN ('MT') ENGINE = InnoDB) */;



# Dump of table sample_runs
# ------------------------------------------------------------

CREATE TABLE `sample_runs` (
  `id` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `sample_id` int(11) NOT NULL,
  `bioproject` varchar(255) DEFAULT '',
  `biosample` varchar(255) NOT NULL DEFAULT '',
  `accession` varchar(255) NOT NULL,
  `libname` varchar(255) DEFAULT '',
  `paired` tinyint(1) NOT NULL,
  PRIMARY KEY (`id`),
  UNIQUE KEY `accession` (`accession`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8;



# Dump of table samples
# ------------------------------------------------------------

CREATE TABLE `samples` (
  `id` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `name` varchar(255) NOT NULL,
  `population` varchar(255) NOT NULL,
  `ancient` tinyint(1) NOT NULL,
  `alias` varchar(255) DEFAULT '',
  `group_a` varchar(255) DEFAULT NULL,
  `group_b` varchar(255) DEFAULT NULL,
  `site` varchar(255) DEFAULT NULL,
  `location` varchar(255) DEFAULT NULL,
  `bp_max` int(11) DEFAULT NULL,
  `bp_min` int(11) DEFAULT NULL,
  `bp_median` int(11) DEFAULT NULL,
  `bin_id` int(11) DEFAULT NULL,
  `age` varchar(255) DEFAULT NULL,
  `period` varchar(255) DEFAULT NULL,
  `lat` varchar(255) DEFAULT NULL,
  `long` varchar(255) DEFAULT NULL,
  `sex` char(1) DEFAULT '',
  `mtdna` varchar(2) DEFAULT NULL,
  `sfs` tinyint(1) DEFAULT NULL,
  `coverage` float DEFAULT NULL,
  `path` varchar(255) DEFAULT NULL,
  PRIMARY KEY (`id`),
  UNIQUE KEY `name` (`name`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8;



# Dump of table selection
# ------------------------------------------------------------

CREATE TABLE `selection` (
  `id` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `population` varchar(255) NOT NULL DEFAULT '',
  `modsnp_id` int(11) unsigned NOT NULL,
  `length` int(11) unsigned NOT NULL,
  `thin` int(11) unsigned NOT NULL,
  `model` decimal(1,1) NOT NULL,
  `no_modern` tinyint(1) NOT NULL,
  `mispolar` tinyint(1) NOT NULL,
  PRIMARY KEY (`id`),
  KEY `modsnp_id` (`modsnp_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8;



# Dump of table selection_ess
# ------------------------------------------------------------

CREATE TABLE `selection_ess` (
  `id` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `selection_id` int(11) NOT NULL,
  `chain` int(11) DEFAULT NULL,
  `lnL` float NOT NULL,
  `pathlnL` float NOT NULL,
  `alpha1` float NOT NULL,
  `alpha2` float NOT NULL,
  `F` float NOT NULL,
  `age` float NOT NULL,
  `end_freq` float NOT NULL,
  `sample_time_0` float DEFAULT NULL,
  `sample_time_1` float DEFAULT NULL,
  `sample_time_2` float DEFAULT NULL,
  `sample_time_3` float DEFAULT NULL,
  `sample_time_4` float DEFAULT NULL,
  `sample_time_5` float DEFAULT NULL,
  `sample_time_6` float DEFAULT NULL,
  `sample_time_7` float DEFAULT NULL,
  `sample_time_8` float DEFAULT NULL,
  `sample_time_9` float DEFAULT NULL,
  `sample_time_10` float DEFAULT NULL,
  `sample_time_11` float DEFAULT NULL,
  `sample_time_12` float DEFAULT NULL,
  `sample_time_13` float DEFAULT NULL,
  `sample_time_14` float DEFAULT NULL,
  `sample_time_15` float DEFAULT NULL,
  `sample_time_16` float DEFAULT NULL,
  `sample_time_17` float DEFAULT NULL,
  `first_nonzero` float NOT NULL,
  PRIMARY KEY (`id`),
  KEY `selection_id` (`selection_id`,`chain`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8;



# Dump of table selection_map
# ------------------------------------------------------------

CREATE TABLE `selection_map` (
  `id` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `selection_id` int(11) NOT NULL,
  `chain` int(11) NOT NULL,
  `lnL` float NOT NULL,
  `pathlnL` float NOT NULL,
  `alpha1` float NOT NULL,
  `alpha2` float NOT NULL,
  `F` float NOT NULL,
  `age` int(11) unsigned NOT NULL,
  `end_freq` float NOT NULL,
  `sample_time_0` int(11) unsigned DEFAULT NULL,
  `sample_time_1` int(11) unsigned DEFAULT NULL,
  `sample_time_2` int(11) unsigned DEFAULT NULL,
  `sample_time_3` int(11) unsigned DEFAULT NULL,
  `sample_time_4` int(11) unsigned DEFAULT NULL,
  `sample_time_5` int(11) unsigned DEFAULT NULL,
  `sample_time_6` int(11) unsigned DEFAULT NULL,
  `sample_time_7` int(11) unsigned DEFAULT NULL,
  `sample_time_8` int(11) unsigned DEFAULT NULL,
  `sample_time_9` int(11) unsigned DEFAULT NULL,
  `sample_time_10` int(11) unsigned DEFAULT NULL,
  `sample_time_11` int(11) unsigned DEFAULT NULL,
  `sample_time_12` int(11) unsigned DEFAULT NULL,
  `sample_time_13` int(11) unsigned DEFAULT NULL,
  `sample_time_14` int(11) unsigned DEFAULT NULL,
  `sample_time_15` int(11) unsigned DEFAULT NULL,
  `sample_time_16` int(11) unsigned DEFAULT NULL,
  `sample_time_17` int(11) unsigned DEFAULT NULL,
  `first_nonzero` int(11) unsigned NOT NULL,
  PRIMARY KEY (`id`),
  UNIQUE KEY `selection_id` (`selection_id`,`chain`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8;



# Dump of table selection_psrf
# ------------------------------------------------------------

CREATE TABLE `selection_psrf` (
  `id` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `selection_id` int(11) NOT NULL,
  `mpsrf` float NOT NULL,
  `lnL` float NOT NULL,
  `pathlnL` float NOT NULL,
  `alpha1` float NOT NULL,
  `alpha2` float NOT NULL,
  `F` float NOT NULL,
  `age` float NOT NULL,
  `end_freq` float NOT NULL,
  `sample_time_0` float DEFAULT NULL,
  `sample_time_1` float DEFAULT NULL,
  `sample_time_2` float DEFAULT NULL,
  `sample_time_3` float DEFAULT NULL,
  `sample_time_4` float DEFAULT NULL,
  `sample_time_5` float DEFAULT NULL,
  `sample_time_6` float DEFAULT NULL,
  `sample_time_7` float DEFAULT NULL,
  `sample_time_8` float DEFAULT NULL,
  `sample_time_9` float DEFAULT NULL,
  `sample_time_10` float DEFAULT NULL,
  `sample_time_11` float DEFAULT NULL,
  `sample_time_12` float DEFAULT NULL,
  `sample_time_13` float DEFAULT NULL,
  `sample_time_14` float DEFAULT NULL,
  `sample_time_15` float DEFAULT NULL,
  `sample_time_16` float DEFAULT NULL,
  `sample_time_17` float DEFAULT NULL,
  `first_nonzero` float NOT NULL,
  PRIMARY KEY (`id`),
  UNIQUE KEY `selection_id` (`selection_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8;



# Dump of table snpchip
# ------------------------------------------------------------

CREATE TABLE `snpchip` (
  `id` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `chip_name` varchar(255) NOT NULL,
  `rsnumber` varchar(255) DEFAULT '',
  `snp_name` varchar(255) NOT NULL,
  PRIMARY KEY (`id`),
  UNIQUE KEY `chip_rsnumber_snp` (`chip_name`,`rsnumber`,`snp_name`),
  KEY `snp_name` (`snp_name`),
  KEY `rsnumber` (`rsnumber`),
  KEY `chip_chrom_site` (`chip_name`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8;



# Dump of table snpchip_axiom
# ------------------------------------------------------------

CREATE TABLE `snpchip_axiom` (
  `id` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `snp_name` varchar(255) DEFAULT NULL,
  `chrom` char(2) DEFAULT NULL,
  `site` int(11) DEFAULT NULL,
  PRIMARY KEY (`id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8;



# Dump of table sweep_snps
# ------------------------------------------------------------

CREATE TABLE `sweep_snps` (
  `id` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `qtl_id` int(11) NOT NULL,
  `chrom` char(2) NOT NULL,
  `site` int(11) NOT NULL,
  `cdf` float NOT NULL,
  `p` float NOT NULL,
  PRIMARY KEY (`id`),
  KEY `qtl_id` (`qtl_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8;



# Dump of table traits
# ------------------------------------------------------------

CREATE TABLE `traits` (
  `id` int(11) unsigned NOT NULL,
  `class` varchar(255) DEFAULT NULL,
  `type` varchar(255) DEFAULT NULL,
  `name` varchar(255) DEFAULT NULL,
  `VT` varchar(255) DEFAULT NULL,
  `CMO` varchar(255) DEFAULT NULL,
  `LPT` varchar(255) DEFAULT NULL,
  PRIMARY KEY (`id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8;




/*!40111 SET SQL_NOTES=@OLD_SQL_NOTES */;
/*!40101 SET SQL_MODE=@OLD_SQL_MODE */;
/*!40014 SET FOREIGN_KEY_CHECKS=@OLD_FOREIGN_KEY_CHECKS */;
/*!40101 SET CHARACTER_SET_CLIENT=@OLD_CHARACTER_SET_CLIENT */;
/*!40101 SET CHARACTER_SET_RESULTS=@OLD_CHARACTER_SET_RESULTS */;
/*!40101 SET COLLATION_CONNECTION=@OLD_COLLATION_CONNECTION */;
