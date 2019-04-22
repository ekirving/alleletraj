# ************************************************************
# Sequel Pro SQL dump
# Version 4541
#
# http://www.sequelpro.com/
# https://github.com/sequelpro/sequelpro
#
# Host: 127.0.0.1 (MySQL 5.6.19)
# Database: alleletraj_horse_equcab2_rel37
# Generation Time: 2019-04-22 08:26:39 +0000
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
  `type` varchar(255) NOT NULL DEFAULT '',
  `rsnumber` varchar(255) NOT NULL DEFAULT '',
  `chrom` char(2) NOT NULL DEFAULT '',
  `site` int(11) NOT NULL,
  `ref` char(1) NOT NULL DEFAULT '',
  `alt` char(1) NOT NULL DEFAULT '',
  `chip_name` varchar(255) DEFAULT '',
  `snp_name` varchar(510) DEFAULT '',
  PRIMARY KEY (`id`),
  UNIQUE KEY `qtl_id` (`qtl_id`,`chrom`,`site`),
  KEY `chrom` (`chrom`,`site`),
  KEY `rsnumber` (`rsnumber`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8;



# Dump of table ascertainment_coverage
# ------------------------------------------------------------

CREATE TABLE `ascertainment_coverage` (
  `sample_id` int(11) unsigned NOT NULL,
  `perc` decimal(24,4) DEFAULT NULL
) ENGINE=InnoDB DEFAULT CHARSET=utf8;



# Dump of table ascertainment_jake
# ------------------------------------------------------------

CREATE TABLE `ascertainment_jake` (
  `id` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `rsnumber` varchar(255) NOT NULL DEFAULT '',
  `cat` varchar(10) NOT NULL DEFAULT '',
  PRIMARY KEY (`id`),
  KEY `rsnumber` (`rsnumber`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8;



# Dump of table ascertainment_old
# ------------------------------------------------------------

CREATE TABLE `ascertainment_old` (
  `id` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `qtl_id` int(11) unsigned DEFAULT NULL,
  `type` varchar(255) NOT NULL DEFAULT '',
  `rsnumber` varchar(255) NOT NULL DEFAULT '',
  `chrom` char(2) NOT NULL DEFAULT '',
  `site` int(11) NOT NULL,
  `ref` char(1) NOT NULL DEFAULT '',
  `alt` char(1) NOT NULL DEFAULT '',
  `chip_name` varchar(255) DEFAULT '',
  `snp_name` varchar(510) DEFAULT '',
  PRIMARY KEY (`id`),
  UNIQUE KEY `qtl_id` (`qtl_id`,`chrom`,`site`),
  KEY `chrom` (`chrom`,`site`),
  KEY `rsnumber` (`rsnumber`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8;



# Dump of table ensembl_genes
# ------------------------------------------------------------

CREATE TABLE `ensembl_genes` (
  `id` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `source` varchar(255) NOT NULL DEFAULT '',
  `gene_id` varchar(255) NOT NULL DEFAULT '',
  `gene_name` varchar(255) DEFAULT '',
  `version` int(11) NOT NULL,
  `biotype` varchar(255) NOT NULL DEFAULT '',
  `chrom` varchar(10) NOT NULL DEFAULT '',
  `start` int(11) NOT NULL,
  `end` int(11) NOT NULL,
  PRIMARY KEY (`id`),
  UNIQUE KEY `gene_id` (`gene_id`),
  KEY `chrom` (`chrom`,`start`,`end`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8;



# Dump of table ensembl_variants
# ------------------------------------------------------------

CREATE TABLE `ensembl_variants` (
  `id` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `dbxref` varchar(255) NOT NULL DEFAULT '',
  `rsnumber` varchar(255) NOT NULL DEFAULT '',
  `type` char(15) NOT NULL DEFAULT '',
  `chrom` varchar(10) NOT NULL DEFAULT '',
  `start` int(11) NOT NULL,
  `end` int(11) NOT NULL,
  `ref` varchar(100) NOT NULL DEFAULT '',
  `alt` varchar(250) NOT NULL DEFAULT '',
  `indel` tinyint(1) DEFAULT NULL,
  KEY `id` (`id`),
  KEY `rsnumber` (`rsnumber`),
  KEY `chrom` (`chrom`,`start`,`end`),
  KEY `type` (`type`),
  KEY `indel` (`indel`),
  KEY `id_2` (`id`,`indel`)
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



# Dump of table intervals
# ------------------------------------------------------------

CREATE TABLE `intervals` (
  `id` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `chrom` char(2) NOT NULL DEFAULT '',
  `start` int(11) NOT NULL,
  `end` int(11) NOT NULL,
  `finished` tinyint(1) NOT NULL DEFAULT '0',
  PRIMARY KEY (`id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8;



# Dump of table intervals_snps
# ------------------------------------------------------------

CREATE TABLE `intervals_snps` (
  `id` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `interval_id` int(11) unsigned NOT NULL,
  `modsnp_id` int(11) unsigned NOT NULL,
  PRIMARY KEY (`id`),
  UNIQUE KEY `interval_id` (`interval_id`,`modsnp_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8;



# Dump of table modern_snps
# ------------------------------------------------------------

CREATE TABLE `modern_snps` (
  `id` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `population` char(4) DEFAULT NULL,
  `chrom` char(2) NOT NULL DEFAULT '',
  `site` int(11) NOT NULL,
  `ancestral` char(1) NOT NULL DEFAULT '',
  `ancestral_count` int(11) NOT NULL,
  `derived` char(1) NOT NULL DEFAULT '',
  `derived_count` int(11) NOT NULL,
  `type` char(2) DEFAULT NULL,
  `daf` float DEFAULT NULL,
  `variant_id` int(11) DEFAULT NULL,
  `snpchip_id` int(11) DEFAULT NULL,
  `gene_id` int(11) unsigned DEFAULT NULL,
  `neutral` tinyint(1) DEFAULT NULL,
  KEY `id` (`id`),
  KEY `population` (`population`,`chrom`,`site`),
  KEY `variant_id` (`variant_id`),
  KEY `snpchip_id` (`snpchip_id`),
  KEY `gene_id` (`gene_id`),
  KEY `neutral` (`neutral`),
  KEY `chrom` (`chrom`,`site`)
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
  UNIQUE KEY `qtl_id` (`qtl_id`,`modsnp_id`),
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
  KEY `trait_fk` (`trait_id`),
  KEY `pubmed_fk` (`pubmed_id`),
  KEY `peak` (`peak`),
  KEY `species` (`chrom`,`start`,`end`),
  KEY `valid` (`valid`),
  KEY `associationType` (`associationType`,`valid`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8;



# Dump of table sample_bins
# ------------------------------------------------------------

CREATE TABLE `sample_bins` (
  `id` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `sample_id` int(11) unsigned DEFAULT NULL,
  `bin` varchar(100) DEFAULT NULL,
  `overlap` int(11) DEFAULT NULL,
  `perct_overlap` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  UNIQUE KEY `sample_id` (`sample_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8;



# Dump of table sample_dates
# ------------------------------------------------------------

CREATE TABLE `sample_dates` (
  `id` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `age` varchar(255) DEFAULT '',
  `confident` varchar(100) DEFAULT '',
  `lower` int(11) DEFAULT NULL,
  `upper` int(11) DEFAULT NULL,
  `median` int(11) DEFAULT NULL,
  PRIMARY KEY (`id`),
  UNIQUE KEY `age` (`age`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8;



# Dump of table sample_dates_c14
# ------------------------------------------------------------

CREATE TABLE `sample_dates_c14` (
  `id` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `accession` varchar(255) DEFAULT '',
  `confident` varchar(100) DEFAULT '',
  `lower` int(11) DEFAULT NULL,
  `upper` int(11) DEFAULT NULL,
  `median` int(11) DEFAULT NULL,
  PRIMARY KEY (`id`),
  UNIQUE KEY `accession` (`accession`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8;



# Dump of table sample_files
# ------------------------------------------------------------

CREATE TABLE `sample_files` (
  `id` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `sample_id` int(11) unsigned DEFAULT NULL,
  `path` varchar(255) DEFAULT NULL,
  PRIMARY KEY (`id`),
  UNIQUE KEY `sample_id` (`sample_id`,`path`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8;



# Dump of table sample_quality
# ------------------------------------------------------------

CREATE TABLE `sample_quality` (
  `id` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `accession` varchar(255) DEFAULT NULL,
  `mapped_reads` int(11) DEFAULT NULL,
  `total_reads` int(11) DEFAULT NULL,
  `mapped` float DEFAULT NULL,
  `mapped_q30` float DEFAULT NULL,
  `mt_reads` int(11) DEFAULT NULL,
  `duplicates` float DEFAULT NULL,
  `readlen_mapped` float DEFAULT NULL,
  `sd_readlen_mapped` float DEFAULT NULL,
  `readlen_all` float DEFAULT NULL,
  `sd_readlen_all` float DEFAULT NULL,
  `internal` varchar(255) DEFAULT NULL,
  `external` varchar(255) DEFAULT NULL,
  `pool` varchar(255) DEFAULT NULL,
  `species` varchar(255) DEFAULT NULL,
  `read_file` varchar(255) DEFAULT NULL,
  PRIMARY KEY (`id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8;



# Dump of table sample_reads
# ------------------------------------------------------------

CREATE TABLE `sample_reads` (
  `id` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `interval_id` int(11) unsigned NOT NULL,
  `sample_id` int(11) unsigned NOT NULL,
  `chrom` varchar(2) NOT NULL DEFAULT '',
  `site` int(11) NOT NULL,
  `base` varchar(1) NOT NULL DEFAULT '',
  `dist` int(11) DEFAULT NULL,
  `mapq` int(11) DEFAULT NULL,
  `baseq` int(11) DEFAULT NULL,
  `genoq` int(11) DEFAULT NULL,
  `quality` tinyint(1) DEFAULT NULL,
  `called` tinyint(1) DEFAULT NULL,
  KEY `id` (`id`),
  KEY `intervalID` (`interval_id`),
  KEY `sampleID` (`sample_id`),
  KEY `chrom` (`chrom`,`site`),
  KEY `snp` (`called`)
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



# Dump of table samples
# ------------------------------------------------------------

CREATE TABLE `samples` (
  `id` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `accession` varchar(255) NOT NULL DEFAULT '',
  `map_reads` int(11) unsigned DEFAULT NULL,
  `map_prcnt` float DEFAULT NULL,
  `age` varchar(255) DEFAULT NULL,
  `age_int` int(11) DEFAULT NULL,
  `period` varchar(255) DEFAULT NULL,
  `location` varchar(255) DEFAULT NULL,
  `country` varchar(255) DEFAULT NULL,
  `status` varchar(255) DEFAULT NULL,
  `gmm_status` varchar(255) DEFAULT NULL,
  `group` varchar(255) DEFAULT NULL,
  `haplogroup` varchar(255) DEFAULT NULL,
  `mc1r_snp` varchar(3) DEFAULT NULL,
  `valid` tinyint(1) DEFAULT NULL,
  PRIMARY KEY (`id`),
  UNIQUE KEY `accession` (`accession`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8;



# Dump of table snpchip
# ------------------------------------------------------------

CREATE TABLE `snpchip` (
  `id` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `chip_name` varchar(255) NOT NULL DEFAULT '',
  `rsnumber` varchar(255) DEFAULT '',
  `chrom` char(2) NOT NULL DEFAULT '',
  `site` int(11) NOT NULL,
  `snp_name` varchar(255) NOT NULL DEFAULT '',
  PRIMARY KEY (`id`),
  UNIQUE KEY `chip_name` (`chip_name`,`rsnumber`,`snp_name`),
  KEY `snp_name` (`snp_name`),
  KEY `rsnumber` (`rsnumber`),
  KEY `chip_name_2` (`chip_name`,`chrom`,`site`),
  KEY `chrom` (`chrom`,`site`)
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
  `chrom` char(2) NOT NULL DEFAULT '',
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
