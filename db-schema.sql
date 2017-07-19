# ************************************************************
# Sequel Pro SQL dump
# Version 4541
#
# http://www.sequelpro.com/
# https://github.com/sequelpro/sequelpro
#
# Host: 127.0.0.1 (MySQL 5.5.55-0ubuntu0.14.04.1)
# Database: allele_trajectory
# Generation Time: 2017-07-19 22:05:50 +0000
# ************************************************************


/*!40101 SET @OLD_CHARACTER_SET_CLIENT=@@CHARACTER_SET_CLIENT */;
/*!40101 SET @OLD_CHARACTER_SET_RESULTS=@@CHARACTER_SET_RESULTS */;
/*!40101 SET @OLD_COLLATION_CONNECTION=@@COLLATION_CONNECTION */;
/*!40101 SET NAMES utf8 */;
/*!40014 SET @OLD_FOREIGN_KEY_CHECKS=@@FOREIGN_KEY_CHECKS, FOREIGN_KEY_CHECKS=0 */;
/*!40101 SET @OLD_SQL_MODE=@@SQL_MODE, SQL_MODE='NO_AUTO_VALUE_ON_ZERO' */;
/*!40111 SET @OLD_SQL_NOTES=@@SQL_NOTES, SQL_NOTES=0 */;


# Dump of table pubmeds
# ------------------------------------------------------------

DROP TABLE IF EXISTS `pubmeds`;

CREATE TABLE `pubmeds` (
  `id` int(11) unsigned NOT NULL,
  `authors` text,
  `year` char(255) DEFAULT NULL,
  `title` text,
  `journal` char(255) DEFAULT NULL,
  PRIMARY KEY (`id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8;



# Dump of table qtls
# ------------------------------------------------------------

DROP TABLE IF EXISTS `qtls`;

CREATE TABLE `qtls` (
  `id` int(11) unsigned NOT NULL,
  `associationType` char(255) DEFAULT NULL,
  `symbol` char(255) DEFAULT NULL,
  `pubmedID` int(11) unsigned DEFAULT NULL,
  `traitID` int(11) unsigned NOT NULL,
  `chromosome` char(2) DEFAULT NULL,
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
  `geneID` char(255) DEFAULT NULL,
  `geneSource` char(255) DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `trait_fk` (`traitID`),
  KEY `pubmed_fk` (`pubmedID`),
  KEY `chromosome` (`chromosome`,`genomeLoc_start`,`genomeLoc_end`),
  CONSTRAINT `pubmed_fk` FOREIGN KEY (`pubmedID`) REFERENCES `pubmeds` (`id`),
  CONSTRAINT `trait_fk` FOREIGN KEY (`traitID`) REFERENCES `traits` (`id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8;



# Dump of table sample_reads
# ------------------------------------------------------------

DROP TABLE IF EXISTS `sample_reads`;

CREATE TABLE `sample_reads` (
  `id` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `sampleID` int(11) unsigned NOT NULL,
  `chrom` varchar(2) NOT NULL DEFAULT '',
  `pos` int(11) NOT NULL,
  `base` varchar(1) NOT NULL DEFAULT '',
  `mapq` int(11) NOT NULL,
  `baseq` int(11) NOT NULL,
  `dist` int(11) NOT NULL,
  `quality` tinyint(1) DEFAULT NULL,
  `random` tinyint(1) DEFAULT NULL,
  `snp` tinyint(1) DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `snp` (`snp`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8;



# Dump of table samples
# ------------------------------------------------------------

DROP TABLE IF EXISTS `samples`;

CREATE TABLE `samples` (
  `id` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `accession` varchar(255) DEFAULT NULL,
  `map_reads` int(11) DEFAULT NULL,
  `map_prcnt` float DEFAULT NULL,
  `map_prcnt_q30` float DEFAULT NULL,
  `age` varchar(255) DEFAULT NULL,
  `location` varchar(255) DEFAULT NULL,
  `country` varchar(255) DEFAULT NULL,
  `status` varchar(255) DEFAULT NULL,
  `path` varchar(255) DEFAULT NULL,
  PRIMARY KEY (`id`),
  UNIQUE KEY `accession` (`accession`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8;



# Dump of table traits
# ------------------------------------------------------------

DROP TABLE IF EXISTS `traits`;

CREATE TABLE `traits` (
  `id` int(11) unsigned NOT NULL,
  `class` char(255) DEFAULT NULL,
  `type` char(255) DEFAULT NULL,
  `name` char(255) DEFAULT NULL,
  `VT` char(255) DEFAULT NULL,
  `CMO` char(255) DEFAULT NULL,
  `LPT` char(255) DEFAULT NULL,
  PRIMARY KEY (`id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8;




/*!40111 SET SQL_NOTES=@OLD_SQL_NOTES */;
/*!40101 SET SQL_MODE=@OLD_SQL_MODE */;
/*!40014 SET FOREIGN_KEY_CHECKS=@OLD_FOREIGN_KEY_CHECKS */;
/*!40101 SET CHARACTER_SET_CLIENT=@OLD_CHARACTER_SET_CLIENT */;
/*!40101 SET CHARACTER_SET_RESULTS=@OLD_CHARACTER_SET_RESULTS */;
/*!40101 SET COLLATION_CONNECTION=@OLD_COLLATION_CONNECTION */;
