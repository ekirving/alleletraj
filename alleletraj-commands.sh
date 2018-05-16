#!/usr/bin/env bash

# mysqldump w/ gzip
mysqldump -u root -p --single-transaction allele_trajectory \
	| gzip > allele_trajectory-fulldump.sql.gz

# import w/ gzip
gunzip < allele_trajectory-fulldump.sql.gz \
	| mysql -u root -p


LOAD DATA LOCAL INFILE '/Users/Evan/Downloads/PigsAgeMapping.tsv' INTO TABLE age_mapping IGNORE 1 LINES (age, confident, lower, upper, median);


create table modern_snps_intervals
select ms.*
from intervals i
join modern_snps ms 
  on ms.species = i.species
 and ms.chrom = i.chrom
 and ms.site between i.start and i.end
where i.species = 'pig';


# run the SNP discovery
nohup python -u build_database.py &> nohup-step1.out &


chrom, quality, random, snp
chrom, baseq, mapq, dist
chrom, site, sample_id


parallel "samtools faidx {} 10:1116-1116" ::: /media/jbod/raid1-sdc1/laurent/full_run_results/Pig/modern/FASTA/*/10.fa > tmp.log; grep -v '>' tmp.log | sort | uniq -c;


age	
confident	
lower
upper
median

nohup mysql < sample_reads_innodb_part.sql &> nohup_innodb_part.out &
nohup mysql < sample_reads_innodb.sql &> nohup_innodb.out &
nohup mysql < sample_reads_mysiam_part.sql &> nohup_mysiam_part.out &
nohup mysql < sample_reads_mysiam.sql &> nohup_mysiam.out &