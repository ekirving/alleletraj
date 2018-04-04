#!/usr/bin/env bash

# mysqldump w/ gzip
mysqldump -u root -p --single-transaction allele_trajectory \
	| gzip > allele_trajectory-fulldump.sql.gz

# import w/ gzip
gunzip < allele_trajectory-fulldump.sql.gz \
	| mysql -u root -p

	
gunzip -c ./cow_9913/BED/*.bed.gz | awk '{if(NR>1) print $4"\t"substr($1,4)"\t"$3}' > dbsnp.cow_9913.data &
gunzip -c ./goat_9925/BED/*.bed.gz | awk '{if(NR>1) print $4"\t"substr($1,4)"\t"$3}' > dbsnp.goat_9925.data &
gunzip -c ./horse_9796/BED/*.bed.gz | awk '{if(NR>1) print $4"\t"substr($1,4)"\t"$3}' > dbsnp.horse_9796.data &
gunzip -c ./pig_9823/BED/*.bed.gz | awk '{if(NR>1) print $4"\t"substr($1,4)"\t"$3}' > dbsnp.pig_9823.data &


LOAD DATA LOCAL INFILE '/home/evan/dnSNP/dbsnp.cow_9913.data' INTO TABLE dbsnp (rsnumber, chrom, site) SET species = 'cow';
LOAD DATA LOCAL INFILE '/home/evan/dnSNP/dbsnp.goat_9925.data' INTO TABLE dbsnp (rsnumber, chrom, site) SET species = 'goat';
LOAD DATA LOCAL INFILE '/home/evan/dnSNP/dbsnp.horse_9796.data' INTO TABLE dbsnp (rsnumber, chrom, site) SET species = 'horse';
LOAD DATA LOCAL INFILE '/home/evan/dnSNP/dbsnp.pig_9823.data' INTO TABLE dbsnp (rsnumber, chrom, site) SET species = 'pig';


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