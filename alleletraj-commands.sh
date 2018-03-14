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




# run the SNP discovery
nohup python -u build_database.py &> nohup-step1.out &


chrom, quality, random, snp
chrom, baseq, mapq, dist
chrom, pos, sampleID


parallel "samtools faidx {} 10:1116-1116" ::: /media/jbod/raid1-sdc1/laurent/full_run_results/Pig/modern/FASTA/*/10.fa > tmp.log; grep -v '>' tmp.log | sort | uniq -c;