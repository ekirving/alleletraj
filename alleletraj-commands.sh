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


bcftools mpileup --region {region} --targets-file {targets} --fasta-ref {ref} {bams} \
	| bcftools call --multiallelic-caller --output-type v \
	| bcftools view --exclude-types indels,bnd,other --exclude INFO/INDEL=1 --output-file {vcf}

bcftools mpileup --region 1:236933-336933 --targets-file vcf/diploid-int1-sample147.tsv.gz --fasta-ref fasta/Sus_scrofa.Sscrofa10.2.dna.toplevel.fa /media/jbod/raid5-sdb/laurent/full_run_results/Pig/KD024/KD024_rmdup.bam /media/jbod/raid5-sdb/laurent/full_run_results/Pig/mtDNA/bam/KD024/KD024_rmdup.bam  \
	| bcftools call --multiallelic-caller --targets-file vcf/diploid-int1-sample147.tsv.gz --constrain alleles --output-type v  \
	| bcftools view --exclude-types indels,bnd,other --exclude INFO/INDEL=1 --output-file vcf/diploid-int1-sample147.vcf

# get the best SNPs for this locus
run_cmd(["printf '{locus}' "
		 "| bedtools intersect -a {snps_file} -b stdin "
		 "| sort -k 5,5 -g "
		 "| head -n 1".format(locus=locus, snps_file=snps_file)])
		 
printf '9\t150219744\t150267729\n' | bedtools intersect -a data/sweep/EUD_Sweep_p001_FINAL_cutoff.bed -b stdin | sort -k 5,5 -g | head -n 5
		 
		 
mysqldump -u root -p allele_trajectory ensembl_variants ensembl_genes | gzip > ensembl.sql.gz




mkfifo --mode=0666 /tmp/SNPchimp_pig

gzip --stdout -d data/SNPchimp/SNPchimp_pig.tsv.gz > /tmp/SNPchimp_pig
