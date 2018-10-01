#!/usr/bin/env bash

# mysqldump w/ gzip
mysqldump -u root -p --single-transaction allele_trajectory \
	| gzip > allele_trajectory-fulldump.sql.gz

# import w/ gzip
gunzip < allele_trajectory-fulldump.sql.gz \
	| mysql -u root -p


#LOAD DATA LOCAL INFILE '/Users/Evan/Downloads/PigsAgeMapping.tsv' INTO TABLE age_mapping IGNORE 1 LINES (age, confident, lower, upper, median);

# SET @i=0;
# SELECT id, name, score, @i:=@i+1 AS rank 
#  FROM ranking 
#  ORDER BY score DESC;

#create table modern_snps_intervals
#select ms.*
#from intervals i
#join modern_snps ms
#  on ms.species = i.species
# and ms.chrom = i.chrom
# and ms.site between i.start and i.end
#where i.species = 'pig';


# run the SNP discovery
nohup python -u build_database.py &> nohup-step1.out &


chrom, quality, random, snp
chrom, baseq, mapq, dist
chrom, site, sample_id


parallel "samtools faidx {} 10:1116-1116" ::: /media/jbod/raid1-sdc1/laurent/full_run_results/Pig/modern/FASTA/*/10.fa > tmp.log; grep -v '>' tmp.log | sort | uniq -c;


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

## get the best SNPs for this locus
#run_cmd(["printf '{locus}' "
#		 "| bedtools intersect -a {snps_file} -b stdin "
#		 "| sort -k 5,5 -g "
#		 "| head -n 1".format(locus=locus, snps_file=snps_file)])
		 
printf '9\t150219744\t150267729\n' | bedtools intersect -a data/sweep/EUD_Sweep_p001_FINAL_cutoff.bed -b stdin | sort -k 5,5 -g | head -n 5
		 
		 
mysqldump -u root -p allele_trajectory ensembl_variants ensembl_genes | gzip > ensembl.sql.gz

# restore dump
gunzip < alleletraj_pig-20180608-1555.sql.gz  | mysql -u root -p -D alleletraj_pig

gunzip < alleletraj_horse-20180623-1004.sql.gz  | mysql -u root -p -D alleletraj_horse


bcftools mpileup --fasta-ref fasta/Sus_scrofa.Sscrofa10.2.dna.toplevel.fa /home/ludo/inbox/BAMslices/Modern/*.bam \
	| bcftools call --multiallelic-caller --output-type v \
	| bcftools view --exclude-types indels,bnd,other --exclude INFO/INDEL=1 --output-file modern_horses.vcf


mkfifo --mode=0666 /tmp/SNPchimp_pig
gzip --stdout -d data/SNPchimp/SNPchimp_pig.tsv.gz > /tmp/SNPchimp_pig


# horses
nohup parallel "~/sratoolkit.2.8.2-1-ubuntu64/bin/fastq-dump --gzip --split-3 --outdir ./fastq {}"  ::: SRR1046129 SRR1769893 SRR1048526 ERR979791 ERR1397960 SRR1046147  SRR2142311 SRR488214 SRR516118 SRR1275408 ERR982704 SRR3017573 ERR979796 ERR982705 ERR982712 ERR982715 ERR982708 ERR982716 ERR982717 ERR979800 ERR868004 SRR1046135 ERR979802 SRR3017575 SRR515216 ERR1021816 ERR1021823 ERR1021824  &> nohup-fastq.out &

# SAMEA2802531 / Equus asinus somalicus
nohup parallel "~/sratoolkit.2.8.2-1-ubuntu64/bin/fastq-dump --gzip --split-3 --outdir ./fastq {}"  ::: ERR650540 ERR650541 ERR650542 ERR650543 ERR650544 ERR650545 ERR650546 ERR650547 ERR650570 ERR650571 ERR650572 ERR650573 ERR650574 ERR650575 ERR650576 ERR650577 ERR650578 ERR650579 ERR650580 ERR650581 ERR650582 ERR650583 ERR650584 ERR650585 ERR650586 ERR650587 ERR650588 ERR650589 ERR650590 ERR650591 ERR650592 ERR650593 ERR650594 ERR650595 ERR650596 ERR650597 ERR650598 ERR650599 ERR650600 ERR650601 ERR650602 ERR650603 ERR650604 ERR650605 ERR650606 ERR650607 ERR650608 ERR650609 ERR650610 ERR650611 ERR650612 ERR650613 ERR650614 ERR650615 ERR650616 ERR650617 ERR650618 ERR650619 ERR650620 ERR650621 ERR650622 ERR650623 ERR650624 ERR650625 ERR650626 ERR650627 ERR650628 ERR650629 ERR650630 ERR650631 ERR650632 ERR650633 ERR650634 ERR650635 ERR650636 ERR650637 ERR650638 ERR650639 ERR650640 ERR650641 ERR650642 ERR650643 ERR650644 ERR650645 ERR650646 ERR650647 ERR650648 ERR650649 ERR650650 ERR650651 ERR650652 ERR650653 ERR650654 ERR650655 ERR650656 ERR650657 ERR650658 ERR650659 ERR650660 ERR650661 ERR650662 ERR650663 ERR650664 ERR650665 ERR650666 ERR650667 ERR650668 ERR650669 ERR650670 ERR650671 ERR650672 ERR650673 ERR650674 ERR650675 ERR650676 ERR650677 ERR650678 ERR650679 ERR650680 ERR650681 ERR650682 ERR650683 ERR650684 ERR650685 ERR650686 ERR650687 ERR650688 ERR650689 ERR650690 ERR650691 ERR650692 ERR650693 ERR650694 ERR650695 ERR650696 ERR650697 ERR650698 ERR650699 ERR650700 ERR650701 ERR650702 ERR650703   &> nohup-somali.out &


# covert PDF to PNG
sips -s format png pdf/all-snps-ages.pdf --out pdf/all-snps-ages.png

# add "chr" to front of chrom names in reference fasta
sed -i -e 's/^>/>chr/' fasta/Equus_caballus.EquCab2.dna.toplevel.fa


# ----------------------------------------------------------------------------------------
# -					sort and index the BAM files that need it							 -
# ----------------------------------------------------------------------------------------

# move files from list into subfolder
cat unsorted.list | xargs mv -t unsorted/

# sort them all, and put in parent folder
parallel "echo {} && samtools sort -O bam {} > ../{}" ::: *.bam

# check the filesize of the sorted files
ls -lh | grep -f unsorted.list

# index them all
parallel -a unsorted.list "echo {} && samtools index {}"

# ----------------------------------------------------------------------------------------

# get the rsnumbers for the pig mummy markers
grep -E 'WU_10.2_18_43089973|ASGA0079708|ALGA0098112|MARC0069672|MARC0089391|WU_10.2_18_43161644|WU_10.2_18_43204315|ASGA0097792|MARC0046857|WU_10.2_18_43287167|H3GA0050799|ASGA0079716|WU_10.2_18_43418494|WU_10.2_18_43464036|ASGA0079719|ASGA0079723|ASGA0079726|ASGA0079728|WU_10.2_18_43630319|WU_10.2_18_43773633|WU_10.2_18_43778188|WU_10.2_18_43803484|ASGA0079732|H3GA0050810|ASGA0079734|ALGA0098146' ~/Downloads/SNPchimp_result_3294637283.tsv

nohup python -u build_database.py &> nohup-pigs.log &
tail -f nohup-pigs.log

