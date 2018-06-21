#!/usr/bin/env bash

# call DOM / Przewalski horses
bcftools mpileup --fasta-ref fasta/Sus_scrofa.Sscrofa10.2.dna.toplevel.fa --bam-list ./data/bam_files_horse_DOM.txt \
	| bcftools call --multiallelic-caller --output-type v \
	| bcftools view --exclude-types indels,bnd,other --exclude INFO/INDEL=1 --output-file modern_DOM.vcf &

# call DOM2 / Domestic horses
bcftools mpileup --fasta-ref fasta/Sus_scrofa.Sscrofa10.2.dna.toplevel.fa --bam-list ./data/bam_files_horse_DOM2.txt \
	| bcftools call --multiallelic-caller --output-type v \
	| bcftools view --exclude-types indels,bnd,other --exclude INFO/INDEL=1 --output-file modern_DOM2.vcf &
