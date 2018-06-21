#!/usr/bin/env bash

# call DOM / Przewalski horses
cmd1=$(cat <<EOF
bcftools mpileup --fasta-ref fasta/Equus_caballus.EquCab2.dna.toplevel.fa --bam-list ./data/bam_files_horse_DOM.txt --regions chr{} \
	| bcftools call --multiallelic-caller --ploidy-file ./data/horses.ploidy --samples-file ./data/horses_DOM.samples --output-type v \
	| bcftools view --types snps --exclude INFO/INDEL=1 --output-file ./vcf/horse_DOM_chr{}.vcf &
EOF
)

parallel ${cmd1} ::: {1..31} X

# call DOM2 / Domestic horses
cmd2=$(cat <<EOF
bcftools mpileup --fasta-ref fasta/Equus_caballus.EquCab2.dna.toplevel.fa --bam-list ./data/bam_files_horse_DOM2.txt --regions chr{} \
	| bcftools call --multiallelic-caller --ploidy-file ./data/horses.ploidy --samples-file ./data/horses_DOM2.samples --output-type v \
	| bcftools view --types snps --exclude INFO/INDEL=1 --output-file ./vcf/horse_DOM2_chr{}.vcf &
EOF
)

parallel ${cmd2} ::: {1..31} X

echo "FINISHED calling modern horses!!"