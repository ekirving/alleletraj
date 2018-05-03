#!/usr/bin/env bash

# bwa aln -l 1024 -n 0.03 -t 20 /home/laurent/REF/SUS/Sus_scrofa.Sscrofa10.2.dna.toplevel.fa AA133_merged_adRm.fq.gz
# bwa aln -l 1024 -n 0.03 -t 20 /home/laurent/REF/SUS/Sus_scrofa.Sscrofa10.2.dna.toplevel.fa AA241_merged_adRm.fq.gz
# bwa aln -l 1024 -n 0.03 -t 20 /home/laurent/REF/SUS/Sus_scrofa.Sscrofa10.2.dna.toplevel.fa AA325_merged_adRm.fq.gz
# bwa aln -l 1024 -n 0.03 -t 20 /home/laurent/REF/SUS/Sus_scrofa.Sscrofa10.2.dna.toplevel.fa AA363_merged_adRm.fq.gz
# bwa aln -l 1024 -n 0.03 -t 20 /home/laurent/REF/SUS/Sus_scrofa.Sscrofa10.2.dna.toplevel.fa AL712_merged_adRm.fq.gz

bwa mem -t 10 -R '@RG\tID:AA133_merged\\tSM:AA133_merged' /home/laurent/REF/SUS/Sus_scrofa.Sscrofa10.2.dna.toplevel.fa media/jbod/raid5-sdb/laurent/screen_results_14062016/raw_split/pig_WG/AA133_merged_adRm.fq.gz \
	| /usr/local/bin/samtools1.3 sort -@ 10 -O bam -o ~/alleletraj/screening/AA133_merged.bam - &

bwa mem -t 10 -R '@RG\tID:AA241_merged\\tSM:AA241_merged' /home/laurent/REF/SUS/Sus_scrofa.Sscrofa10.2.dna.toplevel.fa media/jbod/raid5-sdb/laurent/screen_results_14062016/raw_split/pig_WG/AA241_merged_adRm.fq.gz \
	| /usr/local/bin/samtools1.3 sort -@ 10 -O bam -o ~/alleletraj/screening/AA241_merged.bam - &
	
bwa mem -t 10 -R '@RG\tID:AA325_merged\\tSM:AA325_merged' /home/laurent/REF/SUS/Sus_scrofa.Sscrofa10.2.dna.toplevel.fa media/jbod/raid5-sdb/laurent/screen_results_14062016/raw_split/pig_WG/AA325_merged_adRm.fq.gz \
	| /usr/local/bin/samtools1.3 sort -@ 10 -O bam -o ~/alleletraj/screening/AA325_merged.bam - &

bwa mem -t 10 -R '@RG\tID:AA363_merged\\tSM:AA363_merged' /home/laurent/REF/SUS/Sus_scrofa.Sscrofa10.2.dna.toplevel.fa media/jbod/raid5-sdb/laurent/screen_results_14062016/raw_split/pig_WG/AA363_merged_adRm.fq.gz \
	| /usr/local/bin/samtools1.3 sort -@ 10 -O bam -o ~/alleletraj/screening/AA363_merged.bam - &
	
bwa mem -t 10 -R '@RG\tID:AL712_merged\\tSM:AL712_merged' /home/laurent/REF/SUS/Sus_scrofa.Sscrofa10.2.dna.toplevel.fa media/jbod/raid5-sdb/laurent/screen_results_14062016/raw_split/pig_WG/AL712_merged_adRm.fq.gz \
	| /usr/local/bin/samtools1.3 sort -@ 10 -O bam -o ~/alleletraj/screening/AL712_merged.bam - &
