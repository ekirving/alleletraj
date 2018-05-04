#!/usr/bin/env bash

# bwa aln -l 1024 -n 0.03 -t 20 /home/laurent/REF/SUS/Sus_scrofa.Sscrofa10.2.dna.toplevel.fa AA133_merged_adRm.fq.gz
# bwa aln -l 1024 -n 0.03 -t 20 /home/laurent/REF/SUS/Sus_scrofa.Sscrofa10.2.dna.toplevel.fa AA241_merged_adRm.fq.gz
# bwa aln -l 1024 -n 0.03 -t 20 /home/laurent/REF/SUS/Sus_scrofa.Sscrofa10.2.dna.toplevel.fa AA325_merged_adRm.fq.gz
# bwa aln -l 1024 -n 0.03 -t 20 /home/laurent/REF/SUS/Sus_scrofa.Sscrofa10.2.dna.toplevel.fa AA363_merged_adRm.fq.gz
# bwa aln -l 1024 -n 0.03 -t 20 /home/laurent/REF/SUS/Sus_scrofa.Sscrofa10.2.dna.toplevel.fa AL712_merged_adRm.fq.gz

# realign all the capture data
bwa mem -t 10 -R '@RG\tID:AA133_merged\\tSM:AA133_merged' /home/laurent/REF/SUS/Sus_scrofa.Sscrofa10.2.dna.toplevel.fa /media/jbod/raid5-sdb/laurent/screen_results_14062016/raw_split/pig_WG/AA133_merged_adRm.fq.gz \
	| /usr/local/bin/samtools1.3 sort -@ 10 -O bam -o ~/alleletraj/screening/AA133_merged.bam - &

bwa mem -t 10 -R '@RG\tID:AA241_merged\\tSM:AA241_merged' /home/laurent/REF/SUS/Sus_scrofa.Sscrofa10.2.dna.toplevel.fa /media/jbod/raid5-sdb/laurent/screen_results_14062016/raw_split/pig_WG/AA241_merged_adRm.fq.gz \
	| /usr/local/bin/samtools1.3 sort -@ 10 -O bam -o ~/alleletraj/screening/AA241_merged.bam - &
	
bwa mem -t 10 -R '@RG\tID:AA325_merged\\tSM:AA325_merged' /home/laurent/REF/SUS/Sus_scrofa.Sscrofa10.2.dna.toplevel.fa /media/jbod/raid5-sdb/laurent/screen_results_14062016/raw_split/pig_WG/AA325_merged_adRm.fq.gz \
	| /usr/local/bin/samtools1.3 sort -@ 10 -O bam -o ~/alleletraj/screening/AA325_merged.bam - &

bwa mem -t 10 -R '@RG\tID:AA363_merged\\tSM:AA363_merged' /home/laurent/REF/SUS/Sus_scrofa.Sscrofa10.2.dna.toplevel.fa /media/jbod/raid5-sdb/laurent/screen_results_14062016/raw_split/pig_WG/AA363_merged_adRm.fq.gz \
	| /usr/local/bin/samtools1.3 sort -@ 10 -O bam -o ~/alleletraj/screening/AA363_merged.bam - &
	
bwa mem -t 10 -R '@RG\tID:AL712_merged\\tSM:AL712_merged' /home/laurent/REF/SUS/Sus_scrofa.Sscrofa10.2.dna.toplevel.fa /media/jbod/raid5-sdb/laurent/screen_results_14062016/raw_split/pig_WG/AL712_merged_adRm.fq.gz \
	| /usr/local/bin/samtools1.3 sort -@ 10 -O bam -o ~/alleletraj/screening/AL712_merged.bam - &

# count number of sequences in fastq files
parallel "echo {} && gunzip -c {} | wc -l | awk '{d=\$1; print d/4;}'" ::: /media/jbod/raid5-sdb/laurent/screen_results_14062016/raw_split/pig_WG/*_merged_adRm.fq.gz

# AA133_merged_adRm.fq.gz		    8499
# AA325_merged_adRm.fq.gz		  271797
# AA241_merged_adRm.fq.gz		12623655
# AA363_merged_adRm.fq.gz		39118759
# AL712_merged_adRm.fq.gz		86680674

# count the number of reads which successfully aligned
# exclude:
#	0x4		UNMAP			segment unmapped
#	0x100	SECONDARY		secondary alignment
#	0x800	SUPPLEMENTARY	supplementary alignment
parallel "echo {} && samtools view -F 0x904 -c {}" ::: ~/alleletraj/screening/*.bam

# AA133_merged.bam					 217
# AA325_merged.bam				  180828
# AA241_merged.bam				 5440145
# AA363_merged.bam				35184093
# AL712_merged.bam				55413942

# count those above mapq30
parallel "echo {} && samtools view -F 0x904 -q 30 -c {}" ::: ~/alleletraj/screening/*.bam

# AA133_merged.bam	24
# AA325_merged.bam	98648
# AA241_merged.bam	564814
# AA363_merged.bam	6172412
# AL712_merged.bam	8645762