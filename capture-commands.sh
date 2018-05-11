#!/usr/bin/env bash

# bwa aln -l 1024 -n 0.03 -t 20 /home/laurent/REF/SUS/Sus_scrofa.Sscrofa10.2.dna.toplevel.fa AA133_merged_adRm.fq.gz
# bwa aln -l 1024 -n 0.03 -t 20 /home/laurent/REF/SUS/Sus_scrofa.Sscrofa10.2.dna.toplevel.fa AA241_merged_adRm.fq.gz
# bwa aln -l 1024 -n 0.03 -t 20 /home/laurent/REF/SUS/Sus_scrofa.Sscrofa10.2.dna.toplevel.fa AA325_merged_adRm.fq.gz
# bwa aln -l 1024 -n 0.03 -t 20 /home/laurent/REF/SUS/Sus_scrofa.Sscrofa10.2.dna.toplevel.fa AA363_merged_adRm.fq.gz
# bwa aln -l 1024 -n 0.03 -t 20 /home/laurent/REF/SUS/Sus_scrofa.Sscrofa10.2.dna.toplevel.fa AL712_merged_adRm.fq.gz

# realign all the capture data
bwa mem -t 10 -R '@RG\tID:AA133_merged\tSM:AA133_merged' /home/laurent/REF/SUS/Sus_scrofa.Sscrofa10.2.dna.toplevel.fa /media/jbod/raid5-sdb/laurent/screen_results_14062016/raw_split/pig_WG/AA133_merged_adRm.fq.gz \
	| /usr/local/bin/samtools1.3 sort -@ 10 -O bam -o ~/alleletraj/screening/AA133_merged.bam - &

bwa mem -t 10 -R '@RG\tID:AA241_merged\tSM:AA241_merged' /home/laurent/REF/SUS/Sus_scrofa.Sscrofa10.2.dna.toplevel.fa /media/jbod/raid5-sdb/laurent/screen_results_14062016/raw_split/pig_WG/AA241_merged_adRm.fq.gz \
	| /usr/local/bin/samtools1.3 sort -@ 10 -O bam -o ~/alleletraj/screening/AA241_merged.bam - &
	
bwa mem -t 10 -R '@RG\tID:AA325_merged\tSM:AA325_merged' /home/laurent/REF/SUS/Sus_scrofa.Sscrofa10.2.dna.toplevel.fa /media/jbod/raid5-sdb/laurent/screen_results_14062016/raw_split/pig_WG/AA325_merged_adRm.fq.gz \
	| /usr/local/bin/samtools1.3 sort -@ 10 -O bam -o ~/alleletraj/screening/AA325_merged.bam - &

bwa mem -t 10 -R '@RG\tID:AA363_merged\tSM:AA363_merged' /home/laurent/REF/SUS/Sus_scrofa.Sscrofa10.2.dna.toplevel.fa /media/jbod/raid5-sdb/laurent/screen_results_14062016/raw_split/pig_WG/AA363_merged_adRm.fq.gz \
	| /usr/local/bin/samtools1.3 sort -@ 10 -O bam -o ~/alleletraj/screening/AA363_merged.bam - &
	
bwa mem -t 10 -R '@RG\tID:AL712_merged\tSM:AL712_merged' /home/laurent/REF/SUS/Sus_scrofa.Sscrofa10.2.dna.toplevel.fa /media/jbod/raid5-sdb/laurent/screen_results_14062016/raw_split/pig_WG/AL712_merged_adRm.fq.gz \
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

# call consensus on PCR duplicates (from SAM format)
/usr/local/bin/samtools1.3 view -h screening/AA133_merged.bam | ./FilterUniqueSAMCons.py | /usr/local/bin/samtools1.3 view -b > screening/AA133_merged_rmdup.bam &
/usr/local/bin/samtools1.3 view -h screening/AA325_merged.bam | ./FilterUniqueSAMCons.py | /usr/local/bin/samtools1.3 view -b > screening/AA325_merged_rmdup.bam &
/usr/local/bin/samtools1.3 view -h screening/AA241_merged.bam | ./FilterUniqueSAMCons.py | /usr/local/bin/samtools1.3 view -b > screening/AA241_merged_rmdup.bam &
/usr/local/bin/samtools1.3 view -h screening/AA363_merged.bam | ./FilterUniqueSAMCons.py | /usr/local/bin/samtools1.3 view -b > screening/AA363_merged_rmdup.bam &
/usr/local/bin/samtools1.3 view -h screening/AL712_merged.bam | ./FilterUniqueSAMCons.py | /usr/local/bin/samtools1.3 view -b > screening/AL712_merged_rmdup.bam &

# bug fix
for accession in AA133 AA325 AA241 AA363 AL712; do
  /usr/local/bin/samtools1.3 view -h screening/"$accession"_merged.bam | sed 's/\\t/\t/' | /usr/local/bin/samtools1.3 view -b > screening/"$accession"_merged_fixed.bam &
done

# do a conventional deduplicate
for accession in AA133 AA325 AA241 AA363 AL712; do
  java -Xmx8G -jar /usr/local/picard-tools-2.5.0/picard.jar MarkDuplicates INPUT=screening/"$accession"_merged_fixed.bam OUTPUT=screening/"$accession"_merged_rmdup_picard.bam METRICS_FILE=screening/"$accession"_merged_rmdup_picard.log REMOVE_DUPLICATES=true QUIET=true &
done

# count post depuplicate
parallel "echo {} && samtools view -F 0x904 -c {}" ::: ~/alleletraj/screening/*_merged_rmdup_picard.bam

# AA133_merged_rmdup_picard.bam	179
# AA325_merged_rmdup_picard.bam	173404
# AA241_merged_rmdup_picard.bam	249819
# AA363_merged_rmdup_picard.bam	9536188
# AL712_merged_rmdup_picard.bam	5690136

# count unaligned reads
parallel "echo {} && samtools view -f 0x4 -c {}" ::: ~/alleletraj/screening/*_merged.bam
AA133_merged.bam	8282
AA325_merged.bam	90969
AA241_merged.bam	7183510
AA363_merged.bam	3934666
AL712_merged.bam	31266732

# count unaligned reads after depuplicating
parallel "echo {} && samtools view -f 0x04 -c {}" ::: ~/alleletraj/screening/*_merged_rmdup.bam




# now get the good reads above mapq30
parallel "echo {} && samtools view -F 0x904 -q 30 -c {}" ::: ~/alleletraj/screening/*_merged_rmdup_picard.bam

# AA133_merged_rmdup_picard.bam	24
# AA325_merged_rmdup_picard.bam	97817
# AA241_merged_rmdup_picard.bam	19741
# AA363_merged_rmdup_picard.bam	2533926
# AL712_merged_rmdup_picard.bam	1499640

# index the BAM files
# parallel "samtools index {}" ::: ~/alleletraj/screening/*.bam

parallel "echo {} && samtools view -F 0x904 -c {}" ::: ~/alleletraj/screening/*_merged_rmdup.bam

# AA133_merged_rmdup.bam	178
# AA325_merged_rmdup.bam	172666
# AA241_merged_rmdup.bam	252129
# AA363_merged_rmdup.bam	9516381

parallel "echo {} && samtools view -F 0x904 -q 30 -c {}" ::: ~/alleletraj/screening/*_merged_rmdup.bam

# AA133_merged_rmdup.bam	24
# AA325_merged_rmdup.bam	97824
# AA241_merged_rmdup.bam	20465
# AA363_merged_rmdup.bam	2542448


parallel "echo {} && qualimap bamqc -bam {} -outdir qualimap_results -outformat pdf" ::: ~/alleletraj/screening/*.bam

for accession in AA133 AA325 AA241 AA363 AL712; do
  qualimap bamqc -bam screening/"$accession"_merged_rmdup_picard.bam -outfile qualimap_results/"$accession"_merged_rmdup_picard.pdf &
done

qualimap bamqc -bam screening/AA363_merged_rmdup_picard.bam -outfile qualimap_results/AA363_merged_rmdup_picard.pdf