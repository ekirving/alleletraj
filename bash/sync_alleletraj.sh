#!/usr/bin/env bash

rsync -avz --partial \
           --exclude '.git*' \
           --exclude '*.sh' \
           --exclude '*.txt' \
           --exclude '*.py*' \
           --exclude '*.r' \
           --exclude '*.R' \
           --exclude '*.old*' \
           --exclude '*tmp*' \
           --exclude '*temp*' \
           --exclude '*.ld' \
           --exclude '*.csv' \
           --exclude '*.xlsx' \
           --exclude '*.haplo.gz' \
           --exclude '*.fa' \
           --exclude '*.tped' \
           --exclude '*.vcf.gz' \
           --exclude '*.simple' \
           --exclude '*epoch*/*' \
           --exclude '*.sql' \
           --exclude '*.bak*' \
           --exclude 'data/fastq' \
           --exclude 'data/bam' \
           --exclude 'data/vcf' \
           --exclude '*.gz' \
           --exclude '*.tgz' \
           --exclude 'data/selection/*.log' \
           --exclude 'data/selection/*.traj*' \
           --exclude 'data/selection/*.time*' \
           --exclude 'data/selection/*.param' \
           --exclude 'data/selection.partial/' \
           --exclude 'data/selection.truncated/' \
           --exclude '*SelectionBenchmark*' \
           --exclude 'cattle-Dom-*.input' \
           --exclude 'data/db/' \
           evan@palaeoprime:~/alleletraj/ ~/Dropbox/Code/alleletraj/
