#!/usr/bin/env bash

rsync -avz --partial \
           --progress \
           --stats \
           --exclude '.git' \
           --exclude '.gitignore' \
           --exclude '*.sh' \
           --exclude '*.txt' \
           --exclude '*.py*' \
           --exclude '*.r' \
           --exclude '*.R' \
           --exclude '*.old*' \
           --exclude '*.bak' \
           --exclude '*tmp*' \
           --exclude '*.ld' \
           --exclude '*.csv' \
           --exclude 'screening' \
           --exclude 'fastq' \
           --exclude 'bam' \
           --exclude 'vcf' \
           --exclude 'selection/*.traj' \
           --exclude 'selection/*.time' \
           --exclude '*.gz' \
           evan@palaeoprime:~/alleletraj/ ~/Dropbox/Code/alleletraj/
