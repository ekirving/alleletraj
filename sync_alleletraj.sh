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
           --exclude '*tmp*' \
           --exclude '*.ld' \
           --exclude 'screening' \
           --exclude 'fastq' \
           --exclude 'vcf' \
           --exclude 'selection/*.traj' \
           --exclude 'selection/*.time' \
           evan@palaeoprime:~/alleletraj/ ~/Dropbox/Code/alleletraj/
