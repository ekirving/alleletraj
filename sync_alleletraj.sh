#!/usr/bin/env bash

rsync -avz --partial \
           --progress \
           --stats \
           --exclude '.git' \
           --exclude '.gitignore' \
           --exclude 'nohup*' \
           --exclude '*.sh' \
           --exclude '*.txt' \
           --exclude '*.py*' \
           --exclude '*.r' \
           --exclude '*.R' \
           --exclude '*.old*' \
           --exclude '*tmp*' \
           --exclude '*.ld' \
           --exclude 'screening' \
           --exclude 'vcf' \
           evan@palaeoprime:~/alleletraj/ ~/Dropbox/Code/alleletraj/
