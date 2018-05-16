#!/usr/bin/env bash

rsync -avz --partial \
           --exclude '.git' \
           --exclude '.gitignore' \
           --exclude 'nohup*' \
           --exclude '*.sh' \
           --exclude '*.py*' \
           --exclude '*.r' \
           --exclude '*.R' \
           --exclude '*.old*' \
           --exclude '*tmp*' \
           --exclude '*.sql.gz' \
           --exclude '*.ld' \
           --exclude 'screening' \
           --exclude 'vcf' \
           evan@palaeoprime:~/alleletraj/ ~/Dropbox/Code/alleletraj/
