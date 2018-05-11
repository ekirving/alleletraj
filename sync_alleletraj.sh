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
           evan@palaeoprime:~/alleletraj/ ~/dropbox/code/alleletraj/

