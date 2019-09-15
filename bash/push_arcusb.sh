#!/usr/bin/env bash

rsync -avz --partial --relative \
    ~/alleletraj/./bash \
    ~/alleletraj/./data/selection/ \
    ~/alleletraj/./data/dadi/*.pop \
    ~/alleletraj/./slurm-* \
    --exclude '*.gz' \
    --exclude '*.log' \
    --exclude '*.diag' \
    --exclude '*.ess' \
    --exclude '*.map' \
    --exclude '*.psrf' \
    --exclude '*.done' \
    --exclude '*.trunc' \
    --exclude '*.mpsrf' \
    --exclude '*luigi-tmp*' \
    scro2860@arcus-b.arc.ox.ac.uk:~/alleletraj/

rsync -avz --partial --relative \
    ~/./bin/sr \
    scro2860@arcus-b.arc.ox.ac.uk:~/
