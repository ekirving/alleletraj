#!/usr/bin/env bash

rsync -avz --partial --relative \
    ~/alleletraj/./bash \
    ~/alleletraj/./data/selection \
    ~/alleletraj/./data/dadi/*.pop \
    ~/alleletraj/./slurm-jobs.list \
    scro2860@arcus-b.arc.ox.ac.uk:~/alleletraj/

rsync -avz --partial --relative \
    ~/./bin/sr \
    scro2860@arcus-b.arc.ox.ac.uk:~/
