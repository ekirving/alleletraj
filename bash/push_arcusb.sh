#!/usr/bin/env bash

rsync -avz --partial --relative \
    ~/alleletraj/./bash \
    ~/alleletraj/./data/selection \
    ~/alleletraj/./slurm-jobs.list \
    scro2860@arcus-b.arc.ox.ac.uk:~/alleletraj/
