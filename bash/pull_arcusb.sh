#!/usr/bin/env bash

rsync -avz --partial \
    scro2860@arcus-b.arc.ox.ac.uk:~/alleletraj/data/selection \
    ~/alleletraj/data/selection
