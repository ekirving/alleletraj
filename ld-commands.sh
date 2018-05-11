#!/usr/bin/env bash

# https://www.biostars.org/p/300381/

# extract EUD population
/usr/local/plink-v1.90b4.4/plink --make-bed --keep-fam eud.list --bfile /media/jbod/raid1-sdc1/laurent/full_run_results/Pig/ped/FINAL/merged_Auto_HighCov --out eud

# thin the SNPs
mapthin -b 20 eud.bim eud-thin.bim

# make the bed file
plink --bfile eud --make-bed --extract 'eud-thin.bim' --out 'eud-thin'

# compute the LD
plink --bfile "eud-thin" --r2 --ld-window-r2 0 --ld-window 999999 --ld-window-kb 8000 --out "eud-thin"


# plink --bfile "eud" --r2 --ld-window-r2 0 --ld-window 999999 --ld-window-kb 8000 --out "eud"
