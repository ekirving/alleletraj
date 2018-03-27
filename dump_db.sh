#!/usr/bin/env bash

now=`date '+%Y%m%d-%H%M'`;
filename="allele_trajectory-$now.sql.gz"

mysqldump allele_trajectory \
    --ignore-table=allele_trajectory.dbsnp_complete \
    --ignore-table=allele_trajectory.dbsnp_snpchip \
    | gzip > $filename