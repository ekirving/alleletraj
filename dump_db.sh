#!/usr/bin/env bash

now=`date '+%Y%m%d-%H%M'`;
filename="allele_trajectory-$now.sql.gz"

# dump all the tables except `sample_reads`
mysqldump allele_trajectory \
    --ignore-table=allele_trajectory.sample_reads \
    | gzip > $filename

# now just dump the SNPs from `sample_reads` and append to the gzip file
mysqldump allele_trajectory \
    sample_reads --where 'snp = 1' \
    | gzip >> $filename
