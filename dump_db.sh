#!/usr/bin/env bash

now=`date '+%Y%m%d-%H%M'`;
filename="alleletraj_$1-$now.sql.gz"

# dump all the tables except `sample_reads`
mysqldump alleletraj_"$1" \
    --ignore-table=alleletraj_"$1".sample_reads \
    | gzip > ${filename}

# now just dump the SNPs from `sample_reads` and append to the gzip file
mysqldump alleletraj_"$1" \
    sample_reads --where 'called = 1' \
    | gzip >> ${filename}
