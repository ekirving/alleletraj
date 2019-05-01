#!/usr/bin/env bash

now=`date '+%Y%m%d-%H%M'`;
filename="alleletraj_$1-$now.sql.gz"

# dump all the tables except `sample_reads`
mysqldump alleletraj_"$1" \
    | gzip > ${filename}
