#!/usr/bin/env bash

# get the current time
now=`date '+%Y%m%d-%H%M'`;

sql_filename="$1-$now.sql.gz"
tar_filename="$1-$now.tgz"

# dump all the tables
mysqldump ${1} \
    | gzip > data/db/${sql_filename}

# archive all the db flags
tar czf data/db/${tar_filename} --directory data/db $1