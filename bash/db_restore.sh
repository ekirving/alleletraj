#!/usr/bin/env bash

# trim of the time stamp
db_name="${1:0:${#1}-14}"

sql_filename="$1.sql.gz"
tar_filename="$1.tgz"

# load all the tables
gunzip -c data/db/${sql_filename} | \
    mysql --database ${db_name}

# replace the db flags folder with the version in the archive
rm -Rf data/db/${db_name} && tar xzf data/db/${tar_filename} --directory data/db