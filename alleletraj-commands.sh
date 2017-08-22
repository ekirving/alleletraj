#!/usr/bin/env bash

# mysqldump w/ gzip
mysqldump -u root -p --single-transaction allele_trajectory \
	| gzip > allele_trajectory-fulldump.sql.gz

# import w/ gzip
gunzip < allele_trajectory-fulldump.sql.gz \
	| mysql -u root -p

# run the SNP discovery
nohup python -u build_database.py &> nohup-build.out &