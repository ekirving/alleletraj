#!/bin/bash

shopt -s nullglob

for gz in ${1}*.gz
do
	if ! gzip -t ${gz} 2> /dev/null; then 
		echo "${gz} is corrupt"
# 		gunzip -c ${gz} | head -n -1 | gzip > ${gz}.new && mv -f ${gz}.new ${gz}
	fi
done
