#!/bin/bash

in_file=$1
split_name=$2
n_split=${3:-28}

echo "Initializing"

head -n1 ${in_file} > header.txt
tail -n +2 ${in_file}> tmp

echo "Splitting"
split -d -n l/${n_split} tmp --verbose ${split_name}

ls ${split_name}*| while read x; do cat header.txt ${x} > ${x}.tsv ;done

rm header.txt
rm tmp
rm ${split_name}??

echo "Done"
