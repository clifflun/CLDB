#!/bin/bash

in_file=$1
n_split=${2:-28}

echo "Initializing"

head -n1 ${in_file} > header.txt
tail -n +2 ${in_file}> tmp

echo "Splitting"
split -d -n l/28 tmp --verbose split

ls split*| while read x; do cat header.txt ${x} > ${x}.tsv ;done

rm header.txt
rm tmp
rm split??

echo "Done"
