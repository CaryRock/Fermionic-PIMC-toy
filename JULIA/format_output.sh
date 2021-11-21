#!/bin/sh

name=$1
shift
for file in "$@"
do
    pimcave.py $file | grep $(head -n 1 "$file") | awk '{print $2}' >> "$name".dat
    pimcave.py $file | grep $(head -n 1 "$file") | awk '{print $3}' >> "$name"_err.dat
done
